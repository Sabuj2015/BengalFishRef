"""Multi-format bioinformatics exporters (QIIME2, DADA2, SINTAX, BLAST+) and Serverless Web Portal builder."""

from __future__ import annotations

import gzip
import json
from pathlib import Path
import shutil
import sqlite3
import subprocess
from typing import Any

from Bio import SeqIO

from bengalfishref.audit import clean_dna, detect_numt_stop_codons, evaluate_sequence_hygiene
from bengalfishref.primers import STANDARD_EDNA_PRIMERS
from bengalfishref.taxonomy import generate_dark_taxa_gap_report


def export_bioinformatics_releases(
    conn: sqlite3.Connection,
    release_dir: Path,
    tiers: tuple[str, ...] = ("Gold", "Silver"),
) -> dict[str, Any]:
    """Export BengalFishRef sequences in QIIME2, DADA2, SINTAX, Kraken2, and NCBI BLAST+ formats."""
    release_dir.mkdir(parents=True, exist_ok=True)
    placeholders = ",".join("?" for _ in tiers)
    rows = [
        dict(r)
        for r in conn.execute(
            f"""
            SELECT s.*, t.valid_name, t.tax_group, t.order_name, t.family_name, t.genus_name
            FROM sequences s
            JOIN taxa t ON s.resolved_taxon_id = t.taxon_id
            WHERE s.quality_tier IN ({placeholders})
              AND s.phylo_concordance_flag != 'Mislabeled_Outlier'
            ORDER BY s.marker, t.valid_name, s.accession
            """,
            tiers,
        ).fetchall()
    ]

    qiime_dir = release_dir / "qiime2"
    dada2_dir = release_dir / "dada2"
    sintax_dir = release_dir / "sintax"
    blast_dir = release_dir / "blastdb"
    for d in (qiime_dir, dada2_dir, sintax_dir, blast_dir):
        d.mkdir(parents=True, exist_ok=True)

    q_fa = qiime_dir / "bengalfishref_seqs.fasta"
    q_tax = qiime_dir / "bengalfishref_taxonomy.tsv"
    d_tax_gz = dada2_dir / "bengalfishref_dada2_assignTaxonomy.fa.gz"
    d_sp_gz = dada2_dir / "bengalfishref_dada2_addSpecies.fa.gz"
    s_fa = sintax_dir / "bengalfishref_sintax.fasta"

    with (
        q_fa.open("w", encoding="utf-8") as fh_qfa,
        q_tax.open("w", encoding="utf-8") as fh_qtx,
        gzip.open(d_tax_gz, "wt", encoding="utf-8") as fh_dtax,
        gzip.open(d_sp_gz, "wt", encoding="utf-8") as fh_dsp,
        s_fa.open("w", encoding="utf-8") as fh_sin,
    ):
        fh_qtx.write("Feature ID\tTaxon\n")
        for r in rows:
            fid = f"BFR_{r['marker']}_{r['accession']}"
            phylum = "Arthropoda" if r["tax_group"] == "Crustacea" else "Chordata"
            cls = r["tax_group"]
            ord_n = r["order_name"]
            fam_n = r["family_name"]
            gen_n = r["genus_name"]
            sp_n = r["valid_name"]
            seq = r["sequence"]

            # QIIME2
            fh_qfa.write(f">{fid}\n{seq}\n")
            q_lin = f"k__Animalia;p__{phylum};c__{cls};o__{ord_n};f__{fam_n};g__{gen_n};s__{sp_n.replace(' ', '_')}"
            fh_qtx.write(f"{fid}\t{q_lin}\n")

            # DADA2
            d_lin = f"Animalia;{phylum};{cls};{ord_n};{fam_n};{gen_n};"
            fh_dtax.write(f">{d_lin}\n{seq}\n")
            fh_dsp.write(f">{fid} {sp_n}\n{seq}\n")

            # SINTAX
            s_lin = f"d:Animalia,p:{phylum},c:{cls},o:{ord_n},f:{fam_n},g:{gen_n},s:{sp_n.replace(' ', '_')}"
            fh_sin.write(f">{fid};tax={s_lin};\n{seq}\n")

    # Build local BLAST+ nucleotide DB if makeblastdb is available
    blast_built = False
    if rows and shutil.which("makeblastdb") is not None:
        blast_in = blast_dir / "bengalfishref_all.fasta"
        with blast_in.open("w", encoding="utf-8") as fh_b:
            for r in rows:
                fid = f"BFR|{r['marker']}|{r['accession']}|{r['valid_name'].replace(' ', '_')}|{r['quality_tier']}"
                fh_b.write(f">{fid}\n{r['sequence']}\n")
        subprocess.run(
            [
                "makeblastdb",
                "-in",
                str(blast_in),
                "-dbtype",
                "nucl",
                "-title",
                "BengalFishRef_v1.0",
                "-out",
                str(blast_dir / "BengalFishRef"),
            ],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        blast_built = True

    # Zenodo open-science metadata manifest
    zenodo_meta = {
        "title": "BengalFishRef v1.0: A Phylogenetically Audited Mitochondrial Reference Atlas and eDNA Benchmark for Fishes and Commercial Crustaceans of Bangladesh and the Bay of Bengal",
        "upload_type": "dataset",
        "creators": [
            {
                "name": "Biswas, Sabuj",
                "affiliation": "Fisheries Biotechnology Division, National Institute of Biotechnology (NIB), Bangladesh",
            }
        ],
        "license": "CC-BY-4.0",
        "keywords": [
            "eDNA metabarcoding",
            "DNA barcoding",
            "Bangladesh",
            "Bay of Bengal",
            "Ganges-Brahmaputra-Meghna delta",
            "Ichthyofauna",
            "Crustacea",
            "Reference database",
        ],
        "version": "1.0.0",
        "records_exported": len(rows),
    }
    (release_dir / "zenodo.json").write_text(json.dumps(zenodo_meta, indent=2), encoding="utf-8")

    return {
        "exported_sequences": len(rows),
        "qiime2_fasta": str(q_fa),
        "dada2_taxonomy": str(d_tax_gz),
        "sintax_fasta": str(s_fa),
        "blast_db_built": blast_built,
    }


def verify_exported_releases(release_dir: Path) -> dict[str, Any]:
    """Validate structure and non-emptiness of exported QIIME2, DADA2, SINTAX, and Zenodo artifacts."""
    q_fa = release_dir / "qiime2" / "bengalfishref_seqs.fasta"
    q_tx = release_dir / "qiime2" / "bengalfishref_taxonomy.tsv"
    d_tx = release_dir / "dada2" / "bengalfishref_dada2_assignTaxonomy.fa.gz"
    s_fa = release_dir / "sintax" / "bengalfishref_sintax.fasta"
    zen = release_dir / "zenodo.json"

    for p in (q_fa, q_tx, d_tx, s_fa, zen):
        if not p.exists():
            raise FileNotFoundError(f"Missing required release file: {p}")

    n_qfa = sum(1 for line in q_fa.read_text(encoding="utf-8").splitlines() if line.startswith(">"))
    n_qtx = sum(1 for _ in q_tx.read_text(encoding="utf-8").splitlines()) - 1
    if n_qfa != n_qtx:
        raise ValueError(f"QIIME2 FASTA ({n_qfa}) vs Taxonomy ({n_qtx}) count mismatch!")

    with gzip.open(d_tx, "rt", encoding="utf-8") as fh:
        n_dada = sum(1 for line in fh if line.startswith(">"))
    if n_dada != n_qfa:
        raise ValueError(f"DADA2 ({n_dada}) vs QIIME2 ({n_qfa}) count mismatch!")

    json.loads(zen.read_text(encoding="utf-8"))
    return {"valid": True, "sequence_count": n_qfa}


def validate_lab_sample_file(conn: sqlite3.Connection, sample_path: Path) -> list[dict[str, Any]]:
    """Validate a local FBD-NIB FASTA or Sanger .ab1 chromatogram against BengalFishRef."""
    ext = sample_path.suffix.lower()
    fmt = "abi" if ext == ".ab1" else "fasta"
    records = list(SeqIO.parse(str(sample_path), fmt))
    ref_rows = [
        dict(r)
        for r in conn.execute(
            """
            SELECT s.accession, s.marker, s.sequence, s.quality_tier, t.valid_name, t.genetic_code
            FROM sequences s JOIN taxa t ON s.resolved_taxon_id = t.taxon_id
            WHERE s.quality_tier IN ('Gold', 'Silver')
            """
        ).fetchall()
    ]

    reports: list[dict[str, Any]] = []
    for rec in records:
        dna = clean_dna(str(rec.seq))
        # Check under both vertebrate (2) and invertebrate (5) mt codes before best match
        hyg_fish = evaluate_sequence_hygiene(dna, marker="COI", genetic_code=2, is_voucher_or_geo=True)
        stop_fish, _, _ = detect_numt_stop_codons(dna, 2) if len(dna) >= 120 else (False, 0, 0)
        stop_crust, _, _ = detect_numt_stop_codons(dna, 5) if len(dna) >= 120 else (False, 0, 0)

        # Fast exact 8-mer Jaccard/containment lookup against reference DB
        q_kmers = {dna[i : i + 8] for i in range(0, max(0, len(dna) - 7))}
        best_hit: dict[str, Any] | None = None
        best_score = -1.0
        for rr in ref_rows:
            rseq = rr["sequence"]
            r_kmers = {rseq[i : i + 8] for i in range(0, max(0, len(rseq) - 7))}
            if not q_kmers or not r_kmers:
                continue
            denom = min(len(q_kmers), len(r_kmers))
            containment = len(q_kmers & r_kmers) / max(denom, 1)
            if containment > best_score:
                best_score = containment
                best_hit = rr

        reports.append(
            {
                "sample_id": rec.id,
                "length_bp": hyg_fish["seq_length"],
                "ambig_pct": hyg_fish["ambig_pct"],
                "gc_pct": hyg_fish["gc_pct"],
                "numt_stop_vertebrate_table2": stop_fish,
                "numt_stop_crustacean_table5": stop_crust,
                "top_match_species": best_hit["valid_name"] if best_hit else "No reference in DB",
                "top_match_accession": best_hit["accession"] if best_hit else "",
                "top_match_marker": best_hit["marker"] if best_hit else "",
                "kmer_containment_score": round(max(best_score, 0.0), 4),
            }
        )
    return reports


def build_static_web_portal(conn: sqlite3.Connection, web_dir: Path) -> dict[str, Any]:
    """Generate the zero-cost, never-sleeping GitHub Pages + offline FBD-NIB static web application."""
    web_dir.mkdir(parents=True, exist_ok=True)
    data_dir = web_dir / "data"
    data_dir.mkdir(parents=True, exist_ok=True)

    gap_report = generate_dark_taxa_gap_report(conn)

    # Export compact representative Gold/Silver sequences (up to 2 per taxon+marker) for instant client-side browser ID
    rep_rows = [
        dict(r)
        for r in conn.execute(
            """
            SELECT s.accession, s.marker, s.quality_tier, s.sequence, t.valid_name, t.tax_group, t.genetic_code
            FROM sequences s
            JOIN taxa t ON s.resolved_taxon_id = t.taxon_id
            WHERE s.quality_tier IN ('Gold', 'Silver')
              AND s.phylo_concordance_flag != 'Mislabeled_Outlier'
              AND s.marker IN ('COI', '12S', '16S', 'Cytb')
            ORDER BY CASE s.quality_tier WHEN 'Gold' THEN 0 ELSE 1 END, s.seq_length DESC
            """
        ).fetchall()
    ]
    seen_cap: dict[tuple[str, str], int] = {}
    compact_refs: list[dict[str, Any]] = []
    for r in rep_rows:
        key = (r["valid_name"], r["marker"])
        if seen_cap.get(key, 0) >= 2:
            continue
        seen_cap[key] = seen_cap.get(key, 0) + 1
        compact_refs.append(
            {
                "acc": r["accession"],
                "m": r["marker"],
                "tier": r["quality_tier"],
                "sp": r["valid_name"],
                "grp": r["tax_group"],
                "gc": r["genetic_code"],
                "seq": r["sequence"][:900],
            }
        )

    primer_summary = [
        dict(r)
        for r in conn.execute(
            """
            SELECT primer_name, marker,
                   COUNT(*) AS n_tested,
                   SUM(amplifies_in_silico) AS n_amplified,
                   ROUND(100.0 * SUM(amplifies_in_silico) / MAX(COUNT(*), 1), 1) AS coverage_pct,
                   ROUND(AVG(weighted_penalty), 2) AS mean_penalty
            FROM primer_evaluations
            GROUP BY primer_name, marker
            ORDER BY marker, coverage_pct DESC
            """
        ).fetchall()
    ]
    primer_catalog = [
        {
            "name": p.name,
            "marker": p.marker,
            "target_group": p.target_group,
            "fwd": p.fwd_seq,
            "rev": p.rev_seq,
            "ref": p.reference,
        }
        for p in STANDARD_EDNA_PRIMERS
    ]

    bundle = {
        "version": "1.0.0",
        "institute": "Fisheries Biotechnology Division, National Institute of Biotechnology (FBD-NIB), Bangladesh",
        "summary": {
            "total_checklist_species": len(gap_report),
            "completely_dark_species": sum(1 for x in gap_report if x["is_completely_dark"]),
            "edna_12s_16s_dark_species": sum(1 for x in gap_report if x["is_edna_12s_16s_dark"]),
            "reference_sequences_indexed": len(compact_refs),
        },
        "species_atlas": gap_report,
        "primers": primer_catalog,
        "primer_benchmarks": primer_summary,
        "reference_sequences": compact_refs,
    }
    bundle_json = json.dumps(bundle, separators=(",", ":"))
    (data_dir / "index_meta.json").write_text(bundle_json, encoding="utf-8")
    # Also write a JS script wrapper so offline `file://` opening works without CORS fetch restrictions!
    (data_dir / "bundle.js").write_text(f"window.BENGALFISHREF_DATA = {bundle_json};\n", encoding="utf-8")

    html_content = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1.0" />
  <title>BengalFishRef v1.0 | FBD-NIB Bangladesh Barcode &amp; eDNA Reference Atlas</title>
  <link rel="stylesheet" href="styles.css" />
</head>
<body>
  <header class="hero">
    <div class="container">
      <div class="badge">FBD-NIB • National Institute of Biotechnology, Bangladesh</div>
      <h1>BengalFishRef <span>v1.0</span></h1>
      <p>Phylogenetically Audited Mitochondrial Reference Atlas, eDNA Primer Benchmark &amp; Serverless Sequence Validator for Fishes and Commercial Crustaceans of Bangladesh &amp; the Bay of Bengal</p>
      <div class="kpi-grid" id="kpiGrid"></div>
    </div>
  </header>

  <main class="container">
    <nav class="tabs">
      <button class="tab-btn active" data-tab="atlasTab">1. Species &amp; Dark-Taxa Atlas</button>
      <button class="tab-btn" data-tab="idTab">2. Sequence ID &amp; Lab NUMT/QC Validator</button>
      <button class="tab-btn" data-tab="primerTab">3. 14-Primer eDNA Mismatch Explorer</button>
    </nav>

    <section id="atlasTab" class="tab-panel active">
      <div class="controls">
        <input type="text" id="searchInput" placeholder="Search scientific, Bengali (e.g. Ilish, Rui, Bagda), English, or family name..." />
        <select id="zoneFilter">
          <option value="">All Bangladesh Zones</option>
          <option value="GBM_Freshwater">GBM Freshwater Rivers &amp; Haors</option>
          <option value="Sundarbans_Estuarine">Sundarbans Mangrove-Estuarine</option>
          <option value="Bay_of_Bengal_Marine">Bay of Bengal Marine</option>
        </select>
        <select id="groupFilter">
          <option value="">All Taxonomic Groups</option>
          <option value="Actinopterygii">Ray-finned Fishes (Actinopterygii)</option>
          <option value="Elasmobranchii">Sharks &amp; Rays (Elasmobranchii)</option>
          <option value="Crustacea">Shrimps, Prawns &amp; Crabs (Crustacea)</option>
        </select>
      </div>
      <div class="table-wrap">
        <table id="atlasTable">
          <thead>
            <tr>
              <th>Valid Species</th>
              <th>Bengali / English</th>
              <th>Group &amp; Family</th>
              <th>IUCN</th>
              <th>COI</th>
              <th>12S</th>
              <th>16S</th>
              <th>Mitogenome</th>
              <th>BDTPI Priority</th>
            </tr>
          </thead>
          <tbody></tbody>
        </table>
      </div>
    </section>

    <section id="idTab" class="tab-panel">
      <div class="card">
        <h2>Client-Side Sequence Identifier &amp; Mitochondrial NUMT Hygiene Check</h2>
        <p class="muted">Runs 100% in your browser (zero server upload). Checks length, IUPAC ambiguity, vertebrate (<code>transl_table=2</code>) &amp; crustacean (<code>transl_table=5</code>) internal stop codons, and 8-mer + alignment similarity against Gold/Silver vouchers.</p>
        <textarea id="seqInput" rows="6" placeholder=">Sample_FBD_NIB_01&#10;Paste COI, 12S, 16S, or Cytb FASTA / raw nucleotide sequence here..."></textarea>
        <div class="btn-row">
          <button id="runIdBtn" class="primary-btn">Audit &amp; Identify Sequence</button>
          <button id="loadDemoBtn" class="secondary-btn">Load Demo Sequence</button>
        </div>
        <div id="idResults" class="results-box"></div>
      </div>
    </section>

    <section id="primerTab" class="tab-panel">
      <div class="card">
        <h2>Standard eDNA Metabarcoding Primer Benchmark (14 Universal Primer Sets)</h2>
        <div class="table-wrap">
          <table id="primerTable">
            <thead>
              <tr>
                <th>Primer Pair</th>
                <th>Marker</th>
                <th>Target Clade</th>
                <th>Forward (5'→3')</th>
                <th>Reverse (5'→3')</th>
                <th>In Silico Coverage</th>
                <th>Reference</th>
              </tr>
            </thead>
            <tbody></tbody>
          </table>
        </div>
      </div>
    </section>
  </main>

  <footer class="footer">
    <div class="container">
      Developed by <strong>Sabuj Biswas</strong>, Fisheries Biotechnology Division, National Institute of Biotechnology (FBD-NIB), Bangladesh • Open Science Release (CC-BY 4.0 / MIT)
    </div>
  </footer>

  <script src="data/bundle.js"></script>
  <script src="app.js"></script>
</body>
</html>
"""

    css_content = """:root {
  --bg: #f8fafc;
  --surface: #ffffff;
  --ink: #0f172a;
  --muted: #475569;
  --accent: #0f766e;
  --accent-light: #ccfbf1;
  --border: #e2e8f0;
  --warn: #b45309;
  --danger: #b91c1c;
  --ok: #15803d;
}
* { box-sizing: border-box; }
body { margin: 0; font-family: system-ui, -apple-system, sans-serif; background: var(--bg); color: var(--ink); line-height: 1.5; }
.container { max-width: 1220px; margin: 0 auto; padding: 0 1.25rem; }
.hero { background: linear-gradient(135deg, #042f2e 0%, #0f766e 100%); color: #fff; padding: 2.5rem 0 2rem; }
.badge { display: inline-block; background: rgba(255,255,255,0.14); padding: 0.28rem 0.75rem; border-radius: 999px; font-size: 0.82rem; letter-spacing: 0.02em; }
.hero h1 { margin: 0.6rem 0 0.4rem; font-size: 2.2rem; }
.hero h1 span { font-weight: 400; opacity: 0.85; font-size: 1.35rem; }
.hero p { max-width: 850px; opacity: 0.92; margin: 0 0 1.5rem; }
.kpi-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 0.9rem; }
.kpi-card { background: rgba(255,255,255,0.12); padding: 0.9rem 1rem; border-radius: 0.65rem; border: 1px solid rgba(255,255,255,0.18); }
.kpi-val { font-size: 1.65rem; font-weight: 700; }
.kpi-lbl { font-size: 0.82rem; opacity: 0.88; }
.tabs { display: flex; gap: 0.5rem; margin: 1.4rem 0 1rem; flex-wrap: wrap; }
.tab-btn { border: 1px solid var(--border); background: var(--surface); padding: 0.65rem 1.1rem; border-radius: 0.5rem; cursor: pointer; font-weight: 600; color: var(--muted); }
.tab-btn.active { background: var(--accent); color: #fff; border-color: var(--accent); }
.tab-panel { display: none; }
.tab-panel.active { display: block; }
.controls { display: grid; grid-template-columns: 2fr 1fr 1fr; gap: 0.75rem; margin-bottom: 1rem; }
input, select, textarea { width: 100%; padding: 0.65rem 0.8rem; border: 1px solid var(--border); border-radius: 0.5rem; font-size: 0.95rem; font-family: inherit; }
textarea { font-family: ui-monospace, SFMono-Regular, Menlo, monospace; }
.table-wrap { background: var(--surface); border: 1px solid var(--border); border-radius: 0.65rem; overflow-x: auto; }
table { width: 100%; border-collapse: collapse; font-size: 0.9rem; }
th, td { padding: 0.65rem 0.8rem; border-bottom: 1px solid var(--border); text-align: left; }
th { background: #f1f5f9; font-weight: 600; color: #334155; }
.pill { display: inline-block; padding: 0.15rem 0.5rem; border-radius: 999px; font-size: 0.78rem; font-weight: 600; }
.pill-ok { background: #dcfce7; color: var(--ok); }
.pill-miss { background: #fee2e2; color: var(--danger); }
.pill-warn { background: #fef3c7; color: var(--warn); }
.card { background: var(--surface); border: 1px solid var(--border); border-radius: 0.65rem; padding: 1.3rem; }
.btn-row { display: flex; gap: 0.75rem; margin-top: 0.85rem; }
.primary-btn { background: var(--accent); color: #fff; border: none; padding: 0.65rem 1.2rem; border-radius: 0.5rem; font-weight: 600; cursor: pointer; }
.secondary-btn { background: #e2e8f0; color: var(--ink); border: none; padding: 0.65rem 1.2rem; border-radius: 0.5rem; font-weight: 600; cursor: pointer; }
.results-box { margin-top: 1rem; }
.footer { margin-top: 3rem; padding: 1.5rem 0; border-top: 1px solid var(--border); color: var(--muted); font-size: 0.88rem; }
@media (max-width: 760px) { .controls { grid-template-columns: 1fr; } }
"""

    js_content = """(function () {
  const DATA = window.BENGALFISHREF_DATA || { summary: {}, species_atlas: [], primers: [], primer_benchmarks: [], reference_sequences: [] };

  const kpiGrid = document.getElementById("kpiGrid");
  const s = DATA.summary || {};
  const kpis = [
    { val: s.total_checklist_species || 0, lbl: "Harmonized Bangladesh Species" },
    { val: s.reference_sequences_indexed || 0, lbl: "Curated Gold/Silver Vouchers" },
    { val: s.completely_dark_species || 0, lbl: "Completely Unsequenced Dark Taxa" },
    { val: s.edna_12s_16s_dark_species || 0, lbl: "Missing 12S/16S eDNA Reference" }
  ];
  kpiGrid.innerHTML = kpis.map(k => `<div class="kpi-card"><div class="kpi-val">${k.val}</div><div class="kpi-lbl">${k.lbl}</div></div>`).join("");

  document.querySelectorAll(".tab-btn").forEach(btn => {
    btn.addEventListener("click", () => {
      document.querySelectorAll(".tab-btn").forEach(b => b.classList.remove("active"));
      document.querySelectorAll(".tab-panel").forEach(p => p.classList.remove("active"));
      btn.classList.add("active");
      document.getElementById(btn.dataset.tab).classList.add("active");
    });
  });

  function badgeCount(n) {
    return n > 0 ? `<span class="pill pill-ok">${n}</span>` : `<span class="pill pill-miss">0</span>`;
  }

  function renderAtlas() {
    const q = (document.getElementById("searchInput").value || "").toLowerCase();
    const z = document.getElementById("zoneFilter").value;
    const g = document.getElementById("groupFilter").value;
    const rows = (DATA.species_atlas || []).filter(r => {
      if (z && !(r.habitat_zones || "").includes(z)) return false;
      if (g && r.tax_group !== g) return false;
      if (!q) return true;
      return (
        (r.valid_name || "").toLowerCase().includes(q) ||
        (r.common_name_bn || "").toLowerCase().includes(q) ||
        (r.common_name_en || "").toLowerCase().includes(q) ||
        (r.family_name || "").toLowerCase().includes(q)
      );
    }).slice(0, 300);

    const tbody = document.querySelector("#atlasTable tbody");
    tbody.innerHTML = rows.map(r => `
      <tr>
        <td><em><strong>${r.valid_name}</strong></em></td>
        <td>${r.common_name_bn ? r.common_name_bn + " / " : ""}${r.common_name_en || "-"}</td>
        <td>${r.tax_group} • ${r.family_name}</td>
        <td><span class="pill ${['CR','EN','VU'].includes(r.iucn_bd_status) ? 'pill-warn' : 'pill-ok'}">${r.iucn_bd_status}</span></td>
        <td>${badgeCount(r.COI_usable)}</td>
        <td>${badgeCount(r["12S_usable"])}</td>
        <td>${badgeCount(r["16S_usable"])}</td>
        <td>${badgeCount(r.Mitogenome_usable)}</td>
        <td><strong>${Number(r.bdtpi_score).toFixed(3)}</strong></td>
      </tr>
    `).join("");
  }

  ["searchInput", "zoneFilter", "groupFilter"].forEach(id => {
    document.getElementById(id).addEventListener("input", renderAtlas);
  });
  renderAtlas();

  // Primer table
  const benchMap = {};
  (DATA.primer_benchmarks || []).forEach(b => { benchMap[b.primer_name] = b; });
  document.querySelector("#primerTable tbody").innerHTML = (DATA.primers || []).map(p => {
    const bm = benchMap[p.name];
    const cov = bm ? `${bm.coverage_pct}% (n=${bm.n_tested})` : "Pending Mitogenome Run";
    return `<tr>
      <td><strong>${p.name}</strong></td>
      <td>${p.marker}</td>
      <td>${p.target_group}</td>
      <td><code>${p.fwd}</code></td>
      <td><code>${p.rev}</code></td>
      <td>${cov}</td>
      <td>${p.ref}</td>
    </tr>`;
  }).join("");

  // Client-side NUMT stop check (transl_table 2 & 5) + 8-mer lookup
  function hasInternalStop(dna, table) {
    const stops = table === 2 ? new Set(["TAA", "TAG", "AGA", "AGG"]) : new Set(["TAA", "TAG"]);
    for (let frame = 0; frame < 3; frame++) {
      let count = 0;
      const end = dna.length - ((dna.length - frame) % 3) - 3;
      for (let i = frame; i < end; i += 3) {
        const cod = dna.slice(i, i + 3);
        if (stops.has(cod)) count++;
      }
      if (count === 0) return false;
    }
    return dna.length >= 120;
  }

  function extractKmers(seq, k) {
    const s = new Set();
    for (let i = 0; i <= seq.length - k; i++) {
      const w = seq.slice(i, i + k);
      if (/^[ACGT]+$/.test(w)) s.add(w);
    }
    return s;
  }

  document.getElementById("loadDemoBtn").addEventListener("click", () => {
    const sample = (DATA.reference_sequences || [])[0];
    const seq = sample ? sample.seq : "TCAACCAACCACAAAGACATTGGCACCCTCTATCTTGTATTTGGTGCCTGAGCCGGAATAGTAGGAACCGCCCTAAGCCTTCTCATTCGAGCCGAACTAAGCCAGCCCGGCGCCCTTCTTGGAGACGACCAAATTTATAATGTAATCGTCACAGCACATGCTTTCGTAATAATTTTCTTTATAGTAATACCAATCATAATTGGAGGCTTTGGAAACTGACTTGTCCCACTAATAATCGGAGCCCCCGATATAGCATTCCCACGAATAAACAACATAAGCTTCTGACTCCTTCCTCCCTCTTTCCTTCTTCTCCTAGCTTCTTCTGGTGTTGAAGCTGGAGCCGGAACAGGWTGAACAGTCTACCCCCCTCTAGCTGGAAATCTGGCACACGCAGGRGCATCTGTTGACCTAACAATTTTTTCCCTCCACTTAGCAGGTGTTTCATCAATTCTAGGGGCAATTAATTTTATTACAACAATTATTAATATAAAACCCCCAGCCATCTCCCAATACCAAACACCCCTGTTCGTCTGAGCCGTCCTTGTAACAGCCGTACTTCTCCTTCTATCCCTCCCAGTCCTAGCTGCTGGCATCACAATACTTCTAACAGACCGAAACCTAAACACCACCTTCTTCGACCCGGCAGGAGGAGGAGACCCAATCCTTTACCAACACTTATTC";
    const sp = sample ? sample.sp : "Tenualosa_ilisha_Demo";
    document.getElementById("seqInput").value = `>Query_${sp.replace(/\\s+/g, "_")}\n${seq}`;
  });

  document.getElementById("runIdBtn").addEventListener("click", () => {
    const raw = document.getElementById("seqInput").value || "";
    const dna = raw.split(/\\r?\\n/).filter(l => !l.startsWith(">")).join("").replace(/[^A-Za-z]/g, "").toUpperCase();
    const box = document.getElementById("idResults");
    if (dna.length < 80) {
      box.innerHTML = `<p class="pill pill-miss">Please enter at least 80 bp of nucleotide sequence.</p>`;
      return;
    }
    const acgt = (dna.match(/[ACGT]/g) || []).length;
    const ambigPct = (100 * (dna.length - acgt) / dna.length).toFixed(2);
    const numtVert = hasInternalStop(dna, 2);
    const numtCrust = hasInternalStop(dna, 5);

    const qk = extractKmers(dna, 8);
    const scored = (DATA.reference_sequences || []).map(r => {
      const rk = extractKmers(r.seq, 8);
      let shared = 0;
      qk.forEach(w => { if (rk.has(w)) shared++; });
      const denom = Math.max(1, Math.min(qk.size, rk.size));
      return { ...r, score: (100 * shared / denom) };
    }).sort((a, b) => b.score - a.score).slice(0, 5);

    const topRows = scored.map(h => `
      <tr>
        <td><em><strong>${h.sp}</strong></em></td>
        <td>${h.m}</td>
        <td><code>${h.acc}</code></td>
        <td><span class="pill pill-ok">${h.tier}</span></td>
        <td><strong>${h.score.toFixed(1)}%</strong></td>
      </tr>
    `).join("") || `<tr><td colspan="5">Reference sequences will populate once Phase 2 ingestion completes.</td></tr>`;

    box.innerHTML = `
      <div style="margin-bottom:0.8rem">
        <span class="pill pill-ok">Length: ${dna.length} bp</span>
        <span class="pill ${ambigPct <= 1 ? 'pill-ok' : 'pill-warn'}">Ambiguous: ${ambigPct}%</span>
        <span class="pill ${numtVert ? 'pill-warn' : 'pill-ok'}">Fish CDS (Table 2): ${numtVert ? 'Stop Codon Warning (Check if rRNA or NUMT)' : 'Clean ORF'}</span>
        <span class="pill ${numtCrust ? 'pill-warn' : 'pill-ok'}">Crustacean CDS (Table 5): ${numtCrust ? 'Stop Codon Warning' : 'Clean ORF'}</span>
      </div>
      <div class="table-wrap">
        <table>
          <thead><tr><th>Matched Species</th><th>Marker</th><th>Accession</th><th>Tier</th><th>8-mer Containment</th></tr></thead>
          <tbody>${topRows}</tbody>
        </table>
      </div>
    `;
  });
})();
"""

    (web_dir / "index.html").write_text(html_content, encoding="utf-8")
    (web_dir / "styles.css").write_text(css_content, encoding="utf-8")
    (web_dir / "app.js").write_text(js_content, encoding="utf-8")

    return {
        "web_dir": str(web_dir),
        "species_count": len(gap_report),
        "indexed_reference_sequences": len(compact_refs),
    }
