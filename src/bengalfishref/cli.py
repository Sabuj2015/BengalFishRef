"""Command-line interface for BengalFishRef v1.0."""

from __future__ import annotations

import argparse
from collections import defaultdict
import json
from pathlib import Path
import random
import tempfile
from typing import Any

from bengalfishref.audit import run_mafft_and_k2p_audit
from bengalfishref.benchmark import evaluate_ecosystem_benchmark, split_reference_and_queries
from bengalfishref.export import (
    build_static_web_portal,
    export_bioinformatics_releases,
    validate_lab_sample_file,
    verify_exported_releases,
)
from bengalfishref.ingest import (
    ingest_bold_coi_by_genera,
    ingest_ncbi_for_checklist,
    insert_ingested_sequence,
)
from bengalfishref.primers import (
    STANDARD_EDNA_PRIMERS,
    compute_amplicon_species_resolution,
    design_bengal_degenerate_primer,
    evaluate_primer_pair_on_template,
)
from bengalfishref.schema import connect_db
from bengalfishref.taxonomy import build_master_checklist, generate_dark_taxa_gap_report

WORKSPACE_ROOT = Path("/home/sabuj/Research/01_Active_Projects/Antigravity")
DEFAULT_DB = WORKSPACE_ROOT / "02_Data" / "processed" / "bengalfishref.sqlite"
DEFAULT_CHECKLIST_TSV = WORKSPACE_ROOT / "02_Data" / "reference" / "bangladesh_master_checklist.tsv"
DEFAULT_RELEASE_DIR = WORKSPACE_ROOT / "02_Data" / "processed" / "releases" / "v1.0"
DEFAULT_WEB_DIR = WORKSPACE_ROOT / "03_Pipelines" / "bengalfishref" / "web"


def cmd_build_checklist(args: argparse.Namespace) -> None:
    conn = connect_db(args.db)
    summary = build_master_checklist(
        conn,
        include_remote_reconciliation=not args.offline_only,
        export_tsv=Path(args.export_tsv),
    )
    print(json.dumps(summary, indent=2))


def cmd_ingest_sequences(args: argparse.Namespace) -> None:
    conn = connect_db(args.db)
    markers = tuple(m.strip() for m in args.markers.split(",") if m.strip())
    ncbi_stats = ingest_ncbi_for_checklist(
        conn,
        markers=markers,
        species_batch_size=args.batch_size,
        retmax_per_batch=args.retmax,
        max_batches=args.max_batches,
    )
    bold_stats = {"skipped": True}
    if args.include_bold:
        bold_stats = ingest_bold_coi_by_genera(conn, max_genera=args.max_bold_genera)
    print(json.dumps({"ncbi": ncbi_stats, "bold": bold_stats}, indent=2))


def run_db_phylogenetic_audit(conn: Any, threads: int = 4, run_iqtree: bool = False) -> dict[str, int]:
    """Run family+marker MAFFT + K2P (+ optional IQ-TREE2) audit across stored non-quarantined sequences."""
    rows = [
        dict(r)
        for r in conn.execute(
            """
            SELECT s.seq_uid, s.resolved_taxon_id, s.marker, s.sequence, s.quality_tier,
                   s.quarantine_reason, t.valid_name, t.family_name, t.genus_name
            FROM sequences s
            JOIN taxa t ON s.resolved_taxon_id = t.taxon_id
            WHERE s.quality_tier IN ('Gold', 'Silver')
              AND s.marker IN ('COI', '12S', '16S', 'Cytb')
            """
        ).fetchall()
    ]
    groups: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for r in rows:
        fam = r["family_name"]
        grp_key = f"GENUS_{r['genus_name']}" if fam in ("Actinopterygii", "Elasmobranchii", "Incertae_sedis", "Crustacea") else fam
        groups[(grp_key, r["marker"])].append(r)

    stats = {"families_audited": 0, "sequences_checked": 0, "mislabeled_quarantined": 0, "cryptic_flagged": 0}
    for (g_name, m_name), grp_recs in groups.items():
        if len(grp_recs) < 2:
            continue
        audit_res, _ = run_mafft_and_k2p_audit(grp_recs, threads=threads, run_iqtree=run_iqtree)
        with conn:
            for uid, info in audit_res.items():
                stats["sequences_checked"] += 1
                flag = info["phylo_concordance_flag"]
                if flag == "Mislabeled_Outlier":
                    stats["mislabeled_quarantined"] += 1
                    conn.execute(
                        """
                        UPDATE sequences
                        SET k2p_nearest_conspecific = ?,
                            k2p_nearest_heterospecific = ?,
                            nearest_heterospecific_taxon = ?,
                            phylo_concordance_flag = ?,
                            quality_tier = 'Quarantined',
                            quarantine_reason = CASE
                                WHEN quarantine_reason = '' THEN ?
                                ELSE quarantine_reason || ';' || ?
                            END
                        WHERE seq_uid = ?
                        """,
                        (
                            info["k2p_nearest_conspecific"],
                            info["k2p_nearest_heterospecific"],
                            info["nearest_heterospecific_taxon"],
                            flag,
                            info["quarantine_add"],
                            info["quarantine_add"],
                            uid,
                        ),
                    )
                else:
                    if flag == "Cryptic_Lineage_Candidate":
                        stats["cryptic_flagged"] += 1
                    conn.execute(
                        """
                        UPDATE sequences
                        SET k2p_nearest_conspecific = ?,
                            k2p_nearest_heterospecific = ?,
                            nearest_heterospecific_taxon = ?,
                            phylo_concordance_flag = ?
                        WHERE seq_uid = ?
                        """,
                        (
                            info["k2p_nearest_conspecific"],
                            info["k2p_nearest_heterospecific"],
                            info["nearest_heterospecific_taxon"],
                            flag,
                            uid,
                        ),
                    )
        stats["families_audited"] += 1
    return stats


def run_db_primer_benchmark(conn: Any) -> dict[str, Any]:
    """Evaluate all 14 standard eDNA primer sets against full mitogenomes and flanking locus sequences."""
    rows = [
        dict(r)
        for r in conn.execute(
            """
            SELECT s.seq_uid, s.marker, s.sequence, t.taxon_id, t.valid_name, t.genus_name, t.tax_group
            FROM sequences s
            JOIN taxa t ON s.resolved_taxon_id = t.taxon_id
            WHERE s.quality_tier IN ('Gold', 'Silver')
            """
        ).fetchall()
    ]
    total_evals = 0
    resolution_summary: dict[str, Any] = {}
    degenerate_designs: dict[str, Any] = {}

    for pair in STANDARD_EDNA_PRIMERS:
        cand_rows = [
            r
            for r in rows
            if r["marker"] in ("Mitogenome", pair.marker)
            and (
                (pair.target_group == "Crustacea" and r["tax_group"] == "Crustacea")
                or (pair.target_group == "Elasmobranchii" and r["tax_group"] == "Elasmobranchii")
                or (pair.target_group in ("Teleostei", "Actinopterygii") and r["tax_group"] == "Actinopterygii")
                or (pair.target_group in ("Pisces", "Vertebrata", "Marine_Vertebrata") and r["tax_group"] != "Crustacea")
                or (pair.target_group == "Metazoa")
            )
        ]
        amplicons: list[tuple[int, str, str, str]] = []
        fwd_windows: list[str] = []

        with conn:
            for r in cand_rows:
                ev = evaluate_primer_pair_on_template(pair, r["sequence"])
                if ev is None:
                    continue
                total_evals += 1
                fwd_windows.append(ev["fwd_target_window"])
                if ev["amplifies_in_silico"] and ev["amplicon_seq"]:
                    amplicons.append((int(r["taxon_id"]), r["valid_name"], r["genus_name"], ev["amplicon_seq"]))
                conn.execute(
                    """
                    INSERT OR REPLACE INTO primer_evaluations (
                        primer_name, marker, taxon_id, seq_uid,
                        fwd_mismatches, rev_mismatches, fwd_3prime_mm, rev_3prime_mm,
                        weighted_penalty, delta_tm_est, amplifies_in_silico,
                        amplicon_length, amplicon_seq
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                    """,
                    (
                        pair.name,
                        pair.marker,
                        int(r["taxon_id"]),
                        r["seq_uid"],
                        ev["fwd_mismatches"],
                        ev["rev_mismatches"],
                        ev["fwd_3prime_mm"],
                        ev["rev_3prime_mm"],
                        ev["weighted_penalty"],
                        ev["delta_tm_est"],
                        ev["amplifies_in_silico"],
                        ev["amplicon_length"],
                        ev["amplicon_seq"],
                    ),
                )

        if amplicons:
            resolution_summary[pair.name] = compute_amplicon_species_resolution(amplicons)
        if fwd_windows and pair.name in ("MiFish-U", "Crust16S", "Teleo"):
            degenerate_designs[f"{pair.name}-Bengal"] = design_bengal_degenerate_primer(pair.fwd_seq, fwd_windows)

    return {
        "evaluations_stored": total_evals,
        "amplicon_resolution": resolution_summary,
        "bengal_degenerate_primers": degenerate_designs,
    }


def run_db_classifier_benchmark(conn: Any, run_id: str = "v1.0_bench", threads: int = 4) -> list[dict[str, Any]]:
    """Compare Unfiltered Raw DB vs. Curated BengalFishRef Gold/Silver across 4 Bangladesh ecosystems."""
    all_rows = [
        dict(r)
        for r in conn.execute(
            """
            SELECT s.*, t.valid_name, t.habitat_zones, t.aquaculture_tier
            FROM sequences s
            JOIN taxa t ON s.resolved_taxon_id = t.taxon_id
            WHERE s.marker = 'COI' AND s.seq_length >= 480
            """
        ).fetchall()
    ]
    stored: list[dict[str, Any]] = []
    ecosystems = ["GBM_Freshwater", "Sundarbans_Estuarine", "Bay_of_Bengal_Marine", "Aquaculture_Core"]

    for eco in ecosystems:
        if eco == "Aquaculture_Core":
            eco_rows = [r for r in all_rows if r["aquaculture_tier"] == "Core_Commercial"]
        else:
            eco_rows = [r for r in all_rows if eco in (r["habitat_zones"] or "")]

        if len(eco_rows) < 8:
            continue

        # Build clean ground-truth records from verified non-quarantined sequences
        clean_pool = [
            {
                "seq_uid": r["seq_uid"],
                "true_species": r["valid_name"],
                "reported_label": r["valid_name"],
                "raw_genbank_label": r["reported_taxon"] or r["valid_name"],
                "sequence": r["sequence"],
                "quality_tier": r["quality_tier"],
            }
            for r in eco_rows
            if r["quality_tier"] in ("Gold", "Silver") and r["phylo_concordance_flag"] != "Mislabeled_Outlier"
        ]
        if len({x["true_species"] for x in clean_pool}) < 4:
            continue

        ref_curated, query_templates, _ = split_reference_and_queries(clean_pool, holdout_dark_taxa_frac=0.12)
        if len(query_templates) > 250:
            rng_cap = random.Random(20261117)
            query_templates = rng_cap.sample(query_templates, 250)

        # Raw unfiltered condition keeps exact unharmonized GenBank organism labels (e.g., Hilsa ilisha,
        # Catla catla, Puntius ticto) PLUS quarantined NUMT and mislabeled GenBank accessions!
        quarantined_rows = [
            {
                "seq_uid": r["seq_uid"],
                "true_species": r["valid_name"],
                "reported_label": r["reported_taxon"] or r["valid_name"],
                "sequence": r["sequence"],
                "quality_tier": "Quarantined",
            }
            for r in eco_rows
            if r["quality_tier"] == "Quarantined"
        ]
        ref_unfiltered = [
            {**x, "reported_label": x["raw_genbank_label"]} for x in ref_curated
        ] + quarantined_rows

        for err_rate in (0.0, 0.005, 0.01):
            for cond_name, ref_set in (
                ("Raw_Unfiltered_Public", ref_unfiltered),
                ("BengalFishRef_Curated", ref_curated),
            ):
                metrics_list = evaluate_ecosystem_benchmark(
                    ref_set,
                    query_templates,
                    reads_per_template=2,
                    error_rate=err_rate,
                    threads=threads,
                )
                with conn:
                    for m in metrics_list:
                        rec = {
                            "run_id": run_id,
                            "ecosystem": eco,
                            "marker": "COI",
                            "error_rate": err_rate,
                            "db_condition": cond_name,
                            **m,
                        }
                        stored.append(rec)
                        conn.execute(
                            """
                            INSERT OR REPLACE INTO benchmark_results (
                                run_id, ecosystem, marker, error_rate, db_condition, classifier,
                                n_queries, precision_macro, recall_macro, f1_macro,
                                f1_ci_low, f1_ci_high, false_positive_exotic_rate, latency_ms_per_query
                            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                            """,
                            (
                                run_id,
                                eco,
                                "COI",
                                err_rate,
                                cond_name,
                                m["classifier"],
                                m["n_queries"],
                                m["precision_macro"],
                                m["recall_macro"],
                                m["f1_macro"],
                                m["f1_ci_low"],
                                m["f1_ci_high"],
                                m["false_positive_exotic_rate"],
                                m["latency_ms_per_query"],
                            ),
                        )
    return stored


def cmd_smoke_test(args: argparse.Namespace) -> None:
    """Deterministic offline 6-species end-to-end verification exercising every pipeline component."""
    rng = random.Random(20261116)
    with tempfile.TemporaryDirectory(prefix="bfr_smoke_") as tmpdir:
        tmp = Path(tmpdir)
        db_path = tmp / "smoke.sqlite"
        conn = connect_db(db_path)
        build_master_checklist(conn, include_remote_reconciliation=False, export_tsv=tmp / "checklist.tsv")

        target_names = [
            "Tenualosa ilisha",
            "Labeo rohita",
            "Labeo catla",
            "Oreochromis niloticus",
            "Penaeus monodon",
            "Macrobrachium rosenbergii",
        ]
        taxa_map = {
            r["valid_name"]: dict(r)
            for r in conn.execute(
                f"SELECT * FROM taxa WHERE valid_name IN ({','.join('?' for _ in target_names)})",
                target_names,
            ).fetchall()
        }

        # Generate realistic COI templates (~10% interspecific divergence on a shared ancestral ORF backbone)
        safe_codons = ["GCT", "GCC", "TTC", "CTG", "ATC", "GTC", "TCC", "CCA", "ACT", "TAC", "CAC", "AAC", "GAC", "GAA", "TGC"]
        ancestral = [rng.choice(safe_codons) for _ in range(180)]
        species_templates: dict[str, str] = {}
        for idx_sp, sp in enumerate(target_names):
            codons = list(ancestral)
            for pos in range(idx_sp * 22, idx_sp * 22 + 22):
                codons[pos % 180] = safe_codons[(idx_sp + pos) % len(safe_codons)]
            species_templates[sp] = "".join(codons)

        with conn:
            for sp in target_names:
                base_seq = species_templates[sp]
                sp_tag = sp.replace(" ", "_").upper()
                for rep in range(3):
                    s_list = list(base_seq)
                    s_list[rep * 9 : rep * 9 + 3] = list("GCC")
                    insert_ingested_sequence(
                        conn,
                        source_db="NCBI",
                        accession=f"SMOKE_{sp_tag}_{rep+1}",
                        reported_taxon=sp,
                        taxon_meta=taxa_map[sp],
                        marker="COI",
                        raw_seq="".join(s_list),
                        country_origin="Bangladesh",
                    )

            # Inject 1 genuine NUMT with internal stops across all 6 reading frames -> must be quarantined!
            stop6_unit = "TAACTAGATAATAATTATTATCTAGTTA"
            numt_seq = stop6_unit * 20
            insert_ingested_sequence(
                conn,
                source_db="NCBI",
                accession="SMOKE_NUMT_01",
                reported_taxon="Tenualosa ilisha",
                taxon_meta=taxa_map["Tenualosa ilisha"],
                marker="COI",
                raw_seq=numt_seq,
                country_origin="Global",
            )

            # Inject 1 mislabeled GenBank record (Labeo rohita sequence mislabeled as Labeo catla!)
            insert_ingested_sequence(
                conn,
                source_db="NCBI",
                accession="SMOKE_MISLABEL_01",
                reported_taxon="Labeo catla",
                taxon_meta=taxa_map["Labeo catla"],
                marker="COI",
                raw_seq=species_templates["Labeo rohita"],
                country_origin="Global",
            )

            # Inject synthetic 12S template flanked by exact MiFish-U forward + reverse-complement binding sites
            mifish_fwd = "GTCGGTAAAACTCGTGCCAGC"
            # Reverse primer CATAGTGGGGTATCTAATCCCAGTTTG binds sense strand as reverse complement:
            mifish_rev_rc = "CAAACTGGGATTAGATACCCCACTATG"
            amp_body = "".join(rng.choice("ACGT") for _ in range(170))
            insert_ingested_sequence(
                conn,
                source_db="NCBI",
                accession="SMOKE_12S_ILISHA",
                reported_taxon="Tenualosa ilisha",
                taxon_meta=taxa_map["Tenualosa ilisha"],
                marker="12S",
                raw_seq="AAAA" + mifish_fwd + amp_body + mifish_rev_rc + "TTTT",
                country_origin="Bangladesh",
            )

        audit_stats = run_db_phylogenetic_audit(conn, threads=args.threads, run_iqtree=True)
        primer_stats = run_db_primer_benchmark(conn)
        bench_stats = run_db_classifier_benchmark(conn, run_id="smoke", threads=args.threads)
        rel_dir = tmp / "release"
        exp_stats = export_bioinformatics_releases(conn, rel_dir)
        ver_stats = verify_exported_releases(rel_dir)
        web_stats = build_static_web_portal(conn, tmp / "web")
        gap_count = len(generate_dark_taxa_gap_report(conn))

        # Assert critical scientific invariants
        numt_row = conn.execute("SELECT quality_tier FROM sequences WHERE accession='SMOKE_NUMT_01'").fetchone()
        mislbl_row = conn.execute("SELECT quality_tier, phylo_concordance_flag FROM sequences WHERE accession='SMOKE_MISLABEL_01'").fetchone()
        assert numt_row["quality_tier"] == "Quarantined", "NUMT was not quarantined!"
        assert mislbl_row["quality_tier"] == "Quarantined", "Mislabeled Labeo sequence was not quarantined by MAFFT+K2P!"

        print(
            json.dumps(
                {
                    "smoke_test": "PASSED",
                    "checklist_species": gap_count,
                    "audit": audit_stats,
                    "primers_evaluated": primer_stats["evaluations_stored"],
                    "classifier_runs": len(bench_stats),
                    "exports": ver_stats,
                    "web_portal": web_stats,
                },
                indent=2,
            )
        )


def main() -> None:
    parser = argparse.ArgumentParser(prog="bengalfishref", description="BengalFishRef v1.0 Pipeline CLI (FBD-NIB)")
    sub = parser.add_subparsers(dest="command", required=True)

    p_chk = sub.add_parser("build-checklist", help="Build harmonized Bangladesh & Bay of Bengal checklist")
    p_chk.add_argument("--db", default=str(DEFAULT_DB))
    p_chk.add_argument("--export-tsv", default=str(DEFAULT_CHECKLIST_TSV))
    p_chk.add_argument("--offline-only", action="store_true")
    p_chk.set_defaults(func=cmd_build_checklist)

    p_ing = sub.add_parser("ingest-sequences", help="Resumably ingest NCBI + BOLD public sequences")
    p_ing.add_argument("--db", default=str(DEFAULT_DB))
    p_ing.add_argument("--markers", default="COI,12S,16S,Cytb,Mitogenome")
    p_ing.add_argument("--batch-size", type=int, default=8)
    p_ing.add_argument("--retmax", type=int, default=120)
    p_ing.add_argument("--max-batches", type=int, default=None)
    p_ing.add_argument("--include-bold", action="store_true")
    p_ing.add_argument("--max-bold-genera", type=int, default=None)
    p_ing.set_defaults(func=cmd_ingest_sequences)

    p_aud = sub.add_parser("run-audit", help="Run MAFFT + K2P + IQ-TREE2 phylogenetic mislabeling audit")
    p_aud.add_argument("--db", default=str(DEFAULT_DB))
    p_aud.add_argument("--threads", type=int, default=4)
    p_aud.add_argument("--iqtree", action="store_true")
    p_aud.set_defaults(
        func=lambda a: print(json.dumps(run_db_phylogenetic_audit(connect_db(a.db), a.threads, a.iqtree), indent=2))
    )

    p_prm = sub.add_parser("benchmark-primers", help="Run 14-primer in silico thermodynamic & resolution audit")
    p_prm.add_argument("--db", default=str(DEFAULT_DB))
    p_prm.set_defaults(func=lambda a: print(json.dumps(run_db_primer_benchmark(connect_db(a.db)), indent=2)))

    p_bch = sub.add_parser("benchmark-classifiers", help="Run 4-ecosystem simulated eDNA classifier benchmark")
    p_bch.add_argument("--db", default=str(DEFAULT_DB))
    p_bch.add_argument("--threads", type=int, default=4)
    p_bch.set_defaults(
        func=lambda a: print(json.dumps({"rows": len(run_db_classifier_benchmark(connect_db(a.db), threads=a.threads))}, indent=2))
    )

    p_exp = sub.add_parser("export-all", help="Export QIIME2, DADA2, SINTAX, BLAST+ & Serverless Web Portal")
    p_exp.add_argument("--db", default=str(DEFAULT_DB))
    p_exp.add_argument("--release-dir", default=str(DEFAULT_RELEASE_DIR))
    p_exp.add_argument("--web-dir", default=str(DEFAULT_WEB_DIR))
    p_exp.set_defaults(
        func=lambda a: print(
            json.dumps(
                {
                    "releases": export_bioinformatics_releases(connect_db(a.db), Path(a.release_dir)),
                    "web": build_static_web_portal(connect_db(a.db), Path(a.web_dir)),
                },
                indent=2,
            )
        )
    )

    p_ver = sub.add_parser("verify-exports", help="Verify exported release artifacts")
    p_ver.add_argument("--release-dir", default=str(DEFAULT_RELEASE_DIR))
    p_ver.set_defaults(func=lambda a: print(json.dumps(verify_exported_releases(Path(a.release_dir)), indent=2)))

    p_val = sub.add_parser("validate", help="Validate local FBD-NIB FASTA or .ab1 sequence file")
    p_val.add_argument("sample_file")
    p_val.add_argument("--db", default=str(DEFAULT_DB))
    p_val.add_argument("--output", "-o", default=None, help="Optional JSON output file path")
    def _run_val(a: argparse.Namespace) -> None:
        res = validate_lab_sample_file(connect_db(a.db), Path(a.sample_file))
        txt = json.dumps(res, indent=2)
        if a.output:
            Path(a.output).parent.mkdir(parents=True, exist_ok=True)
            Path(a.output).write_text(txt, encoding="utf-8")
            print(f"Success! Validation report written to: {a.output}")
        else:
            print(txt)
    p_val.set_defaults(func=_run_val)

    p_smk = sub.add_parser("smoke-test", help="Run deterministic offline 6-species end-to-end smoke test")
    p_smk.add_argument("--threads", type=int, default=4)
    p_smk.set_defaults(func=cmd_smoke_test)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
