(function () {
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
    document.getElementById("seqInput").value = `>Query_${sp.replace(/\s+/g, "_")}
${seq}`;
  });

  document.getElementById("runIdBtn").addEventListener("click", () => {
    const raw = document.getElementById("seqInput").value || "";
    const dna = raw.split(/\r?\n/).filter(l => !l.startsWith(">")).join("").replace(/[^A-Za-z]/g, "").toUpperCase();
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
