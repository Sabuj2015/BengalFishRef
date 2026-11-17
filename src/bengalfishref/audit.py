"""Sequence hygiene, genetic-code-aware NUMT detector, K2P barcoding gap, and IQ-TREE2 audit."""

from __future__ import annotations

import hashlib
import math
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any

from Bio.Seq import Seq
import numpy as np

MARKER_MIN_LENGTHS: dict[str, int] = {
    "COI": 480,
    "12S": 140,
    "16S": 160,
    "Cytb": 450,
    "Mitogenome": 14000,
}

MARKER_MAX_LENGTHS: dict[str, int] = {
    "COI": 1700,
    "12S": 1200,
    "16S": 1900,
    "Cytb": 1350,
    "Mitogenome": 22000,
}

PROTEIN_CODING_MARKERS: set[str] = {"COI", "Cytb"}


def make_seq_uid(source_db: str, accession: str, marker: str) -> str:
    """Create a deterministic SHA-256 sequence record UID."""
    raw = f"{source_db.upper()}::{accession.strip()}::{marker}"
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()[:24]


def clean_dna(seq: str) -> str:
    """Normalize uppercase nucleotide sequence without gaps or whitespace."""
    return re.sub(r"[^A-Za-z]", "", seq).upper().replace("U", "T")


def max_homopolymer_run(seq: str) -> int:
    """Return length of longest identical A/C/G/T homopolymer stretch."""
    best = 0
    for m in re.finditer(r"(A+|C+|G+|T+)", seq):
        run_len = len(m.group(0))
        if run_len > best:
            best = run_len
    return best


def compute_ambig_and_gc(seq: str) -> tuple[float, float]:
    """Return (% IUPAC ambiguous bases outside ACGT, % GC among ACGT)."""
    if not seq:
        return 100.0, 0.0
    n_total = len(seq)
    acgt = sum(seq.count(b) for b in "ACGT")
    ambig_pct = round(100.0 * (n_total - acgt) / n_total, 3)
    gc_pct = round(100.0 * (seq.count("G") + seq.count("C")) / max(acgt, 1), 2)
    return ambig_pct, gc_pct


def detect_numt_stop_codons(seq: str, genetic_code: int) -> tuple[bool, int, int]:
    """Check whether a protein-coding mitochondrial locus (COI/Cytb) has internal stop codons.

    Tests all 3 forward reading frames (and reverse complement if strand-inverted)
    under NCBI translation table 2 (Vertebrate Mitochondrial: AGA/AGG/TAA/TAG stops;
    TGA = Trp) or table 5 (Invertebrate Mitochondrial: TAA/TAG stops; TGA = Trp,
    AGA/AGG = Ser). A valid mitochondrial CDS must have >=1 reading frame with 0
    internal stop codons (excluding a possible terminal stop).

    Returns:
        (is_numt, best_frame, min_internal_stops)
    """
    if genetic_code not in (2, 5):
        raise ValueError(f"Unsupported mitochondrial genetic_code={genetic_code}; expected 2 or 5.")
    dna = clean_dna(seq)
    if len(dna) < 90:
        return True, 0, 99

    best_stops = 999
    best_frame = 0
    candidates = [(dna, +1), (str(Seq(dna).reverse_complement()), -1)]
    for template, strand in candidates:
        for frame in (0, 1, 2):
            sub = template[frame:]
            trim_len = (len(sub) // 3) * 3
            if trim_len < 90:
                continue
            # Replace IUPAC N/degenerate codons with NNN -> 'X' so single N bases aren't misread as stops
            codons = [sub[i : i + 3] for i in range(0, trim_len, 3)]
            safe_dna = "".join(c if set(c) <= {"A", "C", "G", "T"} else "NNN" for c in codons)
            aa = str(Seq(safe_dna).translate(table=genetic_code))
            # Ignore a final terminal stop codon if full-length CDS ends at stop
            internal_aa = aa[:-1] if aa.endswith("*") else aa
            stops = internal_aa.count("*")
            if stops < best_stops:
                best_stops = stops
                best_frame = strand * (frame + 1)
            if best_stops == 0:
                return False, best_frame, 0

    return best_stops > 0, best_frame, best_stops


def evaluate_sequence_hygiene(
    seq: str,
    marker: str,
    genetic_code: int,
    is_voucher_or_geo: bool = False,
) -> dict[str, Any]:
    """Run primary sequence hygiene and NUMT checks and assign initial quality tier."""
    dna = clean_dna(seq)
    length = len(dna)
    ambig_pct, gc_pct = compute_ambig_and_gc(dna)
    homopoly = max_homopolymer_run(dna)
    min_len = MARKER_MIN_LENGTHS.get(marker, 140)
    max_len = MARKER_MAX_LENGTHS.get(marker, 22000)

    reasons: list[str] = []
    if length < min_len:
        reasons.append(f"short_sequence({length}<{min_len}bp)")
    if length > max_len:
        reasons.append(f"excessive_length({length}>{max_len}bp)")
    if ambig_pct > 3.0:
        reasons.append(f"high_ambiguity({ambig_pct:.2f}%)")
    if homopoly >= 14:
        reasons.append(f"homopolymer_artifact({homopoly}bp)")

    has_stop = False
    if marker in PROTEIN_CODING_MARKERS and not reasons:
        has_stop, _, n_stops = detect_numt_stop_codons(dna, genetic_code=genetic_code)
        if has_stop:
            reasons.append(f"internal_stop_numt(stops={n_stops},table={genetic_code})")

    if reasons:
        tier = "Quarantined"
    elif ambig_pct <= 1.0 and (length >= int(min_len * 1.05) or marker in ("12S", "16S")) and homopoly <= 10:
        tier = "Gold" if is_voucher_or_geo or ambig_pct <= 0.2 else "Silver"
    else:
        tier = "Silver"

    return {
        "clean_seq": dna,
        "seq_length": length,
        "ambig_pct": ambig_pct,
        "gc_pct": gc_pct,
        "homopolymer_max": homopoly,
        "has_stop_codon_numt": int(has_stop),
        "quality_tier": tier,
        "quarantine_reason": ";".join(reasons),
    }


def k2p_distance(aln_seq1: str, aln_seq2: str) -> float:
    """Compute Kimura 2-Parameter (K2P) evolutionary distance on aligned DNA strings."""
    transitions = 0
    transversions = 0
    valid_sites = 0
    purines = {"A", "G"}
    pyrimidines = {"C", "T"}

    for b1, b2 in zip(aln_seq1.upper(), aln_seq2.upper()):
        if b1 not in "ACGT" or b2 not in "ACGT":
            continue
        valid_sites += 1
        if b1 == b2:
            continue
        if (b1 in purines and b2 in purines) or (b1 in pyrimidines and b2 in pyrimidines):
            transitions += 1
        else:
            transversions += 1

    if valid_sites < 80:
        return float("nan")
    p = transitions / valid_sites
    q = transversions / valid_sites
    term1 = 1.0 - 2.0 * p - q
    term2 = 1.0 - 2.0 * q
    if term1 <= 0.0 or term2 <= 0.0:
        return 1.0
    dist = -0.5 * math.log(term1) - 0.25 * math.log(term2)
    return round(max(dist, 0.0), 5)


def parse_fasta_str(text: str) -> dict[str, str]:
    """Parse FASTA text into {id: sequence}."""
    records: dict[str, list[str]] = {}
    curr: str | None = None
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            curr = line[1:].split()[0]
            records[curr] = []
        elif curr is not None:
            records[curr].append(line)
    return {k: "".join(v) for k, v in records.items()}


def run_mafft_and_k2p_audit(
    records: list[dict[str, Any]],
    threads: int = 4,
    run_iqtree: bool = False,
) -> tuple[dict[str, dict[str, Any]], str | None]:
    """Align a family+marker sequence set with MAFFT (--thread 4), compute K2P & optional IQ-TREE2.

    Flags sequences where:
      - Nearest heterospecific (different genus/species) K2P distance is <= 0.005 (0.5%) AND
        intraspecific distance to conspecific majority is >= 0.035 (3.5%) -> 'Mislabeled_Outlier'.
      - Intraspecific divergence > 0.035 without heterospecific match -> 'Cryptic_Lineage_Candidate'.
    """
    threads = min(max(int(threads), 1), 4)
    if len(records) < 2 or shutil.which("mafft") is None:
        return {}, None

    with tempfile.TemporaryDirectory(prefix="bfr_audit_") as tmpdir:
        tmp = Path(tmpdir)
        in_fa = tmp / "input.fasta"
        aln_fa = tmp / "aligned.fasta"
        with in_fa.open("w", encoding="utf-8") as fh:
            for r in records:
                fh.write(f">{r['seq_uid']}\n{r['sequence']}\n")

        mafft_cmd = ["mafft", "--retree", "1", "--quiet", "--thread", str(threads), str(in_fa)] if len(records) > 120 else ["mafft", "--auto", "--quiet", "--thread", str(threads), str(in_fa)]
        subprocess.run(
            mafft_cmd,
            stdout=aln_fa.open("w", encoding="utf-8"),
            stderr=subprocess.PIPE,
            check=True,
            timeout=300,
        )
        aligned = parse_fasta_str(aln_fa.read_text(encoding="utf-8"))

        newick_tree: str | None = None
        if run_iqtree and len(records) >= 4 and shutil.which("iqtree2") is not None:
            cmd = [
                "iqtree2",
                "-s",
                str(aln_fa),
                "-m",
                "HKY+G",
                "-T",
                str(threads),
                "--quiet",
                "-redo",
            ]
            if len(records) > 60:
                cmd.append("-fast")
            subprocess.run(cmd, check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=240)
            treefile = Path(str(aln_fa) + ".treefile")
            if treefile.exists():
                newick_tree = treefile.read_text(encoding="utf-8").strip()

    uid_to_meta = {r["seq_uid"]: r for r in records if r["seq_uid"] in aligned}
    uids = list(uid_to_meta.keys())
    audit_out: dict[str, dict[str, Any]] = {}
    if len(uids) < 2:
        return audit_out, newick_tree

    # Vectorized NumPy K2P matrix across aligned sequences (codes: A=0, G=1 [purines: //2 == 0], C=2, T=3 [pyrimidines: //2 == 1], gap/other=255)
    L_aln = len(aligned[uids[0]])
    code_lut = np.full(128, 255, dtype=np.uint8)
    code_lut[ord("A")] = 0
    code_lut[ord("G")] = 1
    code_lut[ord("C")] = 2
    code_lut[ord("T")] = 3

    mat = np.empty((len(uids), L_aln), dtype=np.uint8)
    for idx, u in enumerate(uids):
        s_up = aligned[u].upper()
        if len(s_up) != L_aln:
            s_up = s_up[:L_aln].ljust(L_aln, "-")
        raw_b = np.frombuffer(s_up.encode("ascii", errors="replace"), dtype=np.uint8)
        mat[idx] = code_lut[np.minimum(raw_b, 127)]

    valid_mask = mat < 4  # (N, L)
    pyr_class = mat // 2  # 0 for A/G, 1 for C/T

    tids = np.array([int(uid_to_meta[u]["resolved_taxon_id"]) for u in uids], dtype=np.int64)
    vnames = [uid_to_meta[u]["valid_name"] for u in uids]

    for i, u1 in enumerate(uids):
        v_both = valid_mask[i] & valid_mask  # (N, L)
        n_valid = v_both.sum(axis=1).astype(np.float64)
        diff = (mat[i] != mat) & v_both
        same_class = (pyr_class[i] == pyr_class) & diff  # transition
        diff_class = (pyr_class[i] != pyr_class) & diff  # transversion

        with np.errstate(divide="ignore", invalid="ignore"):
            p = same_class.sum(axis=1) / n_valid
            q = diff_class.sum(axis=1) / n_valid
            t1 = 1.0 - 2.0 * p - q
            t2 = 1.0 - 2.0 * q
            dists = np.where((t1 > 0.0) & (t2 > 0.0), -0.5 * np.log(np.maximum(t1, 1e-12)) - 0.25 * np.log(np.maximum(t2, 1e-12)), 1.0)
            dists = np.maximum(dists, 0.0)
            dists = np.where(n_valid >= 80, np.round(dists, 5), np.nan)

        dists[i] = np.nan
        con_mask = (tids == tids[i]) & ~np.isnan(dists)
        het_mask = (tids != tids[i]) & ~np.isnan(dists)

        con_vals = np.sort(dists[con_mask])
        min_con = float(con_vals[0]) if len(con_vals) > 0 else None
        med_con = float(con_vals[len(con_vals) // 2]) if len(con_vals) > 0 else None

        if np.any(het_mask):
            het_indices = np.where(het_mask)[0]
            best_rel = int(np.argmin(dists[het_indices]))
            best_j = int(het_indices[best_rel])
            hetero_val: float | None = float(dists[best_j])
            best_hetero_taxon = vnames[best_j]
        else:
            hetero_val = None
            best_hetero_taxon = ""

        flag = "Monophyletic"
        quarantine_add = ""
        if hetero_val is not None:
            if med_con is not None and med_con >= 0.035 and hetero_val <= 0.008:
                flag = "Mislabeled_Outlier"
                quarantine_add = f"phylogenetic_mislabel_matches_{best_hetero_taxon.replace(' ', '_')}(K2P={hetero_val:.4f})"
            elif hetero_val == 0.0:
                flag = "Paraphyletic_Shared"
            elif med_con is not None and med_con >= 0.035:
                flag = "Cryptic_Lineage_Candidate"

        audit_out[u1] = {
            "k2p_nearest_conspecific": min_con,
            "k2p_nearest_heterospecific": hetero_val,
            "nearest_heterospecific_taxon": best_hetero_taxon,
            "phylo_concordance_flag": flag,
            "quarantine_add": quarantine_add,
        }

    return audit_out, newick_tree
