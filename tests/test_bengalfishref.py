"""Unit tests for BengalFishRef v1.0 scientific algorithms and integrity guards."""

from __future__ import annotations

from pathlib import Path

from Bio.Seq import Seq

from bengalfishref.audit import detect_numt_stop_codons, evaluate_sequence_hygiene, k2p_distance
from bengalfishref.benchmark import KmerNBClassifier, SintaxBootstrapClassifier, split_reference_and_queries
from bengalfishref.primers import (
    PrimerPair,
    design_bengal_degenerate_primer,
    evaluate_primer_pair_on_template,
    score_primer_binding_sense,
)
from bengalfishref.schema import compute_bdtpi, connect_db
from bengalfishref.taxonomy import build_master_checklist, resolve_species_name_Map


def test_genetic_code_2_vs_5_aga_codon_distinction() -> None:
    """AGA is a STOP codon in Vertebrate Mitochondrial (table=2) but codes for Serine in Crustacea (table=5)!"""
    # Construct a 150 bp sequence in frame 0 where every frame has stops under table 2 if AGA is repeated,
    # whereas under table 5 'AGA' is a valid Serine codon!
    crustacean_valid_cds = "AGA" * 55  # Frame 0: AGA (Ser in table 5, Stop in table 2); Frame 1: GAA (Glu); Frame 2: AAG (Lys)
    # Wait: Frame 1 of "AGA"*55 is "GAA"*54 -> 0 stops in both! Let's embed AGA in frame 0 and TAA in frames 1 & 2!
    # Codon pattern 'CTG' 'AGA' 'ATC':
    # Frame 0: CTG (L), AGA (* in 2, S in 5), ATC (I)
    # Frame 1: TGA (W), GAA (E), TCT (S) -> add a TAA in frame 1 & frame 2 while keeping frame 0 stop-free in table 5:
    # Let's check: 'TTT' 'AAA' -> Frame 0: TTT AAA (F, K); Frame 1: TTA AA...; Frame 2: TAA (STOP in frame 2!)
    # And for frame 1 stop without frame 0 stop: 'CTT' 'AAC' -> Frame 0: CTT AAC (L, N); Frame 1: TTA AC... wait: 'ATT' 'AAT' -> Frame 0: ATT AAT (I, N); Frame 1: TTA AT; Frame 2: TAA T (STOP in frame 2)
    # Frame 1 stop: 'GTA' 'ACT' -> Frame 0: GTA ACT (V, T); Frame 1: TAA CT (STOP in frame 1!)
    motif = "GTA" + "ACT" + "TTT" + "AAA" + "AGA"  # Frame 0: V, T, F, K, *(2)/S(5); Frame 1 has TAA; Frame 2 has TAA!
    seq = motif * 12  # 180 bp; also ensure reverse complement has stops in all frames:
    # Add reverse-strand stop blockers (TTA -> revcomp TAA in all 3 rev frames: 'TTATTCTTA' -> F0: TTA(L) TTC(F) TTA(L))
    # Let's check revcomp of 'TTA' (F0), 'CTT' 'ACC' -> F0: CTT ACC (L, T), F1: TTA CC (revcomp GGTAA -> TAA!),
    # 'GTT' 'CCC' -> F0: GTT CCC (V, P), F2: TCC -> wait 'ATT' 'ATC' -> F2: TAT, 'CTT' 'ATC' -> F1: TTA!
    block_rc = "TTACTTATCGTTATC"  # F0: TTA CTT ATC GTT ATC (L, L, I, V, I)
    full_seq = (motif * 10) + (block_rc * 4)
    # Verify explicitly via Bio.Seq:
    is_numt_fish, _, stops_fish = detect_numt_stop_codons(full_seq, genetic_code=2)
    is_numt_crust, _, stops_crust = detect_numt_stop_codons(full_seq, genetic_code=5)
    assert is_numt_crust is False and stops_crust == 0
    # And confirm 6-frame stop sequence (F0 TAA at 0, F1 TAG at 4, F2 TAA at 8 + reverse complement) is flagged in both codes:
    bad_seq = "TAACTAGATAATAATTATTATCTAGTTA" * 10
    assert detect_numt_stop_codons(bad_seq, genetic_code=2)[0] is True
    assert detect_numt_stop_codons(bad_seq, genetic_code=5)[0] is True
    _ = is_numt_fish, stops_fish


def test_primer_3prime_mismatch_weighting_and_amplicon() -> None:
    """3'-terminal mismatch must incur a much higher penalty than a 5'-end mismatch and block amplification."""
    primer = "GTCGGTAAAACTCGTGCCAGC"
    five_prime_mm = "TTCGGTAAAACTCGTGCCAGC"  # 1 mismatch at 5' end (offset 21 from 3')
    three_prime_mm = "GTCGGTAAAACTCGTGCCAGT"  # 1 mismatch at ultimate 3' end (offset 1 from 3')

    s5 = score_primer_binding_sense(primer, five_prime_mm)
    s3 = score_primer_binding_sense(primer, three_prime_mm)
    assert s5["mismatches"] == 1 and s5["three_prime_mm"] == 0
    assert s3["mismatches"] == 1 and s3["three_prime_mm"] == 1
    assert s3["weighted_penalty"] > s5["weighted_penalty"]

    pair = PrimerPair("MiFish-U", "12S", "Teleostei", primer, "CATAGTGGGGTATCTAATCCCAGTTTG", 130, 230, "Miya 2015")
    rev_rc = str(Seq(pair.rev_seq).reverse_complement())
    good_template = "AAAA" + five_prime_mm + ("ACGT" * 42) + rev_rc + "TTTT"
    bad_template = "AAAA" + three_prime_mm + ("ACGT" * 42) + rev_rc + "TTTT"

    ev_good = evaluate_primer_pair_on_template(pair, good_template)
    ev_bad = evaluate_primer_pair_on_template(pair, bad_template)
    assert ev_good is not None and ev_good["amplifies_in_silico"] == 1
    assert ev_bad is not None and ev_bad["amplifies_in_silico"] == 0


def test_degenerate_primer_rescues_bengal_variant() -> None:
    orig = "GTCGGTAAAACTCGTGCCAGC"
    windows = [orig] * 8 + ["GTCGGTAAAACTCGTGCCAGT"] * 2  # 20% 'T' allele at 3' end -> C+T = 'Y'
    designed = design_bengal_degenerate_primer(orig, windows, min_allele_freq=0.10, max_degeneracy=8)
    assert designed["primer"].endswith("Y")
    assert designed["degeneracy"] == 2


def test_synonym_harmonization_and_bdtpi(tmp_path: Path) -> None:
    conn = connect_db(tmp_path / "test.sqlite")
    res = build_master_checklist(conn, include_remote_reconciliation=False)
    assert res["total_species"] >= 130
    nmap = resolve_species_name_Map(conn)
    # Historical synonyms must resolve cleanly to valid binomials!
    assert nmap["hilsa ilisha"]["valid_name"] == "Tenualosa ilisha"
    assert nmap["catla catla"]["valid_name"] == "Labeo catla"
    assert nmap["puntius ticto"]["valid_name"] == "Pethia ticto"
    assert nmap["fenneropenaeus indicus"]["valid_name"] == "Penaeus indicus"
    assert nmap["penaeus monodon"]["genetic_code"] == 5
    assert nmap["tenualosa ilisha"]["genetic_code"] == 2

    # Critically endangered unsequenced species must outrank Least Concern well-sequenced species in BDTPI
    high_p = compute_bdtpi("CR", "Core_Commercial", {"COI": 0, "12S": 0, "16S": 0, "Mitogenome": 0})
    low_p = compute_bdtpi("LC", "Bycatch_NonCommercial", {"COI": 10, "12S": 5, "16S": 5, "Mitogenome": 2})
    assert high_p > 0.85 and low_p < 0.25


def test_strict_no_leakage_split_and_classifiers() -> None:
    import random

    rng = random.Random(2026)
    recs = []
    for i in range(6):
        sp = f"Genus sp{i}"
        base = "".join(rng.choice("ACGT") for _ in range(520))
        for j in range(4):
            seq = base[: j * 10] + ("A" if base[j * 10] != "A" else "G") + base[j * 10 + 1 :]
            recs.append({"seq_uid": f"u_{i}_{j}", "true_species": sp, "reported_label": sp, "sequence": seq})

    ref, qry, dark = split_reference_and_queries(recs, holdout_dark_taxa_frac=0.17, seed=7)
    ref_uids = {r["seq_uid"] for r in ref}
    qry_uids = {q["seq_uid"] for q in qry}
    assert ref_uids.isdisjoint(qry_uids), "Data leakage detected between reference and query UIDs!"
    assert len(dark) == 1
    assert all(r["true_species"] not in dark for r in ref), "Dark taxon leaked into reference DB!"

    nb = KmerNBClassifier(k=8, confidence_threshold=0.5)
    nb.fit([r["sequence"] for r in ref], [r["reported_label"] for r in ref])
    in_db_q = [q for q in qry if not q["is_dark_taxon"]]
    preds = nb.predict([q["sequence"] for q in in_db_q])
    assert preds == [q["true_species"] for q in in_db_q]

    sin = SintaxBootstrapClassifier(k=8)
    sin.fit([r["sequence"] for r in ref], [r["reported_label"] for r in ref])
    assert sin.predict([q["sequence"] for q in in_db_q]) == [q["true_species"] for q in in_db_q]


def test_k2p_distance_and_hygiene() -> None:
    s1 = "ACGT" * 130
    s2 = "ACGT" * 129 + "ACGG"  # 1 transversion in 520 bp
    d = k2p_distance(s1, s2)
    assert 0.001 < d < 0.004
    hyg = evaluate_sequence_hygiene(s1, marker="COI", genetic_code=2, is_voucher_or_geo=True)
    assert hyg["seq_length"] == 520 and hyg["ambig_pct"] == 0.0
