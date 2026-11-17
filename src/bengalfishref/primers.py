"""In silico eDNA primer thermodynamics, 3'-mismatch audit, taxonomic resolution, and degeneracy optimizer."""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Any

from Bio.Seq import Seq
import numpy as np

_ASCII_TO_MASK = np.zeros(128, dtype=np.uint8)
for _ch, _bases in {
    "A": 1, "C": 2, "G": 4, "T": 8,
    "R": 5, "Y": 10, "S": 6, "W": 9, "K": 12, "M": 3,
    "B": 14, "D": 13, "H": 11, "V": 7, "N": 15, "I": 15,
}.items():
    _ASCII_TO_MASK[ord(_ch)] = _bases

IUPAC_DNA: dict[str, set[str]] = {
    "A": {"A"},
    "C": {"C"},
    "G": {"G"},
    "T": {"T"},
    "R": {"A", "G"},
    "Y": {"C", "T"},
    "S": {"G", "C"},
    "W": {"A", "T"},
    "K": {"G", "T"},
    "M": {"A", "C"},
    "B": {"C", "G", "T"},
    "D": {"A", "G", "T"},
    "H": {"A", "C", "T"},
    "V": {"A", "C", "G"},
    "N": {"A", "C", "G", "T"},
    "I": {"A", "C", "G", "T"},  # Inosine universal pairing in jgHCO2198
}

REVERSE_IUPAC: dict[frozenset[str], str] = {
    frozenset(v): k for k, v in IUPAC_DNA.items() if k != "I"
}


@dataclass(frozen=True)
class PrimerPair:
    name: str
    marker: str
    target_group: str
    fwd_seq: str
    rev_seq: str
    min_amplicon: int
    max_amplicon: int
    reference: str


STANDARD_EDNA_PRIMERS: tuple[PrimerPair, ...] = (
    PrimerPair(
        name="MiFish-U",
        marker="12S",
        target_group="Teleostei",
        fwd_seq="GTCGGTAAAACTCGTGCCAGC",
        rev_seq="CATAGTGGGGTATCTAATCCCAGTTTG",
        min_amplicon=130,
        max_amplicon=230,
        reference="Miya et al. (2015) R. Soc. Open Sci.",
    ),
    PrimerPair(
        name="MiFish-E",
        marker="12S",
        target_group="Elasmobranchii",
        fwd_seq="GTTGGTAAATCTCGTGCCAGC",
        rev_seq="CATAGTGGGGTATCTAATCCTAGTTTG",
        min_amplicon=130,
        max_amplicon=230,
        reference="Miya et al. (2015) R. Soc. Open Sci.",
    ),
    PrimerPair(
        name="Teleo",
        marker="12S",
        target_group="Teleostei",
        fwd_seq="ACACCGCCCGTCACTCT",
        rev_seq="CTTCCGGTACACTTACCATG",
        min_amplicon=45,
        max_amplicon=115,
        reference="Valentini et al. (2016) Mol. Ecol.",
    ),
    PrimerPair(
        name="AcMDB07",
        marker="12S",
        target_group="Actinopterygii",
        fwd_seq="GCCTATATACCGCCGTCG",
        rev_seq="GTACACTTACCATGTTACGACTT",
        min_amplicon=240,
        max_amplicon=360,
        reference="Bylemans et al. (2018) Metabarcoding Metagenom.",
    ),
    PrimerPair(
        name="EcoPrimers-12S-V5",
        marker="12S",
        target_group="Vertebrata",
        fwd_seq="ACTGGGATTAGATACCCC",
        rev_seq="TAGAACAGGCTCCTCTAG",
        min_amplicon=70,
        max_amplicon=140,
        reference="Riaz et al. (2011) Nucleic Acids Res.",
    ),
    PrimerPair(
        name="MarVer1",
        marker="12S",
        target_group="Marine_Vertebrata",
        fwd_seq="CGTGCCAGCCACCGCG",
        rev_seq="GGGTATCTAATCCCAGTTTG",
        min_amplicon=150,
        max_amplicon=250,
        reference="Valsecchi et al. (2020) Mol. Ecol. Resour.",
    ),
    PrimerPair(
        name="Tele02",
        marker="12S",
        target_group="Teleostei",
        fwd_seq="AAACTCGTGCCAGCCACC",
        rev_seq="GGGTATCTAATCCCAGTTTG",
        min_amplicon=130,
        max_amplicon=220,
        reference="Taberlet et al. (2018) Oxford Univ. Press",
    ),
    PrimerPair(
        name="Fish16S",
        marker="16S",
        target_group="Pisces",
        fwd_seq="GACCCTATGGAGCTTTAGAC",
        rev_seq="CGCTGTTATCCCTADRGTAACT",
        min_amplicon=160,
        max_amplicon=280,
        reference="Berry et al. (2017) Mol. Ecol. Resour.",
    ),
    PrimerPair(
        name="MarVer3",
        marker="16S",
        target_group="Marine_Vertebrata",
        fwd_seq="AGACGAGAAGACCCTRTG",
        rev_seq="GGATTGCGCTGTTATCCC",
        min_amplicon=190,
        max_amplicon=300,
        reference="Valsecchi et al. (2020) Mol. Ecol. Resour.",
    ),
    PrimerPair(
        name="Crust16S",
        marker="16S",
        target_group="Crustacea",
        fwd_seq="GGGACGATAAGACCCTATA",
        rev_seq="ATTACGCTGTTATCCCTAAAG",
        min_amplicon=130,
        max_amplicon=240,
        reference="Berry et al. (2017) Mol. Ecol. Resour.",
    ),
    PrimerPair(
        name="Deca16S",
        marker="16S",
        target_group="Crustacea",
        fwd_seq="TGCCTGTTTATCAAAAACAT",
        rev_seq="CTCCGGTTTGAACTCAGATCA",
        min_amplicon=420,
        max_amplicon=580,
        reference="Palumbi / Bracken-Grissom Decapod 16S",
    ),
    PrimerPair(
        name="MiDeca",
        marker="16S",
        target_group="Crustacea",
        fwd_seq="GGACGATAAGACCCTATAAA",
        rev_seq="ACGCTGTTATCCCTAAAGT",
        min_amplicon=135,
        max_amplicon=215,
        reference="Komai et al. (2019) Mar. Biodivers.",
    ),
    PrimerPair(
        name="Leray-COI",
        marker="COI",
        target_group="Metazoa",
        fwd_seq="GGWACWGGWTGAACWGTWTAYCCYCC",
        rev_seq="TAIACYTCIGGRTGICCRAARAAYCA",
        min_amplicon=290,
        max_amplicon=340,
        reference="Leray et al. (2013) Front. Zool.",
    ),
    PrimerPair(
        name="Ward-FishCOI",
        marker="COI",
        target_group="Pisces",
        fwd_seq="TCAACCAACCACAAAGACATTGGCAC",
        rev_seq="TAGACTTCTGGGTGGCCAAAGAATCA",
        min_amplicon=620,
        max_amplicon=690,
        reference="Ward et al. (2005) Philos. Trans. R. Soc. B",
    ),
)


def base_matches(primer_base: str, template_base: str) -> bool:
    """Return True if IUPAC primer nucleotide matches IUPAC template nucleotide."""
    p_set = IUPAC_DNA.get(primer_base.upper(), set())
    t_set = IUPAC_DNA.get(template_base.upper(), set())
    return bool(p_set & t_set)


def score_primer_binding_sense(primer_5to3: str, target_sense_5to3: str) -> dict[str, Any]:
    """Score a primer (5'->3') against an equal-length target binding site (same 5'->3' sense).

    Penalizes 3'-pentamer mismatches heavily because Taq/Q5 polymerases stall on 3' mismatches:
      - distance 1 from 3' end (ultimate 3' base): weight 4.0
      - distance 2: weight 3.0
      - distance 3: weight 2.5
      - distance 4..5: weight 2.0
      - internal (>5 from 3' end): weight 1.0
    """
    L = len(primer_5to3)
    if len(target_sense_5to3) != L:
        raise ValueError("Primer and target window must have identical length.")

    total_mm = 0
    three_prime_mm = 0
    weighted_penalty = 0.0
    mismatch_positions_from_3prime: list[int] = []

    for idx, (pb, tb) in enumerate(zip(primer_5to3.upper(), target_sense_5to3.upper())):
        dist_3p = L - idx  # 1 = 3' terminal nucleotide
        if not base_matches(pb, tb):
            total_mm += 1
            mismatch_positions_from_3prime.append(dist_3p)
            if dist_3p == 1:
                three_prime_mm += 1
                weighted_penalty += 4.0
            elif dist_3p == 2:
                three_prime_mm += 1
                weighted_penalty += 3.0
            elif dist_3p == 3:
                three_prime_mm += 1
                weighted_penalty += 2.5
            elif dist_3p in (4, 5):
                three_prime_mm += 1
                weighted_penalty += 2.0
            else:
                weighted_penalty += 1.0

    # Empirical nearest-neighbor approximation: ~2.2 C drop per internal mismatch, ~4.0 C near 3' end
    delta_tm = round(2.0 * (total_mm - three_prime_mm) + 3.8 * three_prime_mm, 2)
    return {
        "mismatches": total_mm,
        "three_prime_mm": three_prime_mm,
        "weighted_penalty": round(weighted_penalty, 2),
        "delta_tm": delta_tm,
        "mismatch_3p_offsets": mismatch_positions_from_3prime,
        "target_window": target_sense_5to3.upper(),
    }


def scan_best_binding_site(primer_5to3: str, template_5to3: str, reverse_primer: bool = False) -> dict[str, Any] | None:
    """Find the lowest-penalty binding site of a primer along a template sequence (5'->3') via vectorized bitmasks."""
    L = len(primer_5to3)
    T = len(template_5to3)
    if T < L:
        return None

    seq_upper = template_5to3.upper()
    t_bytes = np.frombuffer(seq_upper.encode("ascii", errors="replace"), dtype=np.uint8)
    t_masks = _ASCII_TO_MASK[np.minimum(t_bytes, 127)]
    windows = np.lib.stride_tricks.sliding_window_view(t_masks, L)  # shape (T - L + 1, L)

    if not reverse_primer:
        p_str = primer_5to3.upper()
        p_masks = _ASCII_TO_MASK[np.frombuffer(p_str.encode("ascii"), dtype=np.uint8)]
        weights = np.array(
            [4.0 if (L - i) == 1 else 3.0 if (L - i) == 2 else 2.5 if (L - i) == 3 else 2.0 if (L - i) in (4, 5) else 1.0 for i in range(L)],
            dtype=np.float32,
        )
    else:
        # Reverse primer anneals antiparallel; sense-strand match is reverse_complement(primer_5to3),
        # where window index 0 corresponds to the ultimate 3' nucleotide (dist_3p = 1) of the reverse primer!
        p_rc = str(Seq(primer_5to3.upper()).reverse_complement())
        p_masks = _ASCII_TO_MASK[np.frombuffer(p_rc.encode("ascii"), dtype=np.uint8)]
        weights = np.array(
            [4.0 if (i + 1) == 1 else 3.0 if (i + 1) == 2 else 2.5 if (i + 1) == 3 else 2.0 if (i + 1) in (4, 5) else 1.0 for i in range(L)],
            dtype=np.float32,
        )

    mismatches_bool = (windows & p_masks) == 0  # True where IUPAC masks do not overlap
    penalties = mismatches_bool @ weights
    mm_counts = mismatches_bool.sum(axis=1)
    # Lexicographic minimization of (weighted_penalty, total_mismatches)
    composite = penalties * 100.0 + mm_counts
    best_start = int(np.argmin(composite))

    win_str = seq_upper[best_start : best_start + L]
    target_equiv = str(Seq(win_str).reverse_complement()) if reverse_primer else win_str
    scored = score_primer_binding_sense(primer_5to3, target_equiv)
    return {
        **scored,
        "start": best_start,
        "end": best_start + L,
    }


def evaluate_primer_pair_on_template(pair: PrimerPair, template_5to3: str) -> dict[str, Any] | None:
    """Evaluate forward + reverse primer binding, amplification plausibility, and amplicon sequence."""
    fwd = scan_best_binding_site(pair.fwd_seq, template_5to3, reverse_primer=False)
    if fwd is None or fwd["mismatches"] > 6:
        return None

    # Search reverse primer downstream of forward primer 3' end
    search_start = fwd["end"]
    downstream = template_5to3[search_start : search_start + pair.max_amplicon + len(pair.rev_seq) + 120]
    rev = scan_best_binding_site(pair.rev_seq, downstream, reverse_primer=True)
    if rev is None or rev["mismatches"] > 6:
        return None

    amplicon_seq = downstream[: rev["start"]]
    amp_len = len(amplicon_seq)
    if not (max(25, pair.min_amplicon - 30) <= amp_len <= pair.max_amplicon + 50):
        return None

    total_penalty = round(fwd["weighted_penalty"] + rev["weighted_penalty"], 2)
    delta_tm_max = max(fwd["delta_tm"], rev["delta_tm"])

    # Realistic PCR amplification rule: <=3 total mismatches per primer, <=1 mismatch in 3'-pentamer,
    # and no mismatch at the ultimate 3' nucleotide (offset 1).
    ult_3p_clash = (1 in fwd["mismatch_3p_offsets"]) or (1 in rev["mismatch_3p_offsets"])
    amplifies = int(
        not ult_3p_clash
        and fwd["mismatches"] <= 3
        and rev["mismatches"] <= 3
        and fwd["three_prime_mm"] <= 1
        and rev["three_prime_mm"] <= 1
        and total_penalty <= 8.0
    )

    return {
        "primer_name": pair.name,
        "marker": pair.marker,
        "fwd_mismatches": fwd["mismatches"],
        "rev_mismatches": rev["mismatches"],
        "fwd_3prime_mm": fwd["three_prime_mm"],
        "rev_3prime_mm": rev["three_prime_mm"],
        "weighted_penalty": total_penalty,
        "delta_tm_est": delta_tm_max,
        "amplifies_in_silico": amplifies,
        "amplicon_length": amp_len,
        "amplicon_seq": amplicon_seq,
        "fwd_target_window": fwd["target_window"],
        "rev_target_window": rev["target_window"],
    }


def compute_amplicon_species_resolution(
    taxon_amplicons: list[tuple[int, str, str, str]],
) -> dict[str, Any]:
    """Compute species- and genus-level distinguishability of extracted inter-primer amplicons.

    Args:
        taxon_amplicons: list of (taxon_id, valid_name, genus_name, amplicon_seq)
    """
    if not taxon_amplicons:
        return {"n_taxa": 0, "species_resolved_pct": 0.0, "genus_resolved_pct": 0.0, "collisions": []}

    seq_to_taxa: dict[str, set[int]] = defaultdict(set)
    seq_to_names: dict[str, set[str]] = defaultdict(set)
    seq_to_genera: dict[str, set[str]] = defaultdict(set)
    all_taxa: set[int] = set()

    for tid, vname, gname, amp in taxon_amplicons:
        if not amp:
            continue
        all_taxa.add(tid)
        seq_to_taxa[amp].add(tid)
        seq_to_names[amp].add(vname)
        seq_to_genera[amp].add(gname)

    if not all_taxa:
        return {"n_taxa": 0, "species_resolved_pct": 0.0, "genus_resolved_pct": 0.0, "collisions": []}

    unambiguous_species: set[int] = set()
    unambiguous_genus: set[int] = set()
    collisions: list[dict[str, Any]] = []

    # A species is uniquely resolved if none of its haplotypes are shared with another species
    shared_taxa: set[int] = set()
    shared_across_genera: set[int] = set()
    for amp, tids in seq_to_taxa.items():
        if len(tids) > 1:
            shared_taxa.update(tids)
            collisions.append(
                {
                    "species": sorted(seq_to_names[amp]),
                    "genera": sorted(seq_to_genera[amp]),
                    "amplicon_len": len(amp),
                }
            )
            if len(seq_to_genera[amp]) > 1:
                shared_across_genera.update(tids)

    unambiguous_species = all_taxa - shared_taxa
    unambiguous_genus = all_taxa - shared_across_genera

    return {
        "n_taxa": len(all_taxa),
        "n_species_resolved": len(unambiguous_species),
        "species_resolved_pct": round(100.0 * len(unambiguous_species) / len(all_taxa), 2),
        "genus_resolved_pct": round(100.0 * len(unambiguous_genus) / len(all_taxa), 2),
        "collisions": collisions[:25],
    }


def design_bengal_degenerate_primer(
    original_primer: str,
    target_windows: list[str],
    min_allele_freq: float = 0.05,
    max_degeneracy: int = 16,
) -> dict[str, Any]:
    """Design a minimal-degeneracy Bengal-adapted primer rescuing South Asian 3' drop-out clades."""
    L = len(original_primer)
    valid_windows = [w.upper() for w in target_windows if len(w) == L]
    if not valid_windows:
        return {"primer": original_primer, "degeneracy": 1, "modified_positions": []}

    n = len(valid_windows)
    # Prioritize rescuing positions closest to the 3' end first within max_degeneracy budget
    pos_order = list(range(L - 1, -1, -1))
    chosen_sets: list[set[str]] = [set(IUPAC_DNA.get(b, {b})) for b in original_primer.upper()]
    modified: list[dict[str, Any]] = []

    def current_degeneracy(sets: list[set[str]]) -> int:
        deg = 1
        for s in sets:
            deg *= max(len(s), 1)
        return deg

    for idx in pos_order:
        counts = Counter(w[idx] for w in valid_windows if w[idx] in "ACGT")
        needed = set(chosen_sets[idx])
        for base, c in counts.items():
            if c / n >= min_allele_freq and base not in needed:
                trial = [set(s) for s in chosen_sets]
                trial[idx] = needed | {base}
                if current_degeneracy(trial) <= max_degeneracy:
                    needed = trial[idx]
        if needed != chosen_sets[idx]:
            old_code = REVERSE_IUPAC.get(frozenset(chosen_sets[idx]), original_primer[idx])
            new_code = REVERSE_IUPAC.get(frozenset(needed), "N")
            chosen_sets[idx] = needed
            modified.append(
                {
                    "pos_1based_5to3": idx + 1,
                    "offset_from_3prime": L - idx,
                    "old_base": old_code,
                    "new_base": new_code,
                }
            )

    new_primer = "".join(
        REVERSE_IUPAC.get(frozenset(s), original_primer[i]) for i, s in enumerate(chosen_sets)
    )
    return {
        "primer": new_primer,
        "degeneracy": current_degeneracy(chosen_sets),
        "modified_positions": modified,
    }
