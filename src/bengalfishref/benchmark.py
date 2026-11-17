"""Simulated 4-ecosystem Bangladesh eDNA mock-community classifier benchmark (ml-best-practices compliant)."""

from __future__ import annotations

from collections import Counter
import random
import shutil
import subprocess
import tempfile
import time
from pathlib import Path
from typing import Any

import numpy as np
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.naive_bayes import MultinomialNB

ECOSYSTEM_DEFINITIONS: dict[str, str] = {
    "GBM_Freshwater": "Padma-Meghna-Jamuna River & Haor Freshwater Assemblage",
    "Sundarbans_Estuarine": "Sundarbans Mangrove & Brackish Estuarine Assemblage",
    "Bay_of_Bengal_Marine": "Northern Bay of Bengal Coastal, Pelagic & Elasmobranch Assemblage",
    "Aquaculture_Core": "Bangladesh Commercial Polyculture & Prawn/Shrimp Pond Assemblage",
}


def mutate_sequence_reads(seq: str, error_rate: float, rng: random.Random) -> str:
    """Simulate Illumina/ONT sequencing substitution + minor indel errors on an amplicon/barcode."""
    if error_rate <= 0.0:
        return seq
    bases = "ACGT"
    out: list[str] = []
    for ch in seq.upper():
        r = rng.random()
        if r < error_rate * 0.90:
            choices = [b for b in bases if b != ch]
            out.append(rng.choice(choices))
        elif r < error_rate * 0.95:
            # 1 bp deletion
            continue
        elif r < error_rate:
            # 1 bp insertion
            out.append(ch)
            out.append(rng.choice(bases))
        else:
            out.append(ch)
    return "".join(out)


def split_reference_and_queries(
    records: list[dict[str, Any]],
    holdout_dark_taxa_frac: float = 0.12,
    seed: int = 20261116,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], set[str]]:
    """Strictly partition accessions into disjoint Reference (Train) and Query (Test) sets BEFORE featurization.

    - Multi-accession species (>=2 records): 1+ records held out purely as test templates (in-DB species).
    - A subset (`holdout_dark_taxa_frac`) of species is completely withheld from the Reference DB to
      simulate unsequenced 'Dark Taxa' and quantify open-set false-positive over-classification!
    """
    rng = random.Random(seed)
    by_species: dict[str, list[dict[str, Any]]] = {}
    for r in records:
        by_species.setdefault(r["true_species"], []).append(r)

    species_list = sorted(by_species.keys())
    rng.shuffle(species_list)
    n_dark = max(1, int(round(len(species_list) * holdout_dark_taxa_frac))) if len(species_list) >= 6 else 0
    dark_species = set(species_list[:n_dark])

    ref_records: list[dict[str, Any]] = []
    query_templates: list[dict[str, Any]] = []

    for sp in species_list:
        recs = list(by_species[sp])
        rng.shuffle(recs)
        if sp in dark_species:
            # Withheld completely from Reference DB -> true label in closed-set species space is 'UNASSIGNED_DARK_TAXON'
            query_templates.extend({**x, "is_dark_taxon": True} for x in recs[:2])
        elif len(recs) >= 2:
            n_test = max(1, len(recs) // 5)
            for x in recs[:n_test]:
                query_templates.append({**x, "is_dark_taxon": False})
            ref_records.extend(recs[n_test:])
        else:
            # Single-record species stays in reference; also generate mutated amplicon reads only if needed
            ref_records.extend(recs)

    return ref_records, query_templates, dark_species


class KmerNBClassifier:
    """QIIME2 / RDP-style 8-mer Multinomial Naive Bayes taxonomic classifier with open-set confidence cutoff."""

    def __init__(self, k: int = 8, confidence_threshold: float = 0.70) -> None:
        self.k = k
        self.confidence_threshold = confidence_threshold
        self.vectorizer = CountVectorizer(analyzer="char", ngram_range=(k, k), lowercase=False)
        self.model = MultinomialNB(alpha=0.01)
        self.classes_: np.ndarray = np.array([])

    def fit(self, sequences: list[str], labels: list[str]) -> None:
        X = self.vectorizer.fit_transform(sequences)
        self.model.fit(X, labels)
        self.classes_ = self.model.classes_

    def predict(self, sequences: list[str]) -> list[str]:
        if not sequences:
            return []
        X = self.vectorizer.transform(sequences)
        probs = self.model.predict_proba(X)
        preds: list[str] = []
        for row in probs:
            best_idx = int(np.argmax(row))
            if float(row[best_idx]) >= self.confidence_threshold:
                preds.append(str(self.classes_[best_idx]))
            else:
                preds.append("UNASSIGNED")
        return preds


class SintaxBootstrapClassifier:
    """Vectorized SINTAX-style k-mer word-subset bootstrap taxonomic classifier (Edgar, 2016)."""

    def __init__(self, k: int = 8, n_boot: int = 32, subsample: int = 24, cutoff: float = 0.80, seed: int = 42) -> None:
        self.k = k
        self.n_boot = n_boot
        self.subsample = subsample
        self.cutoff = cutoff
        self.seed = seed
        self.vectorizer = CountVectorizer(analyzer="char", ngram_range=(k, k), lowercase=False, binary=True)
        self.ref_matrix_T: Any = None  # shape (n_vocab, n_refs) CSR
        self.ref_labels: list[str] = []

    def fit(self, sequences: list[str], labels: list[str]) -> None:
        self.ref_labels = list(labels)
        if sequences:
            X = self.vectorizer.fit_transform([s.upper() for s in sequences])
            self.ref_matrix_T = X.T.tocsr()

    def predict(self, sequences: list[str]) -> list[str]:
        if not sequences or self.ref_matrix_T is None or not self.ref_labels:
            return ["UNASSIGNED"] * len(sequences)
        Q = self.vectorizer.transform([s.upper() for s in sequences]).tocsr()
        rng = np.random.default_rng(self.seed)
        preds: list[str] = []

        for row_i in range(Q.shape[0]):
            word_indices = Q.indices[Q.indptr[row_i] : Q.indptr[row_i + 1]]
            if len(word_indices) < 8:
                preds.append("UNASSIGNED")
                continue
            sub_n = min(self.subsample, len(word_indices))
            # Sub-slice only the query's present k-mers -> dense (len(word_indices), n_refs) uint8/float32
            sub_ref = self.ref_matrix_T[word_indices, :].toarray()  # shape (W, n_refs)
            # Sample (n_boot, sub_n) indices without replacement per bootstrap replicate
            boot_hits = np.empty((self.n_boot, sub_ref.shape[1]), dtype=np.int16)
            for b in range(self.n_boot):
                chosen = rng.choice(len(word_indices), size=sub_n, replace=False)
                boot_hits[b] = sub_ref[chosen].sum(axis=0)

            best_refs = np.argmax(boot_hits, axis=1)
            best_vals = boot_hits[np.arange(self.n_boot), best_refs]
            boot_votes: Counter[str] = Counter()
            for r_idx, hits in zip(best_refs.tolist(), best_vals.tolist()):
                if hits / sub_n >= 0.45:
                    boot_votes[self.ref_labels[r_idx]] += 1
                else:
                    boot_votes["UNASSIGNED"] += 1
            top_label, top_count = boot_votes.most_common(1)[0]
            if top_label != "UNASSIGNED" and (top_count / self.n_boot) >= self.cutoff:
                preds.append(top_label)
            else:
                preds.append("UNASSIGNED")
        return preds


def run_blast_lca_classifier(
    ref_sequences: list[str],
    ref_labels: list[str],
    query_sequences: list[str],
    threads: int = 4,
    min_species_pident: float = 97.5,
    lca_margin: float = 0.5,
) -> list[str]:
    """Run local NCBI blastn (-num_threads 4) with Top-Hit / 0.5%-margin Species LCA assignment."""
    threads = min(max(int(threads), 1), 4)
    if not ref_sequences or not query_sequences or shutil.which("makeblastdb") is None or shutil.which("blastn") is None:
        return ["UNASSIGNED"] * len(query_sequences)

    with tempfile.TemporaryDirectory(prefix="bfr_blast_") as tmpdir:
        tmp = Path(tmpdir)
        ref_fa = tmp / "ref.fasta"
        qry_fa = tmp / "qry.fasta"
        db_prefix = tmp / "refdb"
        out_tsv = tmp / "hits.tsv"

        id_to_label: dict[str, str] = {}
        with ref_fa.open("w", encoding="utf-8") as fh:
            for i, (seq, lbl) in enumerate(zip(ref_sequences, ref_labels)):
                rid = f"REF_{i}"
                id_to_label[rid] = lbl
                fh.write(f">{rid}\n{seq}\n")

        with qry_fa.open("w", encoding="utf-8") as fh:
            for j, qseq in enumerate(query_sequences):
                fh.write(f">QRY_{j}\n{qseq}\n")

        subprocess.run(
            ["makeblastdb", "-in", str(ref_fa), "-dbtype", "nucl", "-out", str(db_prefix)],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        subprocess.run(
            [
                "blastn",
                "-query",
                str(qry_fa),
                "-db",
                str(db_prefix),
                "-num_threads",
                str(threads),
                "-max_target_seqs",
                "10",
                "-perc_identity",
                "90",
                "-outfmt",
                "6 qseqid sseqid pident bitscore qcovs",
                "-out",
                str(out_tsv),
            ],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

        hits_by_q: dict[str, list[tuple[str, float, float, float]]] = {}
        if out_tsv.exists():
            for line in out_tsv.read_text(encoding="utf-8").splitlines():
                parts = line.strip().split("\t")
                if len(parts) < 5:
                    continue
                qid, sid = parts[0], parts[1]
                pident, bits, qcov = float(parts[2]), float(parts[3]), float(parts[4])
                if qcov >= 80.0:
                    hits_by_q.setdefault(qid, []).append((sid, pident, bits, qcov))

        preds: list[str] = []
        for j in range(len(query_sequences)):
            qid = f"QRY_{j}"
            qhits = hits_by_q.get(qid, [])
            if not qhits:
                preds.append("UNASSIGNED")
                continue
            qhits.sort(key=lambda x: (x[2], x[1]), reverse=True)
            top_sid, top_pident, top_bits, _ = qhits[0]
            if top_pident < min_species_pident:
                preds.append("UNASSIGNED")
                continue
            candidate_species = {
                id_to_label[sid]
                for sid, pid, bits, _ in qhits
                if pid >= top_pident - lca_margin and bits >= top_bits * 0.98
            }
            if len(candidate_species) == 1:
                preds.append(next(iter(candidate_species)))
            else:
                # Ambiguous across multiple species -> conservative LCA drops species call
                preds.append("UNASSIGNED")
        return preds


def compute_metrics_with_bootstrap_ci(
    y_true: list[str],
    y_pred: list[str],
    is_dark_flags: list[bool],
    n_boot: int = 1000,
    seed: int = 2026,
) -> dict[str, Any]:
    """Compute macro Precision, Recall, F1, False-Positive Exotic/Dark Rate, and 95% Bootstrap CIs."""
    n = len(y_true)
    if n == 0:
        return {
            "precision_macro": 0.0,
            "recall_macro": 0.0,
            "f1_macro": 0.0,
            "f1_ci_low": 0.0,
            "f1_ci_high": 0.0,
            "false_positive_exotic_rate": 0.0,
        }

    def _calc(indices: list[int]) -> tuple[float, float, float, float]:
        yt = [y_true[i] for i in indices]
        yp = [y_pred[i] for i in indices]
        dk = [is_dark_flags[i] for i in indices]

        # Evaluate species-level macro precision/recall across true in-DB species
        true_species_set = sorted({t for t, d in zip(yt, dk) if not d and t != "UNASSIGNED"})
        precisions: list[float] = []
        recalls: list[float] = []
        f1s: list[float] = []

        for sp in true_species_set:
            tp = sum(1 for t, p in zip(yt, yp) if t == sp and p == sp)
            fp = sum(1 for t, p in zip(yt, yp) if t != sp and p == sp)
            fn = sum(1 for t, p in zip(yt, yp) if t == sp and p != sp)
            p_val = tp / (tp + fp) if (tp + fp) > 0 else 0.0
            r_val = tp / (tp + fn) if (tp + fn) > 0 else 0.0
            f_val = (2 * p_val * r_val / (p_val + r_val)) if (p_val + r_val) > 0 else 0.0
            precisions.append(p_val)
            recalls.append(r_val)
            f1s.append(f_val)

        p_macro = float(np.mean(precisions)) if precisions else 0.0
        r_macro = float(np.mean(recalls)) if recalls else 0.0
        f_macro = float(np.mean(f1s)) if f1s else 0.0

        # False-positive misassignment / exotic call rate across ALL queries (wrong non-UNASSIGNED species call)
        wrong_species_calls = sum(1 for t, p in zip(yt, yp) if p != "UNASSIGNED" and p != t)
        fp_rate = wrong_species_calls / len(yt)
        return p_macro, r_macro, f_macro, fp_rate

    full_idx = list(range(n))
    p_m, r_m, f_m, fp_r = _calc(full_idx)

    rng = np.random.default_rng(seed)
    boot_f1 = np.empty(n_boot, dtype=float)
    for b in range(n_boot):
        sample_idx = rng.integers(0, n, size=n).tolist()
        _, _, bf, _ = _calc(sample_idx)
        boot_f1[b] = bf

    ci_low, ci_high = np.percentile(boot_f1, [2.5, 97.5])
    return {
        "precision_macro": round(p_m, 4),
        "recall_macro": round(r_m, 4),
        "f1_macro": round(f_m, 4),
        "f1_ci_low": round(float(ci_low), 4),
        "f1_ci_high": round(float(ci_high), 4),
        "false_positive_exotic_rate": round(fp_r, 4),
    }


def evaluate_ecosystem_benchmark(
    ref_records: list[dict[str, Any]],
    query_templates: list[dict[str, Any]],
    reads_per_template: int = 5,
    error_rate: float = 0.005,
    threads: int = 4,
    seed: int = 20261116,
) -> list[dict[str, Any]]:
    """Evaluate all 4 classifiers on a Reference condition vs. simulated noisy query reads."""
    rng = random.Random(seed)
    ref_seqs = [r["sequence"] for r in ref_records]
    # Note: In Raw Unfiltered GenBank condition, reference labels carry 'reported_label' (which may be mislabeled!),
    # whereas ground-truth evaluation always tests against 'true_species'!
    ref_labels = [r["reported_label"] for r in ref_records]

    q_seqs: list[str] = []
    y_true: list[str] = []
    dark_flags: list[bool] = []
    for qt in query_templates:
        for _ in range(reads_per_template):
            q_seqs.append(mutate_sequence_reads(qt["sequence"], error_rate=error_rate, rng=rng))
            is_dk = bool(qt.get("is_dark_taxon", False))
            dark_flags.append(is_dk)
            y_true.append("UNASSIGNED" if is_dk else qt["true_species"])

    results: list[dict[str, Any]] = []
    if not ref_seqs or not q_seqs:
        return results

    # 1. Naive Majority Prior Baseline
    t0 = time.perf_counter()
    most_common_lbl = Counter(ref_labels).most_common(1)[0][0]
    pred_prior = [most_common_lbl] * len(q_seqs)
    lat_prior = 1000.0 * (time.perf_counter() - t0) / len(q_seqs)
    results.append(
        {
            "classifier": "Naive_Prior_Baseline",
            "n_queries": len(q_seqs),
            "latency_ms_per_query": round(lat_prior, 4),
            **compute_metrics_with_bootstrap_ci(y_true, pred_prior, dark_flags),
        }
    )

    # 2. 8-mer Multinomial Naive Bayes (QIIME2/RDP style)
    if len(set(ref_labels)) >= 2:
        nb = KmerNBClassifier(k=8, confidence_threshold=0.70)
        nb.fit(ref_seqs, ref_labels)
        t0 = time.perf_counter()
        pred_nb = nb.predict(q_seqs)
        lat_nb = 1000.0 * (time.perf_counter() - t0) / len(q_seqs)
        results.append(
            {
                "classifier": "Kmer_MultinomialNB",
                "n_queries": len(q_seqs),
                "latency_ms_per_query": round(lat_nb, 4),
                **compute_metrics_with_bootstrap_ci(y_true, pred_nb, dark_flags),
            }
        )

    # 3. SINTAX-style k-mer Bootstrap Classifier
    sintax = SintaxBootstrapClassifier(k=8, n_boot=32, subsample=24, cutoff=0.80, seed=seed)
    sintax.fit(ref_seqs, ref_labels)
    t0 = time.perf_counter()
    pred_sintax = sintax.predict(q_seqs)
    lat_sintax = 1000.0 * (time.perf_counter() - t0) / len(q_seqs)
    results.append(
        {
            "classifier": "SINTAX_Kmer_80",
            "n_queries": len(q_seqs),
            "latency_ms_per_query": round(lat_sintax, 4),
            **compute_metrics_with_bootstrap_ci(y_true, pred_sintax, dark_flags),
        }
    )

    # 4. BLAST+ Top-Hit / Species LCA
    t0 = time.perf_counter()
    pred_blast = run_blast_lca_classifier(ref_seqs, ref_labels, q_seqs, threads=threads)
    lat_blast = 1000.0 * (time.perf_counter() - t0) / len(q_seqs)
    results.append(
        {
            "classifier": "BLAST_LCA_97.5",
            "n_queries": len(q_seqs),
            "latency_ms_per_query": round(lat_blast, 4),
            **compute_metrics_with_bootstrap_ci(y_true, pred_blast, dark_flags),
        }
    )

    return results
