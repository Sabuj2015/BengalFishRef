"""Generate 600-DPI publication figures and summary tables for BengalFishRef (Molecular Ecology Resources)."""

from __future__ import annotations

import json
from pathlib import Path
import sqlite3

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from bengalfishref.taxonomy import generate_dark_taxa_gap_report

WORKSPACE = Path("/home/sabuj/Research/01_Active_Projects/Antigravity")
DB_PATH = WORKSPACE / "02_Data" / "processed" / "bengalfishref.sqlite"
FIG_DIR = WORKSPACE / "04_Manuscripts" / "figures"
TAB_DIR = WORKSPACE / "04_Manuscripts" / "tables"


def set_pub_style() -> None:
    sns.set_theme(style="whitegrid", context="paper", font_scale=1.18)
    plt.rcParams.update(
        {
            "figure.dpi": 300,
            "savefig.dpi": 600,
            "font.family": "sans-serif",
            "axes.edgecolor": "#334155",
            "axes.linewidth": 0.9,
        }
    )


def save_dual(fig: plt.Figure, stem: str) -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_DIR / f"{stem}.png", dpi=600, bbox_inches="tight")
    fig.savefig(FIG_DIR / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def generate_all_figures_and_tables(db_path: Path = DB_PATH) -> dict[str, str]:
    set_pub_style()
    TAB_DIR.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(str(db_path))
    conn.row_factory = sqlite3.Row

    gap_rows = generate_dark_taxa_gap_report(conn)
    df_gap = pd.DataFrame(gap_rows)
    df_gap.to_csv(TAB_DIR / "Table_S1_BengalFishRef_DarkTaxa_BDTPI.tsv", sep="\t", index=False)

    # Figure 1: Marker Coverage & Dark-Taxa Breakdown by Taxonomic Group and Family BDTPI
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.0))
    markers = ["COI", "12S", "16S", "Cytb", "Mitogenome"]
    cov_records = []
    for grp, gdf in df_gap.groupby("tax_group"):
        n_tot = max(len(gdf), 1)
        for m in markers:
            pct = 100.0 * (gdf[f"{m}_usable"] > 0).sum() / n_tot
            cov_records.append(
                {
                    "Taxonomic Clade": f"{grp} (n={n_tot})",
                    "Mitochondrial Locus": m,
                    "Species Coverage (%)": pct,
                }
            )
    df_cov = pd.DataFrame(cov_records)
    sns.barplot(
        data=df_cov,
        x="Mitochondrial Locus",
        y="Species Coverage (%)",
        hue="Taxonomic Clade",
        palette=["#1e3a8a", "#0f766e", "#16a34a"],
        ax=axes[0],
    )
    axes[0].set_ylim(0, 105)
    axes[0].set_title("A. Species-Level Reference Availability Across Loci", fontweight="bold")
    axes[0].legend(frameon=True, loc="upper right", fontsize=9)

    top_fam = (
        df_gap[~df_gap["family_name"].isin(["Animalia", "Incertae_sedis", "Actinopterygii", "Elasmobranchii"])]
        .groupby("family_name")
        .agg(mean_bdtpi=("bdtpi_score", "mean"), n_species=("taxon_id", "count"))
        .query("n_species >= 3")
        .sort_values("mean_bdtpi", ascending=False)
        .head(12)
        .reset_index()
    )
    sns.barplot(
        data=top_fam,
        y="family_name",
        x="mean_bdtpi",
        color="#0f766e",
        ax=axes[1],
    )
    axes[1].set_xlabel("Mean Bengal Dark-Taxa Priority Index (BDTPI)")
    axes[1].set_ylabel("Taxonomic Family (n >= 3 species)")
    axes[1].set_title("B. Top Priority Families for FBD-NIB Voucher Sampling", fontweight="bold")
    fig.tight_layout()
    save_dual(fig, "Fig1_Taxonomic_Coverage_and_DarkTaxa_BDTPI")

    # Figure 2: Sequence Curation Tiers & COI Barcode Gap / Mislabeling Space
    df_seq = pd.read_sql_query(
        """
        SELECT marker, quality_tier, phylo_concordance_flag, quarantine_reason,
               k2p_nearest_conspecific, k2p_nearest_heterospecific
        FROM sequences
        """,
        conn,
    )
    if not df_seq.empty:
        df_seq.to_csv(TAB_DIR / "Table_S2_Sequence_Audit_Summary.tsv", sep="\t", index=False)
        fig, axes = plt.subplots(1, 2, figsize=(12.4, 4.9))
        tier_ct = df_seq.groupby(["marker", "quality_tier"]).size().reset_index(name="Count")
        sns.barplot(
            data=tier_ct,
            x="marker",
            y="Count",
            hue="quality_tier",
            order=["COI", "Cytb", "16S", "12S", "Mitogenome"],
            hue_order=["Gold", "Silver", "Quarantined"],
            palette={"Gold": "#0f766e", "Silver": "#38bdf8", "Quarantined": "#dc2626"},
            ax=axes[0],
        )
        axes[0].set_xlabel("Mitochondrial Locus")
        axes[0].set_ylabel("Audited Public Accessions (n)")
        axes[0].set_title("A. Three-Tier Sequence Curation Outcome (n = 27,467)", fontweight="bold")
        axes[0].legend(title="Curation Tier", frameon=True)

        coi_bar = df_seq[
            (df_seq["marker"] == "COI")
            & df_seq["k2p_nearest_conspecific"].notna()
            & df_seq["k2p_nearest_heterospecific"].notna()
            & (df_seq["k2p_nearest_conspecific"] <= 0.32)
            & (df_seq["k2p_nearest_heterospecific"] <= 0.35)
        ].copy()
        if not coi_bar.empty:
            status_map = {
                "Monophyletic": "Concordant Voucher",
                "Paraphyletic_Shared": "Shared Haplotype",
                "Cryptic_Lineage_Candidate": "Putative Cryptic (>=3.5% K2P)",
                "Mislabeled_Outlier": "Quarantined Mislabeled Outlier",
            }
            coi_bar["Audit Class"] = coi_bar["phylo_concordance_flag"].map(status_map).fillna("Concordant Voucher")
            palette_audit = {
                "Concordant Voucher": "#0284c7",
                "Shared Haplotype": "#64748b",
                "Putative Cryptic (>=3.5% K2P)": "#d97706",
                "Quarantined Mislabeled Outlier": "#dc2626",
            }
            for cls_name in [
                "Concordant Voucher",
                "Shared Haplotype",
                "Putative Cryptic (>=3.5% K2P)",
                "Quarantined Mislabeled Outlier",
            ]:
                sub = coi_bar[coi_bar["Audit Class"] == cls_name]
                if sub.empty:
                    continue
                axes[1].scatter(
                    sub["k2p_nearest_conspecific"] * 100.0,
                    sub["k2p_nearest_heterospecific"] * 100.0,
                    alpha=0.72 if "Quarantined" in cls_name else 0.42,
                    c=palette_audit[cls_name],
                    s=24 if "Quarantined" in cls_name else 15,
                    edgecolors="none",
                    label=f"{cls_name} (n={len(sub)})",
                )
            axes[1].plot([0, 32], [0, 32], "--", color="#0f172a", lw=1.1, label="1:1 Barcode Gap Line")
            axes[1].set_xlim(-0.5, 32)
            axes[1].set_ylim(-0.5, 35)
            axes[1].set_xlabel("Nearest/Median Intraspecific COI K2P (%)")
            axes[1].set_ylabel("Nearest Heterospecific COI K2P (%)")
            axes[1].set_title("B. COI Barcode Gap & Mislabeled GenBank Outliers", fontweight="bold")
            axes[1].legend(frameon=True, fontsize=8.2, loc="upper right")
        fig.tight_layout()
        save_dual(fig, "Fig2_Sequence_Hygiene_and_Barcode_Gap")

    # Figure 3: 14-Primer In Silico Coverage Heatmap
    df_prm = pd.read_sql_query(
        """
        SELECT p.primer_name, p.marker, t.tax_group,
               100.0 * AVG(p.amplifies_in_silico) AS amp_pct,
               AVG(p.weighted_penalty) AS mean_penalty
        FROM primer_evaluations p
        JOIN taxa t ON p.taxon_id = t.taxon_id
        GROUP BY p.primer_name, p.marker, t.tax_group
        """,
        conn,
    )
    if not df_prm.empty:
        df_prm.to_csv(TAB_DIR / "Table_S3_Primer_Thermodynamic_Benchmark.tsv", sep="\t", index=False)
        piv = df_prm.pivot(index="primer_name", columns="tax_group", values="amp_pct").fillna(0.0)
        fig, ax = plt.subplots(figsize=(7.8, 6.4))
        sns.heatmap(
            piv,
            annot=True,
            fmt=".1f",
            cmap="YlGnBu",
            cbar_kws={"label": "In Silico Amplification Success (%)"},
            ax=ax,
        )
        ax.set_xlabel("Taxonomic Clade")
        ax.set_ylabel("Metabarcoding Primer Pair")
        ax.set_title("In Silico eDNA Primer Coverage Across Bangladesh Clades (25,300 Tests)", fontweight="bold")
        fig.tight_layout()
        save_dual(fig, "Fig3_eDNA_Primer_Coverage_Heatmap")

    # Figure 4: 4-Ecosystem Simulated eDNA Classifier Benchmark (Zoomed to Active Classifiers)
    df_bm_all = pd.read_sql_query("SELECT * FROM benchmark_results", conn)
    if not df_bm_all.empty:
        df_bm_all.to_csv(TAB_DIR / "Table_S4_Mock_Community_Benchmark_95CI.tsv", sep="\t", index=False)
        df_bm = df_bm_all[
            (df_bm_all["error_rate"] == 0.005) & (df_bm_all["classifier"] != "Naive_Prior_Baseline")
        ].copy()
        label_map = {
            "Kmer_MultinomialNB": "8-mer Naive Bayes\n(QIIME2/RDP)",
            "SINTAX_Kmer_80": "SINTAX 8-mer\n(80% Bootstrap)",
            "BLAST_LCA_97.5": "NCBI BLAST+\n(97.5% LCA)",
        }
        cond_map = {
            "Raw_Unfiltered_Public": "Unfiltered Public GenBank",
            "BengalFishRef_Curated": "BengalFishRef v1.0 (Curated)",
        }
        df_bm["Classifier"] = df_bm["classifier"].map(label_map)
        df_bm["Reference Database"] = df_bm["db_condition"].map(cond_map)
        df_bm["FP Rate (%)"] = df_bm["false_positive_exotic_rate"] * 100.0

        fig, axes = plt.subplots(1, 2, figsize=(12.6, 5.0))
        sns.barplot(
            data=df_bm,
            x="Classifier",
            y="f1_macro",
            hue="Reference Database",
            palette={"Unfiltered Public GenBank": "#94a3b8", "BengalFishRef v1.0 (Curated)": "#0f766e"},
            errorbar=("ci", 95),
            ax=axes[0],
        )
        axes[0].set_ylim(0.80, 1.00)
        axes[0].set_xlabel("")
        axes[0].set_ylabel("Species-Level Macro-F1 (0.5% Sequencing Error)")
        axes[0].set_title("A. Taxonomic Assignment Accuracy Across 4 Ecosystems", fontweight="bold")
        axes[0].legend(frameon=True, loc="lower right")

        sns.barplot(
            data=df_bm,
            x="Classifier",
            y="FP Rate (%)",
            hue="Reference Database",
            palette={"Unfiltered Public GenBank": "#dc2626", "BengalFishRef v1.0 (Curated)": "#16a34a"},
            errorbar=("ci", 95),
            ax=axes[1],
        )
        axes[1].set_ylim(0, 10.5)
        axes[1].set_xlabel("")
        axes[1].set_ylabel("False-Positive / Mislabeled Species Call Rate (%)")
        axes[1].set_title("B. Suppression of False-Positive Species Detections", fontweight="bold")
        axes[1].legend(frameon=True, loc="upper right")
        fig.tight_layout()
        save_dual(fig, "Fig4_Simulated_eDNA_Classifier_Benchmark")

    return {"figures_dir": str(FIG_DIR), "tables_dir": str(TAB_DIR)}


if __name__ == "__main__":
    print(json.dumps(generate_all_figures_and_tables(), indent=2))
