"""Export publication-formatted Wiley ScholarOne Word (.docx) manuscript & cover letter with line numbers."""

from __future__ import annotations

import json
from pathlib import Path
import re
import sqlite3
from typing import Any

from docx import Document
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, RGBColor

WORKSPACE = Path("/home/sabuj/Research/01_Active_Projects/Antigravity")
DB_PATH = WORKSPACE / "02_Data" / "processed" / "bengalfishref.sqlite"
DRAFT_DIR = WORKSPACE / "04_Manuscripts" / "drafts"
FIG_DIR = WORKSPACE / "04_Manuscripts" / "figures"


def _enable_continuous_line_numbers_and_page_footer(doc: Document) -> None:
    """Configure 1-inch margins, continuous ScholarOne line numbering, and centered footer page numbers."""
    for section in doc.sections:
        section.top_margin = Inches(1.0)
        section.bottom_margin = Inches(1.0)
        section.left_margin = Inches(1.0)
        section.right_margin = Inches(1.0)

        sectPr = section._sectPr
        lnNumType = OxmlElement("w:lnNumType")
        lnNumType.set(qn("w:countBy"), "1")
        lnNumType.set(qn("w:restart"), "continuous")
        sectPr.append(lnNumType)

        footer = section.footer
        p = footer.paragraphs[0]
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        run = p.add_run("Page ")
        run.font.name = "Times New Roman"
        run.font.size = Pt(10)
        fldSimple = OxmlElement("w:fldSimple")
        fldSimple.set(qn("w:instr"), "PAGE")
        p._p.append(fldSimple)


def _set_base_styles(doc: Document, double_spaced: bool = True) -> None:
    normal = doc.styles["Normal"]
    normal.font.name = "Times New Roman"
    normal.font.size = Pt(12)
    normal.font.color.rgb = RGBColor(0, 0, 0)
    pf = normal.paragraph_format
    pf.line_spacing = 2.0 if double_spaced else 1.15
    pf.space_after = Pt(6 if double_spaced else 8)


def _add_formatted_runs(paragraph: Any, text: str, base_size_pt: float = 12.0) -> None:
    """Render inline Markdown (**bold**, *italic*, `code`, $^{sup}$) cleanly into Word runs."""
    # Clean display math wrappers or LaTeX macros lightly for readable Word math strings
    cleaned = (
        text.replace("\\$", "$")
        .replace("\\cdot", "×")
        .replace("\\le", "≤")
        .replace("\\ge", "≥")
        .replace("\\times", "×")
        .replace("\\text{BDTPI}", "BDTPI")
    )
    # Tokenize **bold**, *italic*, `code`, and $^{...}$
    token_re = re.compile(r"(\*\*[^*]+\*\*|\*[^*]+\*|`[^`]+`|\$\^\{[^}]+\}\$)")
    pos = 0
    for m in token_re.finditer(cleaned):
        if m.start() > pos:
            r = paragraph.add_run(cleaned[pos : m.start()])
            r.font.name = "Times New Roman"
            r.font.size = Pt(base_size_pt)
        tok = m.group(0)
        if tok.startswith("**") and tok.endswith("**"):
            inner = tok[2:-2]
            # Check nested italics inside bold (e.g., ***Species*** or **text *Species***)
            sub_parts = re.split(r"(\*[^*]+\*)", inner)
            for sp in sub_parts:
                if not sp:
                    continue
                r = paragraph.add_run(sp[1:-1] if (sp.startswith("*") and sp.endswith("*")) else sp)
                r.bold = True
                if sp.startswith("*") and sp.endswith("*"):
                    r.italic = True
                r.font.name = "Times New Roman"
                r.font.size = Pt(base_size_pt)
        elif tok.startswith("*") and tok.endswith("*"):
            r = paragraph.add_run(tok[1:-1])
            r.italic = True
            r.font.name = "Times New Roman"
            r.font.size = Pt(base_size_pt)
        elif tok.startswith("`") and tok.endswith("`"):
            r = paragraph.add_run(tok[1:-1])
            r.font.name = "Courier New"
            r.font.size = Pt(max(base_size_pt - 1.5, 9.5))
        elif tok.startswith("$^{") and tok.endswith("}$"):
            r = paragraph.add_run(tok[3:-2])
            r.font.superscript = True
            r.font.name = "Times New Roman"
            r.font.size = Pt(base_size_pt)
        pos = m.end()

    if pos < len(cleaned):
        r = paragraph.add_run(cleaned[pos:])
        r.font.name = "Times New Roman"
        r.font.size = Pt(base_size_pt)


def _add_styled_table(doc: Document, caption: str, headers: list[str], rows: list[list[str]]) -> None:
    cap_p = doc.add_paragraph()
    cap_p.paragraph_format.line_spacing = 1.15
    cap_p.paragraph_format.space_before = Pt(12)
    cap_p.paragraph_format.space_after = Pt(4)
    _add_formatted_runs(cap_p, caption, base_size_pt=11.0)

    tbl = doc.add_table(rows=1 + len(rows), cols=len(headers))
    tbl.style = "Table Grid"
    tbl.alignment = WD_TABLE_ALIGNMENT.CENTER

    hdr_cells = tbl.rows[0].cells
    for idx, h in enumerate(headers):
        p = hdr_cells[idx].paragraphs[0]
        p.paragraph_format.line_spacing = 1.05
        p.paragraph_format.space_after = Pt(2)
        r = p.add_run(h)
        r.bold = True
        r.font.name = "Times New Roman"
        r.font.size = Pt(9.5)

    for r_i, row_vals in enumerate(rows):
        cells = tbl.rows[r_i + 1].cells
        for c_i, val in enumerate(row_vals):
            p = cells[c_i].paragraphs[0]
            p.paragraph_format.line_spacing = 1.05
            p.paragraph_format.space_after = Pt(2)
            _add_formatted_runs(p, str(val), base_size_pt=9.5)


def _render_markdown_into_docx(doc: Document, md_text: str) -> None:
    lines = md_text.splitlines()
    i = 0
    while i < len(lines):
        raw = lines[i].rstrip()
        stripped = raw.strip()
        i += 1
        if not stripped or stripped == "---":
            continue
        if stripped.startswith("# "):
            p = doc.add_paragraph()
            p.paragraph_format.line_spacing = 1.3
            p.paragraph_format.space_after = Pt(10)
            r = p.add_run(stripped[2:].strip())
            r.bold = True
            r.font.name = "Times New Roman"
            r.font.size = Pt(15)
        elif stripped.startswith("## "):
            p = doc.add_paragraph()
            p.paragraph_format.space_before = Pt(12)
            p.paragraph_format.space_after = Pt(4)
            r = p.add_run(stripped[3:].strip())
            r.bold = True
            r.font.name = "Times New Roman"
            r.font.size = Pt(13)
        elif stripped.startswith("### "):
            p = doc.add_paragraph()
            p.paragraph_format.space_before = Pt(8)
            p.paragraph_format.space_after = Pt(2)
            r = p.add_run(stripped[4:].strip())
            r.bold = True
            r.italic = True
            r.font.name = "Times New Roman"
            r.font.size = Pt(12)
        elif stripped.startswith("$$") and stripped.endswith("$$"):
            eq_txt = (
                stripped[2:-2]
                .replace("\\text{BDTPI}", "BDTPI")
                .replace("\\cdot", " × ")
                .replace("S_{\\text{cons}}", "S_cons")
                .replace("S_{\\text{comm}}", "S_comm")
                .replace("\\sum_{m \\in \\{\\text{COI}, \\text{12S}, \\text{16S}, \\text{Mitogenome}\\}} w_m D_m", "Σ (w_m × D_m)")
            )
            p = doc.add_paragraph()
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER
            r = p.add_run(eq_txt)
            r.italic = True
            r.font.name = "Times New Roman"
            r.font.size = Pt(11.5)
        elif stripped.startswith("* ") or stripped.startswith("- "):
            p = doc.add_paragraph(style="List Bullet")
            _add_formatted_runs(p, stripped[2:].strip(), base_size_pt=12.0)
        else:
            p = doc.add_paragraph()
            _add_formatted_runs(p, stripped, base_size_pt=12.0)


def build_submission_docx_files(db_path: Path = DB_PATH) -> dict[str, str]:
    """Generate ScholarOne-ready Manuscript (.docx) and Cover Letter (.docx)."""
    conn = sqlite3.connect(str(db_path))
    conn.row_factory = sqlite3.Row

    # 1. Manuscript .docx (double-spaced, continuous line numbers, tables + figures)
    ms_md = (DRAFT_DIR / "BengalFishRef_MER_Manuscript.md").read_text(encoding="utf-8")
    ms_doc = Document()
    _set_base_styles(ms_doc, double_spaced=True)
    _enable_continuous_line_numbers_and_page_footer(ms_doc)
    _render_markdown_into_docx(ms_doc, ms_md)

    ms_doc.add_page_break()
    h_tbl = ms_doc.add_paragraph()
    r_ht = h_tbl.add_run("Tables")
    r_ht.bold = True
    r_ht.font.name = "Times New Roman"
    r_ht.font.size = Pt(14)

    # Table 1 from SQLite
    t1_rows_db = conn.execute(
        """
        SELECT marker,
               SUM(CASE WHEN quality_tier='Gold' THEN 1 ELSE 0 END) AS gold,
               SUM(CASE WHEN quality_tier='Silver' THEN 1 ELSE 0 END) AS silver,
               SUM(CASE WHEN quality_tier='Quarantined' THEN 1 ELSE 0 END) AS quar,
               SUM(has_stop_codon_numt) AS numts,
               SUM(CASE WHEN phylo_concordance_flag='Mislabeled_Outlier' THEN 1 ELSE 0 END) AS mislbl,
               COUNT(*) AS tot,
               COUNT(DISTINCT CASE WHEN quality_tier IN ('Gold','Silver') THEN resolved_taxon_id END) AS n_sp
        FROM sequences
        GROUP BY marker
        ORDER BY tot DESC
        """
    ).fetchall()
    t1_formatted = [
        [
            f"**{r['marker']}**",
            f"{r['n_sp']}",
            f"{r['gold']:,}",
            f"{r['silver']:,}",
            f"{r['quar']:,}",
            f"{r['numts']:,}",
            f"{r['mislbl']:,}",
            f"**{r['tot']:,}**",
        ]
        for r in t1_rows_db
    ]
    _add_styled_table(
        ms_doc,
        "**Table 1.** Multi-locus mitochondrial sequence inventory and two-stage phylogenetic curation summary across 859 documented fish and commercial decapod species in `BengalFishRef v1.0`.",
        ["Locus", "Species Covered", "Gold", "Silver", "Quarantined", "NUMT Stops", "Mislabeled Outliers", "Total Audited"],
        t1_formatted,
    )

    # Table 2: Primer coverage summary
    t2_rows_db = conn.execute(
        """
        SELECT p.primer_name, p.marker,
               ROUND(100.0 * AVG(CASE WHEN t.tax_group='Actinopterygii' THEN p.amplifies_in_silico END), 1) AS act_pct,
               ROUND(100.0 * AVG(CASE WHEN t.tax_group='Elasmobranchii' THEN p.amplifies_in_silico END), 1) AS ela_pct,
               ROUND(100.0 * AVG(CASE WHEN t.tax_group='Crustacea' THEN p.amplifies_in_silico END), 1) AS cru_pct,
               ROUND(AVG(p.weighted_penalty), 2) AS mean_wmm,
               COUNT(*) AS n_eval
        FROM primer_evaluations p
        JOIN taxa t ON p.taxon_id = t.taxon_id
        GROUP BY p.primer_name, p.marker
        ORDER BY p.marker, p.primer_name
        """
    ).fetchall()
    t2_formatted = [
        [
            f"`{r['primer_name']}`",
            r["marker"],
            f"{r['act_pct']}%" if r["act_pct"] is not None else "—",
            f"{r['ela_pct']}%" if r["ela_pct"] is not None else "—",
            f"{r['cru_pct']}%" if r["cru_pct"] is not None else "—",
            f"{r['mean_wmm']:.2f}",
            f"{r['n_eval']:,}",
        ]
        for r in t2_rows_db
    ]
    _add_styled_table(
        ms_doc,
        "**Table 2.** *In silico* amplification success (%) and mean 3'-weighted mismatch penalty (*W_mm*) across 14 universal metabarcoding primer pairs evaluated on `BengalFishRef v1.0` templates (*n* = 25,300 tests).",
        ["Primer Pair", "Locus", "Actinopterygii", "Elasmobranchii", "Crustacea", "Mean W_mm", "Templates Tested"],
        t2_formatted,
    )

    # Table 3: Simulated eDNA Classifier Benchmark (0.5% error)
    t3_rows_db = conn.execute(
        """
        SELECT ecosystem, classifier, db_condition, f1_macro, f1_ci_low, f1_ci_high, false_positive_exotic_rate
        FROM benchmark_results
        WHERE error_rate = 0.005 AND classifier != 'Naive_Prior_Baseline'
        ORDER BY ecosystem, classifier, db_condition DESC
        """
    ).fetchall()
    t3_formatted = [
        [
            r["ecosystem"].replace("_", " "),
            r["classifier"],
            "BengalFishRef Curated" if "BengalFishRef" in r["db_condition"] else "Raw Public GenBank",
            f"{r['f1_macro']:.4f} ({r['f1_ci_low']:.4f}–{r['f1_ci_high']:.4f})",
            f"{100.0 * r['false_positive_exotic_rate']:.2f}%",
        ]
        for r in t3_rows_db
    ]
    _add_styled_table(
        ms_doc,
        "**Table 3.** Four-ecosystem simulated eDNA taxonomic classification benchmark at 0.5% read error (with 12% held-out regional dark taxa and 1,000-replicate 95% bootstrap confidence intervals).",
        ["Ecosystem Assemblage", "Classifier", "Reference Condition", "Macro-F1 (95% CI)", "False-Positive Rate"],
        t3_formatted,
    )

    # Figures section
    ms_doc.add_page_break()
    h_fig = ms_doc.add_paragraph()
    r_hf = h_fig.add_run("Figures and Figure Legends")
    r_hf.bold = True
    r_hf.font.name = "Times New Roman"
    r_hf.font.size = Pt(14)

    fig_specs = [
        (
            "Fig1_Taxonomic_Coverage_and_DarkTaxa_BDTPI.png",
            "**Figure 1.** **(A)** Species-level molecular reference coverage across five mitochondrial loci (`COI`, `12S`, `16S`, `Cytb`, and complete mitogenomes) for 859 documented Actinopterygii (*n* = 734), Elasmobranchii (*n* = 79), and commercial Crustacea (*n* = 46) of Bangladesh and the northern Bay of Bengal. **(B)** Highest-priority families (*n* ≥ 3 species) ranked by mean Bengal Dark-Taxa Priority Index (BDTPI) for targeted voucher sequencing.",
        ),
        (
            "Fig2_Sequence_Hygiene_and_Barcode_Gap.png",
            "**Figure 2.** **(A)** Three-tier curation outcomes (`Gold`, `Silver`, `Quarantined`) across 27,467 public mitochondrial accessions. **(B)** Kimura 2-parameter (K2P) intraspecific versus nearest heterospecific distance distribution across `COI` records, illustrating clear separation of concordant vouchers (*n* = 9,503), putative South Asian cryptic lineages (*n* = 2,087), and quarantined mislabeled GenBank outliers (*n* = 452) clustering near 0% heterospecific divergence.",
        ),
        (
            "Fig3_eDNA_Primer_Coverage_Heatmap.png",
            "**Figure 3.** *In silico* amplification coverage (%) of 14 universal fish and decapod metabarcoding primer sets across 25,300 template evaluations spanning Bangladeshi Actinopterygii, commercial Crustacea, and Elasmobranchii.",
        ),
        (
            "Fig4_Simulated_eDNA_Classifier_Benchmark.png",
            "**Figure 4.** Simulated multi-ecosystem eDNA taxonomic assignment benchmark (0.5% read error, 95% CIs across four Bangladeshi assemblages) comparing unfiltered public GenBank against curated `BengalFishRef v1.0`: **(A)** species-level Macro-F1 and **(B)** suppression of false-positive / mislabeled species calls.",
        ),
    ]
    for fname, legend in fig_specs:
        fpath = FIG_DIR / fname
        if fpath.exists():
            pic_p = ms_doc.add_paragraph()
            pic_p.alignment = WD_ALIGN_PARAGRAPH.CENTER
            pic_p.add_run().add_picture(str(fpath), width=Inches(6.2))
            leg_p = ms_doc.add_paragraph()
            leg_p.paragraph_format.line_spacing = 1.35
            leg_p.paragraph_format.space_after = Pt(14)
            _add_formatted_runs(leg_p, legend, base_size_pt=10.5)

    ms_out = DRAFT_DIR / "BengalFishRef_MER_Submission_Manuscript.docx"
    ms_doc.save(str(ms_out))

    # 2. Cover Letter .docx (clean single-spaced formal formatting)
    cl_md = (DRAFT_DIR / "Cover_Letter_Molecular_Ecology_Resources.md").read_text(encoding="utf-8")
    cl_doc = Document()
    _set_base_styles(cl_doc, double_spaced=False)
    for section in cl_doc.sections:
        section.top_margin = Inches(1.0)
        section.bottom_margin = Inches(1.0)
        section.left_margin = Inches(1.0)
        section.right_margin = Inches(1.0)
    _render_markdown_into_docx(cl_doc, cl_md)
    cl_out = DRAFT_DIR / "Cover_Letter_Molecular_Ecology_Resources.docx"
    cl_doc.save(str(cl_out))

    return {"manuscript_docx": str(ms_out), "cover_letter_docx": str(cl_out)}


if __name__ == "__main__":
    print(json.dumps(build_submission_docx_files(), indent=2))
