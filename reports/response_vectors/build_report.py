"""Build the initial companion report for the BulkFormer document.

Requires python-docx and matplotlib. This is a methods/implementation report;
it does not present the illustrative figure as a benchmark result.
Run: uv run --with python-docx python reports/response_vectors/build_report.py --out-dir results/report
"""
from __future__ import annotations

import argparse
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from docx import Document
from docx.shared import Inches, Pt, RGBColor
from docx.oxml import OxmlElement
from docx.oxml.ns import qn


def schematic(path):
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 2.65))
    genes = ["A", "B", "C", "D", "E", "F"]
    values = ([1.4, 0.2, 0.8, -0.3, np.nan, np.nan], [80, 50, 20, np.nan, 5, np.nan])
    colors = ["#4C78A8", "#DD8452"]
    for ax, vals, color in zip(axes, values, colors):
        vals = np.array(vals)
        good = np.isfinite(vals)
        ax.bar(np.arange(6)[good], vals[good], color=color, width=0.62)
        ax.set_xticks(range(6), genes)
        ax.axhline(0, color="#444444", linewidth=0.7)
        for i in np.where(~good)[0]:
            ax.text(i, 0.05, "NA", transform=ax.get_xaxis_transform(), ha="center", color="#666666", fontsize=10)
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_xlabel("Gene")
    axes[0].set_title("Observed response", fontsize=12, fontweight="bold")
    axes[0].set_ylabel("Oriented shrunken log2FC")
    axes[1].set_title("Predicted targeting", fontsize=12, fontweight="bold")
    axes[1].set_ylabel("Score (higher = stronger)")
    fig.tight_layout(w_pad=3)
    fig.savefig(path, dpi=220, facecolor="white")
    plt.close(fig)


def build(out: pathlib.Path):
    out.mkdir(parents=True, exist_ok=True)
    image = out / "vector_example.png"
    schematic(image)
    doc = Document()
    section = doc.sections[0]
    section.page_width, section.page_height = Inches(8.5), Inches(11)
    section.top_margin = section.bottom_margin = Inches(1)
    section.left_margin = section.right_margin = Inches(0.866)
    for name, size in (("Normal", 11), ("Title", 16), ("Heading 1", 13), ("Heading 2", 11)):
        style = doc.styles[name]
        style.font.name = "Arial"
        style.font.size = Pt(size)
        style.font.color.rgb = RGBColor(0, 0, 0)
        style.paragraph_format.space_after = Pt(7)
        style.paragraph_format.line_spacing = 1.08
        if name != "Normal":
            style.font.bold = True
            style.paragraph_format.keep_with_next = True
    doc.styles["Title"].paragraph_format.space_after = Pt(9)
    doc.styles["Heading 1"].paragraph_format.space_before = Pt(9)
    for style in doc.styles:
        for border in list(style.element.iter(qn("w:pBdr"))):
            border.getparent().remove(border)

    def p(text, style=None):
        return doc.add_paragraph(text, style)

    def h(text):
        doc.add_heading(text, level=1)

    def caption(number, text):
        para = p("")
        para.add_run(f"Figure {number}. ").bold = True
        para.add_run(text)

    doc.add_heading("Gene level response vectors and miRNA specificity", level=0)
    p("Version of 2 November 2026. Companion to the expression-space and BulkFormer report. "
      "Methods and first implementation; the biological analysis has not yet been run.")
    h("Why keep the genes")
    p("The earlier report found that the direction of an experiment's whole-transcriptome movement "
      "can be poorly reproducible even when the perturbed miRNA's seed family can be recovered. "
      "This motivates comparing predictions with the response of each gene, keeping gene identity "
      "instead of averaging the transcriptome into one sample embedding.")
    h("What we have implemented")
    p("We have implemented a separate analysis that reads the new benchmark tables and predictor "
      "scores. Its initial scope is the 21 OE test comparisons and the fixed list of 19,226 genes. "
      "Flagged comparisons are retained. KO and set-aside comparisons are skipped.")
    p("Each experiment has an observed response vector: one oriented shrunken DESeq2 log2FC per "
      "gene. For OE, we reverse the sign, so a positive value means repression. Each predictor "
      "has a targeting vector for the same miRNA, in the same gene order; higher scores mean "
      "stronger predicted targeting. Scores are matched by gene ID, not by their input row order.")
    para = p("")
    para.paragraph_format.keep_with_next = True
    para.add_run().add_picture(str(image), width=Inches(6.65))
    caption(1, "Illustrative values, not experimental results. Genes A-D have usable responses. "
            "D has no prediction but stays in the comparison at the bottom of the ranking. E has "
            "a prediction but no usable response, so it cannot contribute to alignment. F has neither. "
            "The two panels use different units.")

    doc.add_page_break()
    h("How we compare the vectors")
    p("For experiment e and gene g, the observed value is y(e,g) = -log2FC_shrunk(e,g). "
      "This is the fitted treatment effect, including the batch term when the experiment's design "
      "has one. It is different from the earlier descriptive movement, which subtracts the mean "
      "log-transformed TPM of controls from that of treated samples.")
    p("Most targeting scores do not predict a fold change in expression. We therefore compare "
      "rankings, rather than treating the score as a predicted log2FC. A positive Spearman "
      "correlation means that genes predicted to be stronger targets tend to be more repressed.")
    p("We also export centered, unit-length rank vectors. Their dot product equals Spearman "
      "correlation on the same genes. This gives the vector comparison a precise geometric "
      "interpretation without introducing a new embedding or fitting a new model.")
    h("Missing values and coverage")
    p("A missing or ineligible response stays missing; it is not a zero effect. All unscored "
      "genes form one tie below the scored genes, following the benchmark's full-pair-set policy. "
      "The exported files distinguish response usability from predictor coverage. Equal scores "
      "receive average ranks, so gene names do not break ties.")
    p("Alignment is unavailable when fewer than three genes have usable responses, the prediction "
      "ranking is constant, or the observed response is constant. These cases are reported as NA "
      "with a reason; an unavailable alignment is not evidence of zero biological association.")
    h("Is the agreement specific to the perturbed miRNA")
    p("For the same observed response, we calculate alignment with the predictor's scores for "
      "each miRNA in the supplied registry. The perturbed miRNA is the cognate. Registry miRNAs "
      "with a different seed at nucleotides 2-8 are its decoys; same-seed members are excluded.")
    p("We report the cognate's descriptive rank against usable decoys, the fraction of decoys "
      "with better alignment, and its z-score relative to their mean and sample SD. Ties receive "
      "half weight: rank = 1 + better decoys + half of tied decoys. A rank of 1 means no decoy "
      "does as well. The z-score is NA with fewer than two usable decoys or no decoy variation.")
    h("Accounting for control expression")
    p("We repeat alignment within tenths of control expression, using the benchmark's existing "
      "definition. Bins need at least 20 usable genes and a defined correlation; correlations are "
      "averaged with weights equal to their gene counts. Specificity is then recalculated from "
      "these alignments. This checks whether agreement survives comparison among genes of "
      "similar baseline expression; it does not establish that all expression-related confounding "
      "has been removed.")

    doc.add_page_break()
    h("What can be concluded so far")
    p("The vector construction and specificity calculations are implemented. All 40 Python tests "
      "passed at this version, including a synthetic end-to-end run and checks of gene alignment, "
      "missing values, same-seed exclusion and ties. The exported rank-vector dot product was "
      "checked against Spearman correlation.")
    p("These are software checks, not biological results. We have not yet run this command on "
      "the full benchmark inputs. There is therefore no new evidence here that one predictor "
      "outperforms another, or that the vector analysis improves on the existing benchmark.")
    h("The next analysis")
    p("The first biological run should show alignment and specificity together, per experiment "
      "and predictor. An experiment with positive alignment but many better decoys has general "
      "agreement without strong miRNA specificity. A cognate that ranks well despite weak "
      "overall alignment may carry a small specific signature within a much larger response.")
    p("The next results figure will compare overall and within-expression-tenth alignment, with "
      "flagged experiments identified. A second figure will place cognate alignment against "
      "the fraction of better decoys, reporting the usable decoy count and coverage alongside "
      "the comparison. Numerical values for every experiment will remain available in the tables.")
    h("What this does not yet measure")
    p("Specificity is relative to the available registry and predictor coverage. It is not a "
      "calibrated p-value or a rank among all human seed families. Different-seed families may "
      "have several registry members, so decoys are not independent family-level observations. "
      "Comparisons between tools must account for their different usable candidate panels.")
    p("We have not yet added replicate stability, controls-only contrasts or a seed-count baseline. "
      "Stability needs sample-level inputs or refitted effects from separate subsets; controls-only "
      "contrasts must preserve batch and donor structure. A seed baseline needs the matching "
      "3'UTR site information. These checks are needed before interpreting alignment as a robust "
      "miRNA-specific effect. Response vectors may include downstream and technical effects; "
      "agreement does not establish direct binding.")
    h("Implementation and reproducibility")
    p("The analysis starts from PR #141 at commit 39be410 and is implemented on "
      "feature/response-vector-analysis at commit 2cca628. Run: "
      "uv run funmirbench-vectors --config benchmark.yaml. The normal benchmark is unchanged.")
    p("Each run saves the vector files, candidate alignments, specificity table, selected and "
      "skipped comparisons, registries and configuration under results/<date_time>/response_vectors/. "
      "This companion report will be updated with the biological run's figures and findings.")

    footer = section.footer.paragraphs[0]
    footer.alignment = 2
    footer.add_run("Response vectors | ").font.size = Pt(9)
    fld = OxmlElement("w:fldSimple")
    fld.set(qn("w:instr"), "PAGE")
    footer._p.append(fld)
    doc.core_properties.title = "Gene level response vectors and miRNA specificity"
    doc.core_properties.subject = "Companion report for expression-space and BulkFormer analyses"
    for para in doc.paragraphs:
        for border in list(para._p.iter(qn("w:pBdr"))):
            border.getparent().remove(border)
    target = out / "Response_vectors_and_miRNA_specificity.docx"
    doc.save(target)
    print(target)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=pathlib.Path, required=True)
    build(parser.parse_args().out_dir)
