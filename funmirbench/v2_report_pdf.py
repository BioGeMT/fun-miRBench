"""REPORT.pdf of the v2 evaluation, with the pages and style of the run report (run_report_pdf.py)."""
from __future__ import annotations

import pathlib

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

from funmirbench.benchmark_reports import _format_summary_value
from funmirbench.run_report_pdf import (
    INVENTORY_WRAP_SCALE,
    REPORT_TABLE_SIZE,
    _bullet_block,
    _draw_basic_table,
    _header,
    _new_page,
    _paginate_table_rows,
    _plot_panel_title,
    _save_page,
    _table_row_weight,
)
from funmirbench.v2_plots import drawn_tools

LOW_COVERAGE = 0.25
INVENTORY_TOP, INVENTORY_ROWS, ROW_HEIGHT = 0.87, 24.0, 0.0345  # rows per page and their height, as shares of a page
PAGE_TOP, PAGE_BOTTOM, BLOCK_TITLE = 0.975, 0.02, 0.05  # plot page, as shares of the page height
PLOTS = (("coverage", "Coverage", "Genes each tool scores and the responsive and no-response pairs it covers."),
         ("oe_all_genes", "All Genes", "AUROC, Spearman with the shrunken log2FC, and the experiments each tool wins."),
         ("recovery", "Recovery", "Recall by prediction budget and each experiment's best precision at 50% recall."))


def tool_means(summary: pd.DataFrame, tools: pd.DataFrame) -> pd.DataFrame:
    """Mean coverage, AUROC and Spearman over all comparisons; one row per drawn tool, indexed by its label."""
    labels = drawn_tools(tools)
    rows = summary[(summary.subset == "all") & summary.predictor.isin(labels)]
    wide = rows.pivot(index="predictor", columns="metric", values="mean")
    return wide.loc[list(labels), ["coverage", "auroc_d05", "spearman_shrunk"]].rename(index=labels)


def key_results(means: pd.DataFrame) -> list[str]:
    auroc, spearman = means.auroc_d05.idxmax(), means.spearman_shrunk.idxmax()
    covered, sparse = means.coverage.idxmax(), means.coverage.idxmin()
    return [f"Highest mean AUROC: {auroc} ({_format_summary_value(means.auroc_d05[auroc])}).",
            f"Highest mean Spearman: {spearman} ({_format_summary_value(means.spearman_shrunk[spearman])}).",
            f"Highest coverage: {covered} ({_format_summary_value(means.coverage[covered], percent=True)}).",
            f"Lowest coverage: {sparse} ({_format_summary_value(means.coverage[sparse], percent=True)})."]


def coverage_notes(means: pd.DataFrame, n_auroc: int, n_comparisons: int, min_per_state: int) -> list[str]:
    notes = ["Coverage is score availability; a gene a predictor does not score ranks last.",
             f"AUROC on {n_auroc} of {n_comparisons} datasets (at least {min_per_state} responsive and "
             f"{min_per_state} no-response genes).",
             f"Spearman on all {n_comparisons} datasets."]
    low = means[means.coverage < LOW_COVERAGE].coverage
    if len(low):
        notes.append(f"Low coverage (<{LOW_COVERAGE:.0%}): "
                     + "; ".join(f"{t}: {_format_summary_value(v, percent=True)}" for t, v in low.items()) + ".")
    return notes


def draw_inventory(pdf, rows: list[list[str]], label: str) -> None:
    widths = [0.56, 0.44]
    pages = _paginate_table_rows(rows, col_widths=widths, wrap_scale=INVENTORY_WRAP_SCALE, max_weight=INVENTORY_ROWS)
    for number, page_rows in enumerate(pages, start=1):
        fig, ax = _new_page()
        _header(ax, "Dataset Inventory", label + (f" - page {number} of {len(pages)}" if len(pages) > 1 else ""))
        weight = 1 + sum(_table_row_weight(row, col_widths=widths, wrap_scale=INVENTORY_WRAP_SCALE) for row in page_rows)
        height = min(INVENTORY_TOP - 0.02, ROW_HEIGHT * weight)
        _draw_basic_table(ax, page_rows, columns=["Experiment", "Experiment ID"], col_widths=widths,
                          bbox=[0.06, INVENTORY_TOP - height, 0.88, height], font_size=REPORT_TABLE_SIZE,
                          wrap_scale=INVENTORY_WRAP_SCALE)
        _save_page(pdf, fig)


def draw_plots(pdf, plots: pathlib.Path) -> None:
    """All plots on one page at full width, the spare height shared between them."""
    fig, _ = _new_page()
    page_w, page_h = fig.get_size_inches()
    width = 0.92
    images = [plt.imread(plots / f"{stem}.png") for stem, _, _ in PLOTS]
    heights = [image.shape[0] / image.shape[1] * width * page_w / page_h for image in images]
    spare = PAGE_TOP - PAGE_BOTTOM - sum(heights) - BLOCK_TITLE * len(PLOTS)
    gap = max(0.0, spare) / (len(PLOTS) - 1)
    y = PAGE_TOP
    for (_, title, caption), image, height in zip(PLOTS, images, heights):
        _plot_panel_title(fig, title=title, subtitle=caption, title_y=y, subtitle_y=y - 0.023)
        y -= BLOCK_TITLE
        ax = fig.add_axes([0.04, y - height, width, height])
        ax.imshow(image, interpolation="antialiased")
        ax.axis("off")
        y -= height + gap
    _save_page(pdf, fig)


def write_report(path: pathlib.Path, *, comps: pd.DataFrame, set_aside: pd.DataFrame, summary: pd.DataFrame,
                 tools: pd.DataFrame, n_genes: int, min_per_state: int, plots: pathlib.Path) -> None:
    means = tool_means(summary, tools)
    auroc_rows = summary[(summary.subset == "all") & (summary.metric == "auroc_d05")]
    n_auroc = int(auroc_rows.n_comparisons.max())
    rows = [[f"{c.mirna_name} · {c.cell_line}", c.experiment_id] for c in comps.itertuples(index=False)]
    rows += [[f"{c.mirna_name} · {c.cell_line} (set aside: {c.flags.replace(';', ', ').replace('_', ' ')})",
              c.experiment_id] for c in set_aside.itertuples(index=False)]
    label = f"Evaluated experiments ({len(comps)})" + (f", skipped ({len(set_aside)})" if len(set_aside) else "")
    with PdfPages(path) as pdf:
        fig, ax = _new_page()
        _header(ax, "fun-miRBench Benchmark Report", f"{len(comps)} datasets, {len(means)} predictors")
        _draw_basic_table(
            ax,
            [["Datasets evaluated", str(len(comps))],
             ["Datasets skipped", str(len(set_aside))],
             ["Predictors evaluated", str(len(means))],
             ["Responsive genes", "oriented shrunken log2FC > 0; FSOS s-value < 0.05 against d = 0.5"],
             ["No-response genes", "whole shrunken 95% interval of the oriented change below d = 0.5"]],
            columns=["Run setting", "Value"], col_widths=[0.25, 0.75], bbox=[0.04, 0.605, 0.92, 0.25],
            font_size=REPORT_TABLE_SIZE)
        next_y = _bullet_block(
            ax, "Evaluation design",
            ["Higher aligned scores mean stronger predicted targeting.",
             f"Every predictor is scored on the same {n_genes:,} genes; a gene it does not score ranks last.",
             "Unresolved genes are left out of the AUROC.",
             "Overexpression (OE) datasets only."],
            x=0.04, y=0.54, width=0.92)
        next_y = _bullet_block(
            ax, "How to read this report",
            ["The cross-dataset table is the primary summary.",
             "Low-coverage predictors stay in the table and the plots."],
            x=0.04, y=next_y - 0.020, width=0.92)
        _draw_basic_table(
            ax,
            [["per_comparison.tsv", "Every metric per dataset and predictor"],
             ["summary.tsv", "Means over datasets with 95% intervals"],
             ["paired_vs_reference.tsv", "Each predictor against seed + expression"],
             ["decoy.tsv", "Each predictor's scores for unrelated miRNAs"],
             ["coverage.tsv, recovery.tsv", "Tables behind the plots"],
             ["precision_at_recall.tsv", "Table behind the recovery plot"]],
            columns=["Primary file", "Purpose"], col_widths=[0.50, 0.50],
            bbox=[0.04, 0.03, 0.92, next_y - 0.02 - 0.03], font_size=10.2, wrap_scale=135)
        _save_page(pdf, fig)

        fig, ax = _new_page()
        _header(ax, "Cross-Dataset Predictor Summary", "Mean performance across datasets, every predictor on all genes.")
        next_y = _bullet_block(ax, "Key results", key_results(means), x=0.04, y=0.84, width=0.92)
        next_y = _bullet_block(ax, "Coverage notes", coverage_notes(means, n_auroc, len(comps), min_per_state),
                               x=0.04, y=next_y - 0.018, width=0.92)
        _draw_basic_table(
            ax,
            [[tool, _format_summary_value(r.coverage, percent=True), _format_summary_value(r.auroc_d05),
              _format_summary_value(r.spearman_shrunk)] for tool, r in means.iterrows()],
            columns=["Predictor", "Coverage", "AUROC", "Spearman"], col_widths=[0.40, 0.20, 0.20, 0.20],
            bbox=[0.035, 0.03, 0.93, next_y - 0.02 - 0.03], font_size=REPORT_TABLE_SIZE)
        _save_page(pdf, fig)

        draw_inventory(pdf, rows, label)

        with plt.rc_context({"savefig.dpi": 300}):  # the plots are PNGs; the default would blur them in the PDF
            draw_plots(pdf, plots)
