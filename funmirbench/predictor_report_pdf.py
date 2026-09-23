"""Publication-quality per-predictor PDF reports."""

from __future__ import annotations

import pathlib

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

from funmirbench.evaluate import REPORT_PAGE_SIZE, _is_publication_tool
from funmirbench.experiment_labels import (
    format_experiment_id_label,
    format_experiment_label,
)


BLUE = "#17324D"
MUTED = "#5B6577"
RULE = "#D8DEE9"
BOX_FACE = "#F5F8FC"
BOX_EDGE = "#D8E2EF"
TABLE_HEADER = "#E9F1FB"
TABLE_ALT = "#F9FBFD"
REPORT_TITLE_SIZE = 23
REPORT_SUBTITLE_SIZE = 12.2
REPORT_SECTION_SIZE = 13.4
REPORT_TABLE_SIZE = 10.2
REPORT_PLOT_TITLE_SIZE = 13.4
REPORT_PLOT_SUBTITLE_SIZE = 10.4


def _metric_value(value, *, percent=False):
    if value is None or pd.isna(value):
        return "NA"
    value = float(value)
    return f"{value:.1%}" if percent else f"{value:.3f}"


def _short_path(path):
    if path is None:
        return "NA"
    path = pathlib.Path(str(path))
    return path.name or str(path)


def _gt_rule_text(fdr_threshold, effect_threshold):
    effect_text = f"effect > {float(effect_threshold)}"
    if fdr_threshold is None:
        return effect_text
    return f"FDR < {float(fdr_threshold)}; {effect_text}"


def _new_page():
    fig = plt.figure(figsize=REPORT_PAGE_SIZE)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis("off")
    fig.patch.set_facecolor("white")
    return fig, ax


def _save_page(pdf, fig):
    fig.patch.set_facecolor("white")
    pdf.savefig(fig, facecolor="white")
    plt.close(fig)


def _header(ax, experiment_label, predictor_label, dataset_id):
    ax.text(
        0.06,
        0.965,
        experiment_label,
        fontsize=18.0,
        fontweight="bold",
        color=BLUE,
        va="top",
        ha="left",
        family="DejaVu Sans",
    )
    ax.text(
        0.06,
        0.922,
        predictor_label,
        fontsize=REPORT_TITLE_SIZE,
        fontweight="bold",
        color=BLUE,
        va="top",
        ha="left",
        family="DejaVu Sans",
    )
    ax.text(
        0.06,
        0.882,
        format_experiment_id_label(dataset_id),
        fontsize=10.2,
        color=MUTED,
        va="top",
        ha="left",
        family="DejaVu Sans",
    )
    ax.add_line(plt.Line2D([0.06, 0.94], [0.858, 0.858], color=RULE, linewidth=1.2))


def _metric_card(ax, label, value, *, x, y):
    ax.text(
        x,
        y,
        f"{label}\n{value}",
        fontsize=11.4,
        fontweight="bold",
        color=BLUE,
        va="top",
        ha="left",
        bbox={"boxstyle": "round,pad=0.42", "facecolor": BOX_FACE, "edgecolor": BOX_EDGE},
    )


def _key_value_table(ax, rows, *, bbox):
    table = ax.table(
        cellText=rows,
        colLabels=["Metric", "Value"],
        colWidths=[0.58, 0.42],
        cellLoc="left",
        colLoc="left",
        bbox=bbox,
    )
    table.auto_set_font_size(False)
    table.set_fontsize(REPORT_TABLE_SIZE)
    table.scale(1.0, 1.28)
    for (row, col), cell in table.get_celld().items():
        cell.PAD = 0.075
        cell.set_edgecolor("#D9E2EC")
        cell.set_linewidth(0.7)
        if row == 0:
            cell.set_facecolor(TABLE_HEADER)
            cell.set_text_props(weight="bold", color=BLUE)
        else:
            cell.set_facecolor("#FFFFFF" if row % 2 else TABLE_ALT)
            if col == 1:
                cell.set_text_props(ha="right")
    return table


def _wide_table(ax, rows, *, title, columns, col_widths, bbox, font_size=REPORT_TABLE_SIZE, right_align_columns=None):
    right_align_columns = set(right_align_columns or [])
    ax.text(
        bbox[0],
        bbox[1] + bbox[3] + 0.018,
        title,
        fontsize=REPORT_SECTION_SIZE,
        fontweight="bold",
        color="#2F5D8C",
        va="bottom",
        ha="left",
    )
    table = ax.table(
        cellText=rows,
        colLabels=columns,
        colWidths=col_widths,
        cellLoc="left",
        colLoc="left",
        bbox=bbox,
    )
    table.auto_set_font_size(False)
    table.set_fontsize(font_size)
    for (row, col), cell in table.get_celld().items():
        cell.PAD = 0.075
        cell.set_edgecolor("#D9E2EC")
        cell.set_linewidth(0.7)
        if row == 0:
            cell.set_facecolor(TABLE_HEADER)
            cell.set_text_props(weight="bold", color=BLUE)
        else:
            cell.set_facecolor("#FFFFFF" if row % 2 else TABLE_ALT)
            if col in right_align_columns:
                cell.set_text_props(ha="right")
    return table


def _plot_panel_title(fig, *, title, subtitle, title_y, subtitle_y, x=0.04):
    fig.text(
        x,
        title_y,
        title,
        fontsize=REPORT_PLOT_TITLE_SIZE,
        fontweight="bold",
        color=BLUE,
        va="top",
        ha="left",
        family="DejaVu Sans",
    )
    fig.text(
        x,
        subtitle_y,
        subtitle,
        fontsize=REPORT_PLOT_SUBTITLE_SIZE,
        color=MUTED,
        va="top",
        ha="left",
        family="DejaVu Sans",
    )


def _plot_panel_subtitle(label):
    subtitles = {
        "Score vs expected effect": "Scores compared with perturbation-aware effect.",
        "GSEA enrichment": "Running enrichment of GT positives along ranked genes.",
        "Precision-recall": "Precision and recall over score thresholds.",
        "ROC": "True-positive rate versus false-positive rate.",
        "Top 10% GT-positive heatmap": "Top GT positives ordered by perturbation-aware effect.",
        "Precision-recall by predictor": "Predictors evaluated on their own scored genes.",
        "ROC by predictor": "Predictors evaluated on their own scored genes.",
        "Top-prediction effect CDFs": "Effect distributions among top-ranked predictions.",
    }
    return subtitles.get(str(label), "Dataset-level diagnostic plot.")


def _plot_grid_page(
    pdf,
    *,
    experiment_label,
    predictor_label,
    dataset_id,
    plot_paths,
):
    fig, ax = _new_page()
    _header(ax, experiment_label, predictor_label, dataset_id)
    panels = [
        {"title_y": 0.825, "subtitle_y": 0.800, "image_box": [0.04, 0.445, 0.92, 0.325]},
        {"title_y": 0.405, "subtitle_y": 0.380, "image_box": [0.04, 0.025, 0.92, 0.325]},
    ]
    for (label, path), panel in zip(plot_paths, panels):
        path = pathlib.Path(path)
        if not path.is_file():
            continue
        _plot_panel_title(
            fig,
            title=str(label),
            subtitle=_plot_panel_subtitle(label),
            title_y=panel["title_y"],
            subtitle_y=panel["subtitle_y"],
        )
        image = plt.imread(path)
        box = panel["image_box"]
        image_ax = fig.add_axes(box)
        image_ax.imshow(image, interpolation="antialiased")
        image_ax.axis("off")
    _save_page(pdf, fig)


def _chunks(items, size):
    for start in range(0, len(items), size):
        yield items[start : start + size]


def _dataset_context_plots(plots_dir):
    paths = []
    heatmap_path = plots_dir / "heatmaps" / "top_10pct_positive_genes.png"
    if heatmap_path.is_file():
        paths.append(("Top 10% GT-positive heatmap", heatmap_path))

    comparison_dir = plots_dir / "comparisons"
    preferred_names = [
        "precision_recall_all_scored.png",
        "roc_all_scored.png",
        "top_100_effect_cdfs.png",
    ]
    seen = {path for _, path in paths}
    for filename in preferred_names:
        path = comparison_dir / filename
        if path.is_file() and path not in seen:
            paths.append((_comparison_plot_label(path), path))
            seen.add(path)
    if comparison_dir.is_dir():
        for path in sorted(comparison_dir.glob("*.png")):
            if "_common" in path.stem or path in seen:
                continue
            paths.append((_comparison_plot_label(path), path))
            seen.add(path)
    return paths


def _comparison_plot_label(path):
    labels = {
        "precision_recall_all_scored": "Precision-recall by predictor",
        "roc_all_scored": "ROC by predictor",
        "top_100_effect_cdfs": "Top-prediction effect CDFs",
    }
    return labels.get(path.stem, path.stem.replace("_", " ").title())


def write_publication_predictor_reports(
    *,
    reports_dir,
    plots_dir,
    dataset_id,
    mirna,
    cell_line,
    perturbation,
    geo_accession,
    de_table_path,
    predictor_output_paths,
    metric_rows,
    tool_labels,
    fdr_threshold,
    effect_threshold,
    common_prediction_summary=None,
    skipped_tool_rows=None,
):
    del common_prediction_summary
    reports_dir = pathlib.Path(reports_dir)
    plots_dir = pathlib.Path(plots_dir)
    written = []
    report_rows = [row for row in metric_rows if _is_publication_tool(row.get("tool_id"))]
    evaluated_tool_ids = {str(row.get("tool_id")) for row in report_rows}
    skipped_rows = [
        row
        for row in (skipped_tool_rows or [])
        if _is_publication_tool(row.get("tool_id"))
        and str(row.get("tool_id")) not in evaluated_tool_ids
    ]
    if not report_rows:
        report_rows = list(metric_rows)
    rows_to_write = [("evaluated", row) for row in report_rows] + [("skipped", row) for row in skipped_rows]
    gt_rule = _gt_rule_text(fdr_threshold, effect_threshold)
    for row_status, row in rows_to_write:
        tool_id = str(row.get("tool_id"))
        label = str(tool_labels.get(tool_id, tool_id))
        report_path = reports_dir / f"{dataset_id}__{tool_id}_evaluation_report.pdf"
        predictor_dir = plots_dir / "predictors" / tool_id
        plot_paths = [
            ("Score vs expected effect", predictor_dir / "score_vs_expected_effect.png"),
            ("GSEA enrichment", predictor_dir / "gsea_enrichment.png"),
            ("Precision-recall", predictor_dir / "precision_recall_curve.png"),
            ("ROC", predictor_dir / "roc_curve.png"),
        ]

        with PdfPages(report_path) as pdf:
            experiment_label = format_experiment_label(
                geo_accession=geo_accession,
                mirna=mirna,
                cell_line=cell_line,
                perturbation=perturbation,
            )
            fig, ax = _new_page()
            _header(ax, experiment_label, label, dataset_id)
            if row_status == "skipped":
                rows_total = int(row.get("rows_total", 0) or 0)
                rows_scored = int(row.get("rows_scored", 0) or 0)
                coverage = rows_scored / rows_total if rows_total else float("nan")
                skip_reason = str(row.get("skip_reason") or "Predictor could not be evaluated.")
                cards = [
                    ("Status", "Skipped", 0.06),
                    ("Coverage", _metric_value(coverage, percent=True), 0.29),
                    ("Rows scored", f"{rows_scored:,}", 0.52),
                    ("GT pos. scored", f"{int(row.get('positives_scored', 0) or 0):,}", 0.75),
                ]
                for label_text, value, x in cards:
                    _metric_card(ax, label_text, value, x=x, y=0.81)
                _key_value_table(
                    ax,
                    [
                        ["Rows total", f"{rows_total:,}"],
                        ["Rows scored", f"{rows_scored:,}"],
                        ["Rows missing score", f"{int(row.get('rows_missing_score', 0) or 0):,}"],
                        ["GT positives total", f"{int(row.get('positives_total', 0) or 0):,}"],
                        ["GT positives scored", f"{int(row.get('positives_scored', 0) or 0):,}"],
                    ],
                    bbox=[0.06, 0.50, 0.40, 0.18],
                )
                _wide_table(
                    ax,
                    [
                        ["Skip reason", skip_reason],
                        ["Diagnostic plots", "Not generated for this dataset"],
                    ],
                    title="Evaluation status",
                    columns=["Item", "Value"],
                    col_widths=[0.32, 0.68],
                    bbox=[0.52, 0.52, 0.42, 0.13],
                    font_size=REPORT_TABLE_SIZE,
                )
                _wide_table(
                    ax,
                    [
                        ["GT positives", gt_rule],
                        ["Effect sign", "-logFC for Overexpression; +logFC for Knockout/Knockdown"],
                        ["Metric scope", "Usable ground truth plus positive and background scored genes"],
                    ],
                    title="Evaluation rule",
                    columns=["Rule", "Value"],
                    col_widths=[0.22, 0.78],
                    bbox=[0.04, 0.30, 0.92, 0.12],
                    font_size=REPORT_TABLE_SIZE,
                )
                _wide_table(
                    ax,
                    [
                        ["Experiment ID", dataset_id],
                        ["GEO accession", geo_accession or "NA"],
                        ["DE table file", _short_path(de_table_path)],
                        ["Predictor file", _short_path(predictor_output_paths.get(tool_id))],
                    ],
                    title="Provenance",
                    columns=["Field", "Value"],
                    col_widths=[0.22, 0.78],
                    bbox=[0.04, 0.10, 0.92, 0.11],
                    font_size=REPORT_TABLE_SIZE,
                )
                _save_page(pdf, fig)
                written.append(report_path)
                continue

            cards = [
                ("Coverage", _metric_value(row.get("coverage"), percent=True), 0.06),
                ("Positive cov.", _metric_value(row.get("positive_coverage"), percent=True), 0.29),
                ("APS", _metric_value(row.get("aps")), 0.52),
                ("AUROC", _metric_value(row.get("auroc")), 0.75),
            ]
            for label_text, value, x in cards:
                _metric_card(ax, label_text, value, x=x, y=0.81)

            metric_table_rows = [
                ["Coverage", _metric_value(row.get("coverage"), percent=True)],
                ["Positive coverage", _metric_value(row.get("positive_coverage"), percent=True)],
                ["APS", _metric_value(row.get("aps"))],
                ["PR-AUC", _metric_value(row.get("pr_auc"))],
                ["AUROC", _metric_value(row.get("auroc"))],
                ["Spearman", _metric_value(row.get("spearman"))],
                ["Pearson", _metric_value(row.get("pearson"))],
                ["GSEA ES", _metric_value(row.get("gsea_es"))],
            ]
            coverage_rows = [
                ["Rows total", f"{int(row.get('rows_total', 0)):,}"],
                ["Rows scored", f"{int(row.get('rows_scored', 0)):,}"],
                ["Rows missing score", f"{int(row.get('rows_missing_score', 0)):,}"],
                ["GT positives total", f"{int(row.get('positives_total', 0)):,}"],
                ["GT positives scored", f"{int(row.get('positives_scored', 0)):,}"],
            ]
            _key_value_table(ax, metric_table_rows, bbox=[0.04, 0.54, 0.42, 0.20])
            _key_value_table(ax, coverage_rows, bbox=[0.52, 0.56, 0.42, 0.17])
            _wide_table(
                ax,
                [
                    ["GT positives", gt_rule],
                    ["Effect sign", "-logFC for Overexpression; +logFC for Knockout/Knockdown"],
                    ["Score direction", "Higher aligned scores mean stronger predicted targeting"],
                    ["Metric scope", "Usable ground truth and an available score"],
                ],
                title="Evaluation rule",
                columns=["Rule", "Value"],
                col_widths=[0.22, 0.78],
                bbox=[0.04, 0.365, 0.92, 0.115],
                font_size=REPORT_TABLE_SIZE,
            )
            _wide_table(
                ax,
                [
                    ["Experiment ID", dataset_id],
                    ["GEO accession", geo_accession or "NA"],
                    ["DE table file", _short_path(de_table_path)],
                    ["Predictor file", _short_path(predictor_output_paths.get(tool_id))],
                ],
                title="Provenance",
                columns=["Field", "Value"],
                col_widths=[0.22, 0.78],
                bbox=[0.04, 0.175, 0.92, 0.105],
                font_size=REPORT_TABLE_SIZE,
            )
            _save_page(pdf, fig)
            for page_index, plot_chunk in enumerate(_chunks(plot_paths, 2), start=1):
                suffix = "" if len(plot_paths) <= 2 else f" ({page_index})"
                _plot_grid_page(
                    pdf,
                    experiment_label=experiment_label,
                    predictor_label=f"{label}{suffix}",
                    dataset_id=dataset_id,
                    plot_paths=plot_chunk,
                )
            dataset_context = _dataset_context_plots(plots_dir)
            for page_index, context_chunk in enumerate(_chunks(dataset_context, 2), start=1):
                suffix = "" if len(dataset_context) <= 2 else f" ({page_index})"
                _plot_grid_page(
                    pdf,
                    experiment_label=experiment_label,
                    predictor_label=f"All selected predictors{suffix}",
                    dataset_id=dataset_id,
                    plot_paths=context_chunk,
                )
        written.append(report_path)
    return written
