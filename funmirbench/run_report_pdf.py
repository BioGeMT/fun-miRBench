"""Publication-focused PDF report helpers."""

from __future__ import annotations

import pathlib
import textwrap

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

from funmirbench.evaluate import REPORT_PAGE_SIZE, _publication_tool_ids, _tool_label
from funmirbench.benchmark_reports import (
    MIN_HEADLINE_COVERAGE,
    _format_summary_value,
    _load_cross_dataset_summary,
)


PUBLICATION_BLUE = "#17324D"
PUBLICATION_MUTED = "#5B6577"
PUBLICATION_RULE = "#D8DEE9"
PUBLICATION_TABLE_HEADER = "#E9F1FB"
PUBLICATION_TABLE_ALT = "#F9FBFD"
REPORT_TITLE_SIZE = 23
REPORT_SUBTITLE_SIZE = 12.2
REPORT_SECTION_SIZE = 13.4
REPORT_BODY_SIZE = 11.2
REPORT_TABLE_SIZE = 10.2
REPORT_PLOT_TITLE_SIZE = 13.4
REPORT_PLOT_SUBTITLE_SIZE = 10.4


def _new_page():
    fig = plt.figure(figsize=REPORT_PAGE_SIZE)
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    ax.axis("off")
    fig.patch.set_facecolor("white")
    return fig, ax


def _save_page(pdf, fig):
    fig.patch.set_facecolor("white")
    pdf.savefig(fig, facecolor="white")
    plt.close(fig)


def _header(ax, title, subtitle=None):
    ax.text(0.06, 0.955, title, fontsize=REPORT_TITLE_SIZE, fontweight="bold", color=PUBLICATION_BLUE, va="top", ha="left", family="DejaVu Sans")
    if subtitle:
        ax.text(0.06, 0.916, subtitle, fontsize=REPORT_SUBTITLE_SIZE, color=PUBLICATION_MUTED, va="top", ha="left", family="DejaVu Sans")
    ax.add_line(plt.Line2D([0.06, 0.94], [0.895, 0.895], color=PUBLICATION_RULE, linewidth=1.2))


def _text_block(ax, title, lines, *, x, y, width, body_size=9.6, title_size=11.8):
    ax.text(x, y, title, fontsize=title_size, fontweight="bold", color="#2F5D8C", va="top", ha="left", family="DejaVu Sans")
    current_y = y - 0.032
    wrap_width = max(30, int(width * 108))
    for line in lines:
        for chunk in textwrap.wrap(str(line), width=wrap_width) or [""]:
            ax.text(x, current_y, chunk, fontsize=body_size, color="#22303C", va="top", ha="left", family="DejaVu Sans")
            current_y -= 0.025
        current_y -= 0.006
    return current_y


def _bullet_block(ax, title, lines, *, x, y, width, body_size=REPORT_BODY_SIZE, title_size=REPORT_SECTION_SIZE):
    ax.text(x, y, title, fontsize=title_size, fontweight="bold", color="#2F5D8C", va="top", ha="left", family="DejaVu Sans")
    current_y = y - 0.038
    wrap_width = max(30, int(width * 130))
    for line in lines:
        chunks = textwrap.wrap(str(line), width=wrap_width) or [""]
        for index, chunk in enumerate(chunks):
            prefix = "- " if index == 0 else "  "
            ax.text(
                x,
                current_y,
                prefix + chunk,
                fontsize=body_size,
                color="#22303C",
                va="top",
                ha="left",
                family="DejaVu Sans",
            )
            current_y -= 0.030
        current_y -= 0.007
    return current_y


def _summary_box(ax, label, value, *, x, y):
    ax.text(
        x,
        y,
        f"{label}\n{value}",
        fontsize=10.1,
        fontweight="bold",
        color=PUBLICATION_BLUE,
        va="top",
        ha="left",
        family="DejaVu Sans",
        bbox={"boxstyle": "round,pad=0.42", "facecolor": "#F5F8FC", "edgecolor": "#D8E2EF"},
    )


def _plot_panel_title(fig, *, title, subtitle, title_y, subtitle_y, x=0.04):
    fig.text(
        x,
        title_y,
        title,
        fontsize=REPORT_PLOT_TITLE_SIZE,
        fontweight="bold",
        color=PUBLICATION_BLUE,
        va="top",
        ha="left",
        family="DejaVu Sans",
    )
    fig.text(
        x,
        subtitle_y,
        subtitle,
        fontsize=REPORT_PLOT_SUBTITLE_SIZE,
        color=PUBLICATION_MUTED,
        va="top",
        ha="left",
        family="DejaVu Sans",
    )


def _draw_basic_table(ax, rows, *, columns, col_widths, bbox, font_size=REPORT_TABLE_SIZE, wrap_scale=96):
    def _wrap_cell(value, width_hint):
        text = str(value)
        max_chars = max(16, int(width_hint * wrap_scale))
        return "\n".join(
            textwrap.wrap(
                text,
                width=max_chars,
                break_long_words=True,
                break_on_hyphens=False,
            )
        )

    wrapped_rows = [[_wrap_cell(value, width) for value, width in zip(row, col_widths)] for row in rows]
    wrapped_columns = [_wrap_cell(value, width) for value, width in zip(columns, col_widths)]
    table = ax.table(
        cellText=wrapped_rows,
        colLabels=wrapped_columns,
        colWidths=col_widths,
        cellLoc="left",
        colLoc="left",
        bbox=bbox,
    )
    table.auto_set_font_size(False)
    table.set_fontsize(font_size)
    table.scale(1.0, 1.35)
    row_weights = [1.0]
    for row in wrapped_rows:
        row_weights.append(max(1.0, max(str(value).count("\n") + 1 for value in row) * 0.82))
    weight_total = sum(row_weights)
    for (row, col), cell in table.get_celld().items():
        cell.PAD = 0.095
        cell.set_height(bbox[3] * row_weights[row] / weight_total)
        cell.set_edgecolor("#D9E2EC")
        cell.set_linewidth(0.7)
        if row == 0:
            cell.set_facecolor(PUBLICATION_TABLE_HEADER)
            cell.set_text_props(weight="bold", color=PUBLICATION_BLUE)
        else:
            cell.set_facecolor("#FFFFFF" if row % 2 else PUBLICATION_TABLE_ALT)
    return table


def _format_summary_table(summary_df):
    display_df = summary_df[["tool_id", "coverage_mean", "positive_coverage_mean", "aps_mean", "pr_auc_mean", "spearman_mean", "auroc_mean"]].copy()
    display_df.columns = ["Predictor", "Coverage", "Positive cov.", "APS", "PR-AUC", "Spearman", "AUROC"]
    percent_columns = {"Coverage", "Positive cov."}
    for column in display_df.columns[1:]:
        percent = column in percent_columns
        display_df[column] = display_df[column].map(lambda value: _format_summary_value(value, percent=percent))
    return display_df


def _draw_summary_table(ax, summary_df, *, bbox):
    display_df = _format_summary_table(summary_df)
    table = ax.table(
        cellText=display_df.values.tolist(),
        colLabels=display_df.columns.tolist(),
        colWidths=[0.19, 0.12, 0.17, 0.11, 0.12, 0.14, 0.12],
        cellLoc="center",
        colLoc="center",
        bbox=bbox,
    )
    table.auto_set_font_size(False)
    table.set_fontsize(REPORT_TABLE_SIZE)
    table.scale(1.0, 1.75)
    for (row, col), cell in table.get_celld().items():
        cell.PAD = 0.06
        cell.set_edgecolor("#D9E2EC")
        cell.set_linewidth(0.7)
        if row == 0:
            cell.set_facecolor(PUBLICATION_TABLE_HEADER)
            cell.set_text_props(weight="bold", color=PUBLICATION_BLUE, fontsize=REPORT_TABLE_SIZE - 0.6)
        else:
            cell.set_facecolor("#FFFFFF" if row % 2 else PUBLICATION_TABLE_ALT)
            if col == 0:
                cell.set_text_props(ha="left")
    return table


def _load_common_prediction_summary(combined_outputs):
    path = combined_outputs.get("tables", {}).get("common_prediction_summary")
    if not path:
        return None
    path = pathlib.Path(path)
    if not path.is_file():
        return None
    try:
        return pd.read_csv(path, sep="\t")
    except pd.errors.EmptyDataError:
        return None


def _draw_common_prediction_page(pdf, combined_outputs):
    summary = _load_common_prediction_summary(combined_outputs)
    if summary is None or summary.empty:
        return
    selected = summary[summary["summary_type"].isin(["report_common_set", "all_real_predictors_common_set"])].copy()
    if selected.empty:
        return
    def _format_tools(value):
        return " + ".join(_tool_label(tool_id.strip()) for tool_id in str(value).split(",") if tool_id.strip())

    selected["tools"] = selected["tools"].map(_format_tools)
    selected["Common predictions"] = selected.apply(
        lambda row: f"{int(row['rows_common']):,}/{int(row['rows_total']):,} ({float(row['percent_common']):.1%})",
        axis=1,
    )
    selected["Set"] = selected["summary_type"].map(
        {
            "report_common_set": "Report common set",
            "all_real_predictors_common_set": "All predictors",
        }
    )
    display = selected[["dataset_id", "Set", "tools", "Common predictions"]].copy()
    display.columns = ["Dataset", "Set", "Predictors", "Common predictions"]
    fig, ax = _new_page()
    _header(
        ax,
        "Common Prediction Coverage",
        "Genes with non-missing scores for each predictor set; correlation heatmaps are omitted.",
    )
    _draw_basic_table(
        ax,
        display.values.tolist(),
        columns=display.columns.tolist(),
        col_widths=[0.30, 0.18, 0.34, 0.18],
        bbox=[0.04, 0.10, 0.92, 0.74],
        font_size=REPORT_TABLE_SIZE,
    )
    _save_page(pdf, fig)


def _coverage_note_bullets(summary_df):
    if summary_df is None or summary_df.empty:
        return ["Cross-dataset predictor summary is unavailable for this run."]
    bullets = [
        "Coverage is score availability.",
        f"Headline rankings require >= {MIN_HEADLINE_COVERAGE:.0%} mean coverage.",
        "Positive coverage is score availability among GT positives.",
    ]
    sparse = summary_df[summary_df["coverage_mean"].astype(float) < 0.25].copy()
    if not sparse.empty:
        sparse_labels = [
            f"{row.tool_id}: {_format_summary_value(row.coverage_mean, percent=True)}"
            for row in sparse.itertuples(index=False)
        ]
        bullets.append("Low coverage (<25%): " + "; ".join(sparse_labels) + ".")
    coverage_gap = (
        summary_df["positive_coverage_mean"].astype(float)
        - summary_df["coverage_mean"].astype(float)
    )
    if coverage_gap.notna().any():
        row = summary_df.loc[coverage_gap.idxmax()]
        bullets.append(
            "Largest GT-positive coverage enrichment: "
            f"{row['tool_id']} {_format_summary_value(row['positive_coverage_mean'], percent=True)} vs "
            f"{_format_summary_value(row['coverage_mean'], percent=True)} overall."
        )
    return bullets


def _key_result_bullets(summary_df):
    if summary_df is None or summary_df.empty:
        return []
    work = summary_df.copy()
    for column in ["coverage_mean", "positive_coverage_mean", "aps_mean", "auroc_mean"]:
        work[column] = pd.to_numeric(work[column], errors="coerce")
    eligible = work[work["coverage_mean"] >= MIN_HEADLINE_COVERAGE].copy()
    bullets = []
    if eligible.empty:
        bullets.append(f"No predictor reached >= {MIN_HEADLINE_COVERAGE:.0%} mean coverage for headline ranking.")
    else:
        aps_row = eligible.sort_values(["aps_mean", "auroc_mean"], ascending=False).iloc[0]
        bullets.append(
            f"Coverage-qualified APS leader: {aps_row['tool_id']} "
            f"(APS {_format_summary_value(aps_row['aps_mean'])})."
        )
        auroc_row = eligible.sort_values(["auroc_mean", "aps_mean"], ascending=False).iloc[0]
        bullets.append(
            f"Coverage-qualified AUROC leader: {auroc_row['tool_id']} "
            f"(AUROC {_format_summary_value(auroc_row['auroc_mean'])})."
        )
    if work["positive_coverage_mean"].notna().any():
        positive_row = work.sort_values("positive_coverage_mean", ascending=False).iloc[0]
        bullets.append(
            f"Highest GT-positive coverage: {positive_row['tool_id']} "
            f"({_format_summary_value(positive_row['positive_coverage_mean'], percent=True)})."
        )
    if work["coverage_mean"].notna().any():
        sparse_row = work.sort_values("coverage_mean", ascending=True).iloc[0]
        bullets.append(
            f"Lowest overall coverage: {sparse_row['tool_id']} "
            f"({_format_summary_value(sparse_row['coverage_mean'], percent=True)})."
        )
    return bullets


def _plot_items(combined_outputs):
    metric_labels = {
        "coverage": "Coverage",
        "positive_coverage": "Positive Coverage",
        "aps": "APS",
        "pr_auc": "PR-AUC",
        "spearman": "Spearman",
        "auroc": "AUROC",
    }
    descriptions = {
        "predictor_combination_expanded_frontier": (
            "Predictor-Combination Frontier",
            "Mean positive coverage versus mean APS.",
        ),
    }
    for metric_name, metric_label in metric_labels.items():
        descriptions[f"cross_dataset_{metric_name}_distribution"] = (
            f"Cross-Dataset {metric_label} Distribution",
            "Per-dataset values summarized for each predictor.",
        )
    descriptions.update(
        {
            "positive_background_local_rank_distributions": (
                "Local Rank Distributions",
                "GT positives and background genes by dataset-local rank.",
            ),
            "positive_background_local_rank_counts": (
                "Local Rank Counts",
                "Binned gene counts by dataset-local rank.",
            ),
            "positive_background_global_rank_distributions": (
                "Global Rank Distributions",
                "GT positives and background genes by predictor-file rank.",
            ),
            "positive_background_global_rank_counts": (
                "Global Rank Counts",
                "Binned gene counts by predictor-file rank.",
            ),
            "positive_recovery_fraction_by_prediction_count": (
                "GT-Positive Recovery by Prediction Count",
                "Mean fraction of positives recovered as top predictions increase.",
            ),
        }
    )
    for key, (title, caption) in descriptions.items():
        path = combined_outputs.get("plots", {}).get(key)
        if path:
            path = pathlib.Path(path)
            if path.is_file():
                yield key, title, caption, path


def _draw_plot_pair_page(pdf, items):
    fig, _ = _new_page()
    panels = [
        {"title_y": 0.965, "subtitle_y": 0.940, "image_box": [0.04, 0.505, 0.92, 0.410]},
        {"title_y": 0.485, "subtitle_y": 0.460, "image_box": [0.04, 0.025, 0.92, 0.410]},
    ]
    for (_, title, caption, path), panel in zip(items, panels):
        _plot_panel_title(
            fig,
            title=title,
            subtitle=caption,
            title_y=panel["title_y"],
            subtitle_y=panel["subtitle_y"],
        )
        box = panel["image_box"]
        image_ax = fig.add_axes(box)
        image_ax.imshow(plt.imread(path), interpolation="antialiased")
        image_ax.axis("off")
    _save_page(pdf, fig)


def _draw_rank_pair_page(pdf, combined_outputs, *, rank_type):
    if rank_type == "local":
        density_key = "positive_background_local_rank_distributions"
        count_key = "positive_background_local_rank_counts"
        density_title = "Local Rank Distributions"
        density_caption = "GT positives and background genes by dataset-local rank."
        count_title = "Local Rank Counts"
        count_caption = "Binned gene counts by dataset-local rank."
    else:
        density_key = "positive_background_global_rank_distributions"
        count_key = "positive_background_global_rank_counts"
        density_title = "Global Rank Distributions"
        density_caption = "GT positives and background genes by predictor-file rank."
        count_title = "Global Rank Counts"
        count_caption = "Binned gene counts by predictor-file rank."

    plots = combined_outputs.get("plots", {})
    density_path = pathlib.Path(plots.get(density_key, ""))
    count_path = pathlib.Path(plots.get(count_key, ""))
    if not density_path.is_file() or not count_path.is_file():
        return set()

    _draw_plot_pair_page(
        pdf,
        [
            (density_key, density_title, density_caption, density_path),
            (count_key, count_title, count_caption, count_path),
        ],
    )
    return {density_key, count_key}


def _run_subtitle(dataset_count, predictor_count, config_path):
    return (
        f"{dataset_count} datasets, {predictor_count} predictors, "
        f"generated from {pathlib.Path(config_path).name}"
    )


def _dataset_report_rows(dataset_outputs):
    return [
        {
            "Dataset": item["dataset_id"],
            "miRNA": item["mirna"],
            "Perturbation": item["perturbation"],
            "Cell line": item["cell_line"],
        }
        for item in dataset_outputs
    ]


def write_publication_run_pdf_report(
    out_dir,
    *,
    config_path,
    dataset_outputs,
    tool_ids,
    metric_tables,
    combined_outputs,
    fdr_threshold,
    effect_threshold,
    predictor_top_fraction,
    skipped_datasets=None,
):
    del metric_tables
    skipped_datasets = list(skipped_datasets or [])
    summary_df = _load_cross_dataset_summary(combined_outputs)
    report_path = pathlib.Path(out_dir) / "REPORT.pdf"
    display_tool_ids = _publication_tool_ids(tool_ids) or list(tool_ids)

    with PdfPages(report_path) as pdf:
        fig, ax = _new_page()
        _header(
            ax,
            "fun-miRBenvh Benchmark Report",
            _run_subtitle(len(dataset_outputs), len(display_tool_ids), config_path),
        )
        _draw_basic_table(
            ax,
            [
                ["Datasets evaluated", str(len(dataset_outputs))],
                ["Datasets skipped", str(len(skipped_datasets))],
                ["Predictors evaluated", str(len(display_tool_ids))],
                [
                    "GT positives",
                    (
                        f"effect > {float(effect_threshold)}; sign-aware (-logFC Overexpression, +logFC Knockout/Knockdown)"
                        if fdr_threshold is None
                        else f"FDR < {float(fdr_threshold)}; effect > {float(effect_threshold)}; sign-aware (-logFC Overexpression, +logFC Knockout/Knockdown)"
                    ),
                ],
                ["Top fraction", f"{predictor_top_fraction:.0%} exact top-k per predictor"],
            ],
            columns=["Run setting", "Value"],
            col_widths=[0.25, 0.75],
            bbox=[0.04, 0.605, 0.92, 0.25],
            font_size=REPORT_TABLE_SIZE,
        )
        next_y = _bullet_block(
            ax,
            "Evaluation design",
            [
                "Higher aligned scores mean stronger predicted targeting.",
                "Dataset plots use dataset-local, tie-aware dense ranks.",
                "Cross-dataset global ranks use each full standardized predictor file.",
            ],
            x=0.04,
            y=0.54,
            width=0.92,
        )
        _bullet_block(
            ax,
            "How to read this report",
            [
                "The cross-dataset table is the primary summary.",
                "Sparse predictors remain in tables, not headline rankings.",
                "Common-prediction percentages report predictor overlap.",
                "Coverage diagnostics are retained in the tables.",
            ],
            x=0.04,
            y=next_y - 0.020,
            width=0.92,
        )
        _draw_basic_table(
            ax,
            [
                ["cross_dataset_predictor_summary.tsv", "Rankings and coverage diagnostics"],
                ["common_prediction_summary.tsv", "Common scored-gene overlap by dataset"],
                ["predictor_combination_summary.tsv", "Predictor and rank-mean combination metrics"],
                ["datasets/<dataset_id>/", "Joined tables, plots, and predictor reports"],
            ],
            columns=["Primary file", "Purpose"],
            col_widths=[0.43, 0.57],
            bbox=[0.04, 0.045, 0.92, 0.155],
            font_size=10.2,
            wrap_scale=135,
        )
        _save_page(pdf, fig)

        fig, ax = _new_page()
        _header(ax, "Cross-Dataset Predictor Summary", "Coverage-aware interpretation of mean performance across selected datasets.")
        if summary_df is not None and not summary_df.empty:
            next_y = _bullet_block(
                ax,
                "Key results",
                _key_result_bullets(summary_df),
                x=0.04,
                y=0.84,
                width=0.92,
            )
            _bullet_block(
                ax,
                "Coverage notes",
                _coverage_note_bullets(summary_df),
                x=0.04,
                y=next_y - 0.018,
                width=0.92,
            )
            _draw_summary_table(ax, summary_df, bbox=[0.035, 0.075, 0.93, 0.305])
        else:
            _text_block(ax, "No summary table", ["Cross-dataset summary table is unavailable for this run."], x=0.06, y=0.83, width=0.88)
        _save_page(pdf, fig)

        fig, ax = _new_page()
        _header(ax, "Dataset Inventory", "Evaluated and skipped experiments for this benchmark run.")
        if dataset_outputs:
            dataset_df = pd.DataFrame(_dataset_report_rows(dataset_outputs))
            _draw_basic_table(
                ax,
                dataset_df.values.tolist(),
                columns=dataset_df.columns.tolist(),
                col_widths=[0.35, 0.25, 0.15, 0.25],
                bbox=[0.06, 0.61, 0.88, 0.20],
                font_size=REPORT_TABLE_SIZE,
            )
        else:
            _text_block(
                ax,
                "Evaluated datasets",
                ["No datasets were evaluated."],
                x=0.06,
                y=0.80,
                width=0.88,
            )

        if skipped_datasets:
            skipped_rows = [
                [
                    item["dataset_id"],
                    item["mirna"],
                    item["reason"],
                ]
                for item in skipped_datasets
            ]
            _draw_basic_table(
                ax,
                skipped_rows,
                columns=["Skipped dataset", "miRNA", "Reason"],
                col_widths=[0.30, 0.22, 0.48],
                bbox=[0.06, 0.31, 0.88, 0.22],
                font_size=REPORT_TABLE_SIZE,
            )

        _bullet_block(
            ax,
            "Included figure families",
            [
                "Combination plots compare predictors and rank-mean combinations.",
                "Common-prediction tables report predictor overlap.",
                "Rank pages pair distributions with binned count views.",
                "Rank-count panels use log-scaled y-axes.",
                "Dataset reports include heatmaps, CDFs, and predictor diagnostics.",
            ],
            x=0.06,
            y=0.24,
            width=0.88,
        )
        _save_page(pdf, fig)

        _draw_common_prediction_page(pdf, combined_outputs)

        paired_rank_keys = set()
        paired_rank_keys |= _draw_rank_pair_page(pdf, combined_outputs, rank_type="local")
        paired_rank_keys |= _draw_rank_pair_page(pdf, combined_outputs, rank_type="global")

        pending_plot_items = []
        for key, title, caption, path in _plot_items(combined_outputs):
            if key in paired_rank_keys:
                continue
            pending_plot_items.append((key, title, caption, path))
            if len(pending_plot_items) == 2:
                _draw_plot_pair_page(pdf, pending_plot_items)
                pending_plot_items = []
        if pending_plot_items:
            _draw_plot_pair_page(pdf, pending_plot_items)

    return report_path
