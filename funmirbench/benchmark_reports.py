"""README and PDF reporting helpers for FuNmiRBench benchmark runs."""

from __future__ import annotations

import pathlib
import textwrap

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

from funmirbench.evaluate import (
    REPORT_PAGE_SIZE,
    describe_gt_rule,
    _publication_tool_ids,
)
from funmirbench.cross_dataset import write_cross_dataset_summaries, write_metric_tables


MIN_HEADLINE_COVERAGE = 0.10


def _relative_display_path(path, *, base_dir):
    resolved = pathlib.Path(path).expanduser().resolve()
    try:
        return str(resolved.relative_to(base_dir.resolve()))
    except ValueError:
        return str(resolved)


def _format_summary_value(value, *, percent=False):
    if pd.isna(value):
        return "NA"
    number = float(value)
    if percent:
        return f"{number:.1%}"
    return f"{number:.3f}"


def _gt_threshold_box_text(fdr_threshold, effect_threshold):
    if fdr_threshold is None:
        return f"no FDR\neffect > {float(effect_threshold)}"
    return f"FDR < {float(fdr_threshold)}\neffect > {float(effect_threshold)}"


def _init_run_layout(base_dir):
    datasets_dir = base_dir / "datasets"
    plots_dir = base_dir / "plots"
    combined_plots_dir = plots_dir / "combined"
    tables_dir = base_dir / "tables"
    per_experiment_tables_dir = tables_dir / "per_experiment"
    combined_tables_dir = tables_dir / "combined"
    for path in (
        datasets_dir,
        plots_dir,
        combined_plots_dir,
        tables_dir,
        per_experiment_tables_dir,
        combined_tables_dir,
    ):
        path.mkdir(parents=True, exist_ok=True)
    return {
        "datasets_dir": datasets_dir,
        "plots_dir": plots_dir,
        "combined_plots_dir": combined_plots_dir,
        "tables_dir": tables_dir,
        "per_experiment_tables_dir": per_experiment_tables_dir,
        "combined_tables_dir": combined_tables_dir,
    }


def _load_cross_dataset_summary(combined_outputs):
    summary_path = combined_outputs.get("tables", {}).get("cross_dataset_predictor_summary")
    if not summary_path:
        return None
    summary_file = pathlib.Path(summary_path)
    if not summary_file.is_file():
        return None
    summary_df = pd.read_csv(summary_file, sep="\t")
    if summary_df.empty:
        return summary_df
    return summary_df.sort_values(
        ["aps_mean", "auroc_mean", "positive_coverage_mean", "coverage_mean"],
        ascending=[False, False, False, False],
    ).reset_index(drop=True)


def _cross_dataset_markdown_table(summary_df):
    lines = [
        "| Predictor | Mean coverage | Mean positive coverage | Mean APS | Mean PR-AUC | Mean Spearman | Mean AUROC |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in summary_df.itertuples(index=False):
        lines.append(
            "| "
            + " | ".join(
                [
                    str(row.tool_id),
                    _format_summary_value(row.coverage_mean, percent=True),
                    _format_summary_value(row.positive_coverage_mean, percent=True),
                    _format_summary_value(row.aps_mean),
                    _format_summary_value(row.pr_auc_mean),
                    _format_summary_value(row.spearman_mean),
                    _format_summary_value(row.auroc_mean),
                ]
            )
            + " |"
        )
    return lines


def _combination_summary_report_df(summary_path, *, max_rows=12):
    if not summary_path:
        return pd.DataFrame()
    summary_path = pathlib.Path(summary_path)
    if not summary_path.is_file():
        return pd.DataFrame()
    summary_df = pd.read_csv(summary_path, sep="\t")
    if summary_df.empty:
        return pd.DataFrame()
    work = summary_df.replace([float("inf"), float("-inf")], pd.NA).dropna(
        subset=["combination_id", "combination_size", "positive_coverage_mean", "aps_mean"]
    )
    if work.empty:
        return pd.DataFrame()
    work = work.sort_values(
        ["aps_mean", "positive_coverage_mean", "coverage_mean"],
        ascending=[False, False, False],
    ).head(max_rows)
    display_df = work[
        [
            "combination_id",
            "combination_size",
            "coverage_mean",
            "positive_coverage_mean",
            "aps_mean",
            "pr_auc_mean",
            "auroc_mean",
        ]
    ].copy()
    display_df.columns = [
        "Predictor or combination",
        "Size",
        "Mean coverage",
        "Mean positive coverage",
        "Mean APS",
        "Mean PR-AUC",
        "Mean AUROC",
    ]
    display_df["Predictor or combination"] = display_df["Predictor or combination"].map(
        lambda value: str(value).replace("+", " + ")
    )
    for column in ["Mean coverage", "Mean positive coverage"]:
        display_df[column] = display_df[column].map(lambda value: _format_summary_value(value, percent=True))
    for column in ["Mean APS", "Mean PR-AUC", "Mean AUROC"]:
        display_df[column] = display_df[column].map(lambda value: _format_summary_value(value))
    return display_df


def _combination_markdown_table(display_df):
    lines = [
        "| Predictor or combination | Size | Mean coverage | Mean positive coverage | Mean APS | Mean PR-AUC | Mean AUROC |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in display_df.itertuples(index=False):
        lines.append(
            "| "
            + " | ".join(str(value) for value in row)
            + " |"
        )
    return lines


def _best_headline_predictor(summary_df, *, min_coverage=MIN_HEADLINE_COVERAGE):
    if summary_df is None or summary_df.empty:
        return None

    coverage = pd.to_numeric(summary_df["coverage_mean"], errors="coerce")
    eligible = summary_df[coverage >= float(min_coverage)].copy()
    if eligible.empty:
        return None
    return eligible.sort_values(["aps_mean", "auroc_mean"], ascending=False).iloc[0]


def _report_takeaways(summary_df):
    if summary_df is None or summary_df.empty:
        return []
    takeaways = []
    specs = [
        ("aps_mean", "Highest mean APS"),
        ("auroc_mean", "Highest mean AUROC"),
        ("positive_coverage_mean", "Highest mean positive coverage"),
        ("spearman_mean", "Highest mean Spearman"),
    ]
    for metric_col, label in specs:
        row = summary_df.sort_values(metric_col, ascending=False).iloc[0]
        formatted = _format_summary_value(
            row[metric_col],
            percent=metric_col.endswith("coverage_mean"),
        )
        takeaways.append(f"{label}: {row['tool_id']} ({formatted})")
    sparse_row = summary_df.sort_values("coverage_mean", ascending=True).iloc[0]
    takeaways.append(
        "Lowest mean overall coverage: "
        f"{sparse_row['tool_id']} ({_format_summary_value(sparse_row['coverage_mean'], percent=True)})"
    )
    return takeaways


def _coverage_analysis_lines(summary_df):
    if summary_df is None or summary_df.empty:
        return ["Cross-dataset predictor summary is unavailable for this run."]
    lines = []
    sparse = summary_df[summary_df["coverage_mean"].astype(float) < 0.25].copy()
    if not sparse.empty:
        sparse_labels = [
            f"{row.tool_id} ({_format_summary_value(row.coverage_mean, percent=True)} coverage; "
            f"{int(getattr(row, 'aps_count', 0))} evaluated datasets)"
            for row in sparse.itertuples(index=False)
        ]
        lines.append("Sparse predictors: " + ", ".join(sparse_labels) + ". Treat their metrics as subset-specific.")
    coverage_gap = (
        summary_df["positive_coverage_mean"].astype(float)
        - summary_df["coverage_mean"].astype(float)
    )
    if coverage_gap.notna().any():
        row = summary_df.loc[coverage_gap.idxmax()]
        lines.append(
            "Largest positive-coverage enrichment: "
            f"{row['tool_id']} scores positives at "
            f"{_format_summary_value(row['positive_coverage_mean'], percent=True)} versus "
            f"{_format_summary_value(row['coverage_mean'], percent=True)} overall coverage."
        )
    best = _best_headline_predictor(summary_df)
    if best is not None:
        lines.append(
            f"Best mean APS with >= {MIN_HEADLINE_COVERAGE:.0%} mean coverage: "
            f"{best['tool_id']} ({_format_summary_value(best['aps_mean'])}; "
            f"AUROC {_format_summary_value(best['auroc_mean'])})."
        )
    else:
        lines.append(
            f"No predictor reached >= {MIN_HEADLINE_COVERAGE:.0%} mean coverage for headline ranking."
        )
    return lines


def write_run_readme(
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
    protein_coding_filter=None,
):
    summary_df = _load_cross_dataset_summary(combined_outputs)
    display_tool_ids = _publication_tool_ids(tool_ids) or list(tool_ids)
    relative_metric_tables = {
        key: _relative_display_path(path, base_dir=out_dir)
        for key, path in metric_tables.items()
    }
    relative_combined_outputs = {
        section: {
            key: _relative_display_path(path, base_dir=out_dir)
            for key, path in values.items()
        }
        for section, values in combined_outputs.items()
    }
    lines = [
        "# FuNmiRBench Run README",
        "",
        "## Summary",
        f"- Config: `{config_path}`",
        f"- Run directory: `{out_dir}`",
        f"- Datasets: `{len(dataset_outputs)}`",
        f"- Predictors: `{', '.join(display_tool_ids)}`",
        "",
        "## Evaluation Settings",
        f"- GT positive threshold: {describe_gt_rule(fdr_threshold, effect_threshold, markdown=True)}",
        (
            f"- Predictor agreement top fraction: `{predictor_top_fraction:.0%}`"
            " (exact top-k per predictor, deterministic tie-break)"
        ),
        "- Score handling: predictors are first aligned so that higher always means stronger",
        "- Per-dataset heatmaps and agreement plots: dataset-local tie-aware dense ranking over scored rows",
        "- Cross-dataset rank-distribution plots: global tie-aware dense ranking over each predictor's full standardized file",
        "- Combined PR/ROC/GSEA comparison plots: computed on the common set of genes scored by all compared predictors",
    ]
    if protein_coding_filter and protein_coding_filter.get("enabled"):
        lines.extend(
            [
                (
                    "- Protein-coding filter: enabled "
                    f"(Ensembl release {protein_coding_filter.get('ensembl_release', 'NA')}; "
                    f"{protein_coding_filter.get('gene_count', 'NA')} genes)"
                ),
                f"- Protein-coding GTF: `{protein_coding_filter.get('gtf_path', 'NA')}`",
                f"- Protein-coding gene cache: `{protein_coding_filter.get('gene_cache', 'NA')}`",
            ]
        )
    else:
        lines.append("- Protein-coding filter: disabled")
    lines.extend(
        [
            "",
            "## Cross-Dataset Summary",
            (
                "Coverage, positive coverage, and mean metric performance are summarized numerically below. "
                "The PDF report carries the same summary and analysis text, so coverage-vs-performance is reported "
                "as a table instead of a separate scatter plot."
            ),
            "",
        ]
    )
    if summary_df is not None and not summary_df.empty:
        lines.extend(_cross_dataset_markdown_table(summary_df))
        lines.extend(["", "### Coverage Notes"])
        lines.extend([f"- {line}" for line in _coverage_analysis_lines(summary_df)])
    else:
        lines.append("Cross-dataset summary table is unavailable for this run.")
    combination_table = _combination_summary_report_df(
        combined_outputs.get("tables", {}).get("predictor_combination_summary")
    )
    if not combination_table.empty:
        lines.extend(
            [
                "",
                "## Predictor-Combination Summary",
                "Top single predictors and rank-mean combinations by mean APS. "
                "The complete table is saved at `tables/combined/predictor_combination_summary.tsv`.",
                "",
            ]
        )
        lines.extend(_combination_markdown_table(combination_table))
    lines.extend(
        [
            "",
            "## Where To Start",
            "- Read `REPORT.pdf` for the main run-level summary with combined plots and explanatory notes",
            "- Open `tables/combined/cross_dataset_predictor_summary.tsv` for the compact numeric cross-dataset summary",
            "- Browse `plots/combined/` for cross-dataset comparison figures",
            "- Browse `datasets/<dataset_id>/` for dataset-specific tables, plots, and reports",
            "",
            "## Layout",
            "- `datasets/<dataset_id>/joined.tsv`: joined DE + predictor score table for one dataset",
            "- `datasets/<dataset_id>/plots/predictors/<tool_id>/`: per-tool visual outputs",
            "- `datasets/<dataset_id>/plots/comparisons/`: multi-predictor comparison plots for that dataset",
            "- `datasets/<dataset_id>/plots/heatmaps/`: dataset-level heatmaps",
            "- `datasets/<dataset_id>/reports/`: per-dataset Markdown/PDF reports and correlation TSVs",
            "- `tables/per_experiment/`: metric tables across datasets, one row per experiment",
            "- `tables/combined/`: cross-dataset predictor summary table",
            "- `plots/combined/metrics/`, `plots/combined/ranks/`, `plots/combined/combinations/`: cross-dataset comparison plots grouped by theme",
            "- `summary.json`: machine-readable run summary",
            "",
            "## Datasets",
            "| Dataset | miRNA | Perturbation | Cell line | Joined table |",
            "| --- | --- | --- | --- | --- |",
        ]
    )
    for item in dataset_outputs:
        lines.append(
            "| "
            + " | ".join(
                [
                    f"`{item['dataset_id']}`",
                    f"`{item['mirna']}`",
                    f"`{item['perturbation']}`",
                    f"`{item['cell_line']}`",
                    f"`{_relative_display_path(item['joined_tsv'], base_dir=out_dir)}`",
                ]
            )
            + " |"
        )
    lines.extend(
        [
            "",
            "## Key Files",
            "- `REPORT.pdf`: polished run-level summary with the main cross-dataset table and selected plots",
            "- `tables/combined/cross_dataset_predictor_summary.tsv`: exact numeric cross-dataset summary used in the report",
            "",
            "### Per-Experiment Tables",
        ]
    )
    metric_descriptions = {
        "coverage": "fraction of genes with predictor scores",
        "positive_coverage": "fraction of GT-positive genes that were scored",
        "aps": "average precision score per dataset",
        "pr_auc": "precision-recall AUC per dataset",
        "spearman": "score vs expected-effect Spearman correlation per dataset",
        "auroc": "AUROC per dataset",
    }
    for key, path in relative_metric_tables.items():
        lines.append(f"- `{path}`: {metric_descriptions.get(key, key)}")
    lines.extend(
        [
            "",
            "### Combined Plots",
        ]
    )
    combined_plot_descriptions = {
        "positive_background_local_rank_distributions": (
            "dataset-local rank separation of GT positives from background genes"
        ),
        "positive_background_local_rank_counts": (
            "dataset-local rank counts of GT positives and background genes on a log scale"
        ),
        "positive_background_global_rank_distributions": (
            "predictor-global rank separation of GT positives from background genes"
        ),
        "positive_background_global_rank_counts": (
            "predictor-global rank counts of GT positives and background genes on a log scale"
        ),
        "positive_recovery_fraction_by_prediction_count": (
            "GT-positive recovery fraction with endpoint labels"
        ),
        "predictor_combination_expanded_frontier": (
            "coverage/APS frontier for predictors and their combinations"
        ),
    }
    for key, path in relative_combined_outputs.get("plots", {}).items():
        if key.startswith("cross_dataset_") and key.endswith("_distribution"):
            metric_name = key[len("cross_dataset_") : -len("_distribution")]
            lines.append(f"- `{path}`: cross-dataset distribution of `{metric_name}` across the selected datasets")
        else:
            lines.append(f"- `{path}`: {combined_plot_descriptions.get(key, key)}")
    lines.extend(
        [
            "",
        ]
    )
    readme_path = out_dir / "README.md"
    readme_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return readme_path


def write_run_pdf_report(
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
):
    del metric_tables
    summary_df = _load_cross_dataset_summary(combined_outputs)
    report_path = out_dir / "REPORT.pdf"

    with PdfPages(report_path) as pdf:
        def new_page():
            page_fig = plt.figure(figsize=REPORT_PAGE_SIZE)
            page_ax = page_fig.add_axes([0.0, 0.0, 1.0, 1.0])
            page_ax.axis("off")
            page_fig.patch.set_facecolor("white")
            return page_fig, page_ax

        def save_page(fig):
            fig.set_size_inches(*REPORT_PAGE_SIZE, forward=True)
            fig.patch.set_facecolor("white")
            pdf.savefig(fig, facecolor="white")
            plt.close(fig)

        def add_header(ax, title, subtitle=None):
            ax.text(
                0.06,
                0.95,
                title,
                fontsize=20,
                fontweight="bold",
                color="#17324D",
                va="top",
                ha="left",
                family="DejaVu Sans",
            )
            if subtitle:
                ax.text(
                    0.06,
                    0.915,
                    subtitle,
                    fontsize=10.5,
                    color="#5B6577",
                    va="top",
                    ha="left",
                    family="DejaVu Sans",
                )
            ax.add_line(
                plt.Line2D([0.06, 0.94], [0.895, 0.895], color="#D8DEE9", linewidth=1.4)
            )

        def add_block(ax, title, lines, *, x, y, width):
            ax.text(
                x,
                y,
                title,
                fontsize=11.5,
                fontweight="bold",
                color="#2F5D8C",
                va="top",
                ha="left",
                family="DejaVu Sans",
            )
            current_y = y - 0.03
            for line in lines:
                wrapped = textwrap.wrap(line, width=max(26, int(width * 95))) or [""]
                for chunk in wrapped:
                    ax.text(
                        x,
                        current_y,
                        chunk,
                        fontsize=9.5,
                        color="#22303C",
                        va="top",
                        ha="left",
                        family="DejaVu Sans",
                    )
                    current_y -= 0.024
                current_y -= 0.004
            return current_y

        fig, ax = new_page()
        add_header(
            ax,
            "FuNmiRBench Run Report",
            f"{len(dataset_outputs)} datasets | {len(tool_ids)} predictors | generated from {config_path.name}",
        )
        summary_boxes = [
            ("Datasets", str(len(dataset_outputs))),
            ("Predictors", str(len(tool_ids))),
            ("GT Threshold", _gt_threshold_box_text(fdr_threshold, effect_threshold)),
            ("Top Fraction", f"{predictor_top_fraction:.0%}\nexact top-k"),
        ]
        x_positions = [0.06, 0.29, 0.52, 0.75]
        for (label, value), x in zip(summary_boxes, x_positions):
            ax.text(
                x,
                0.84,
                f"{label}\n{value}",
                fontsize=10.5,
                fontweight="bold",
                color="#17324D",
                va="top",
                ha="left",
                family="DejaVu Sans",
                bbox={
                    "boxstyle": "round,pad=0.45",
                    "facecolor": "#F5F8FC",
                    "edgecolor": "#D8E2EF",
                },
            )
        left_bottom = add_block(
            ax,
            "Evaluation Settings",
            [
                f"GT positives: {describe_gt_rule(fdr_threshold, effect_threshold)}",
                "Predictor scores are aligned so that higher always means stronger before evaluation.",
                "Per-dataset heatmaps and agreement plots use dataset-local tie-aware dense ranks.",
                "Combined PR/ROC/GSEA plots use only the common set of genes scored by all compared predictors.",
                "Top-prediction effect CDFs are optional diagnostics and are not written by default.",
            ],
            x=0.06,
            y=0.69,
            width=0.40,
        )
        right_bottom = add_block(
            ax,
            "What This Report Emphasizes",
            [
                "Cross-dataset coverage and performance are summarized numerically in the predictor table on the next page.",
                "The exact numeric source for the report is tables/combined/cross_dataset_predictor_summary.tsv.",
                "Sparse predictors can be skipped for datasets where their scored subset lacks positives or background genes.",
            ],
            x=0.54,
            y=0.69,
            width=0.38,
        )
        add_block(
            ax,
            "Key Files",
            [
                "REPORT.pdf: this run-level summary",
                "tables/combined/cross_dataset_predictor_summary.tsv: exact mean/median/std values by predictor",
                "tables/per_experiment/: per-dataset APS, PR-AUC, AUROC, coverage, positive coverage, and Spearman tables",
                "datasets/<dataset_id>/: joined tables, plots, and per-predictor reports",
            ],
            x=0.06,
            y=min(left_bottom, right_bottom) - 0.04,
            width=0.86,
        )
        save_page(fig)

        fig, ax = new_page()
        add_header(ax, "Cross-Dataset Summary", "Mean values are across the selected datasets only.")
        current_y = 0.85
        takeaways = _report_takeaways(summary_df)
        if takeaways:
            current_y = add_block(ax, "Quick Takeaways", takeaways, x=0.06, y=current_y, width=0.88) - 0.02
        coverage_lines = _coverage_analysis_lines(summary_df)
        if coverage_lines:
            current_y = add_block(ax, "Coverage Analysis", coverage_lines, x=0.06, y=current_y, width=0.88) - 0.03
        if summary_df is not None and not summary_df.empty:
            display_df = summary_df[
                [
                    "tool_id",
                    "coverage_mean",
                    "positive_coverage_mean",
                    "aps_mean",
                    "pr_auc_mean",
                    "spearman_mean",
                    "auroc_mean",
                ]
            ].copy()
            display_df.columns = [
                "Predictor",
                "Mean coverage",
                "Mean positive coverage",
                "Mean APS",
                "Mean PR-AUC",
                "Mean Spearman",
                "Mean AUROC",
            ]
            for column in display_df.columns[1:]:
                percent = "coverage" in column.lower()
                display_df[column] = display_df[column].map(
                    lambda value: _format_summary_value(value, percent=percent)
                )
            table_top = min(0.54, current_y)
            table_bottom = 0.18
            table = ax.table(
                cellText=display_df.values.tolist(),
                colLabels=display_df.columns.tolist(),
                cellLoc="center",
                colLoc="center",
                bbox=[0.04, table_bottom, 0.92, max(0.24, table_top - table_bottom)],
            )
            table.auto_set_font_size(False)
            table.set_fontsize(8.2)
            table.scale(1.0, 1.28)
            for (row, col), cell in table.get_celld().items():
                if row == 0:
                    cell.set_facecolor("#E9F1FB")
                    cell.set_edgecolor("#D8E2EF")
                    cell.set_text_props(weight="bold", color="#17324D")
                else:
                    cell.set_edgecolor("#E1E8F0")
                    cell.set_facecolor("#FFFFFF" if row % 2 else "#F9FBFD")
            footer = (
                "Coverage and positive coverage are table-first diagnostics. "
                "This avoids over-reading a scatter when sparse predictors are evaluated on very different subsets. "
                "Use the TSV for count/median/std/min/max."
            )
            for i, line in enumerate(textwrap.wrap(footer, width=112)):
                ax.text(
                    0.06,
                    0.125 - i * 0.024,
                    line,
                    fontsize=8.8,
                    color="#22303C",
                    va="top",
                    ha="left",
                    family="DejaVu Sans",
                )
        save_page(fig)

        fig, ax = new_page()
        add_header(ax, "Dataset Inventory", "Datasets included in this benchmark run.")
        dataset_df = pd.DataFrame(
            [
                {
                    "Dataset": item["dataset_id"],
                    "miRNA": item["mirna"],
                    "Perturbation": item["perturbation"],
                    "Cell line": item["cell_line"],
                }
                for item in dataset_outputs
            ]
        )
        table = ax.table(
            cellText=dataset_df.values.tolist(),
            colLabels=dataset_df.columns.tolist(),
            cellLoc="left",
            colLoc="left",
            bbox=[0.06, 0.56, 0.88, 0.24],
        )
        table.auto_set_font_size(False)
        table.set_fontsize(9.2)
        table.scale(1.0, 1.45)
        for (row, col), cell in table.get_celld().items():
            if row == 0:
                cell.set_facecolor("#E9F1FB")
                cell.set_edgecolor("#D8E2EF")
                cell.set_text_props(weight="bold", color="#17324D")
            else:
                cell.set_edgecolor("#E1E8F0")
                cell.set_facecolor("#FFFFFF" if row % 2 else "#F9FBFD")
        add_block(
            ax,
            "Included Combined Figures",
            [
                "plots/combined/metrics/cross_dataset_<metric>_distribution.png: one figure per metric showing how that metric varies across the selected datasets",
                "positive_background_local_rank_distributions.png: whether positives rank above background on the dataset-local rank scale",
                "positive_background_local_rank_counts.png: binned dataset-local counts for positives and background genes on a log scale",
                "positive_background_global_rank_distributions.png: whether positives rank above background on the predictor-global rank scale",
                "positive_background_global_rank_counts.png: binned predictor-global counts for positives and background genes on a log scale",
                "positive_recovery_fraction_by_prediction_count.png: GT-positive recovery fraction with endpoint labels",
                "predictor_combination_expanded_frontier.png: essential coverage/APS candidates for the best singles and combinations",
            ],
            x=0.06,
            y=0.45,
            width=0.88,
        )
        save_page(fig)

        combination_table = _combination_summary_report_df(
            combined_outputs.get("tables", {}).get("predictor_combination_summary")
        )
        if not combination_table.empty:
            fig, ax = new_page()
            add_header(
                ax,
                "Predictor-Combination Summary",
                "Top single predictors and rank-mean combinations by mean APS.",
            )
            table = ax.table(
                cellText=combination_table.values.tolist(),
                colLabels=combination_table.columns.tolist(),
                cellLoc="center",
                colLoc="center",
                bbox=[0.04, 0.24, 0.92, 0.56],
            )
            table.auto_set_font_size(False)
            table.set_fontsize(7.6)
            table.scale(1.0, 1.22)
            for (row, col), cell in table.get_celld().items():
                if row == 0:
                    cell.set_facecolor("#E9F1FB")
                    cell.set_edgecolor("#D8E2EF")
                    cell.set_text_props(weight="bold", color="#17324D")
                else:
                    cell.set_edgecolor("#E1E8F0")
                    cell.set_facecolor("#FFFFFF" if row % 2 else "#F9FBFD")
                    if col == 0:
                        cell.set_text_props(ha="left")
            footer = (
                "This page prints the leading rows from tables/combined/predictor_combination_summary.tsv. "
                "The TSV contains the complete set of evaluated combinations."
            )
            for i, line in enumerate(textwrap.wrap(footer, width=112)):
                ax.text(
                    0.06,
                    0.17 - i * 0.024,
                    line,
                    fontsize=8.8,
                    color="#22303C",
                    va="top",
                    ha="left",
                    family="DejaVu Sans",
                )
            save_page(fig)

        plot_items = []
        plot_descriptions = {}
        for metric_name in ["coverage", "positive_coverage", "aps", "pr_auc", "spearman", "auroc"]:
            plot_descriptions[f"cross_dataset_{metric_name}_distribution"] = (
                f"Cross-dataset {metric_name.upper()} distribution",
                f"Spread of {metric_name.upper()} across the selected datasets for every predictor.",
            )
        plot_descriptions.update({
            "positive_background_local_rank_distributions": (
                "Positive vs background local rank distributions",
                "Dataset-local rank distributions aggregated across datasets, split into GT positives and background genes. "
                "Stronger predictors should push positives higher than background."
            ),
            "positive_background_local_rank_counts": (
                "Positive vs background local rank counts",
                "Binned dataset-local rank counts for background genes and GT positives. "
                "The log y-axis keeps smaller GT-positive counts readable next to large background counts."
            ),
            "positive_background_global_rank_distributions": (
                "Positive vs background global rank distributions",
                "Predictor-global rank distributions aggregated across datasets, split into GT positives and background genes. "
                "This keeps each predictor on the rank scale of its full standardized file."
            ),
            "positive_background_global_rank_counts": (
                "Positive vs background global rank counts",
                "Binned predictor-global rank counts for background genes and GT positives. "
                "The log y-axis keeps smaller GT-positive counts readable next to large background counts."
            ),
            "positive_recovery_fraction_by_prediction_count": (
                "GT-positive recovery fraction by prediction count",
                "Cumulative curves showing the mean fraction of GT-positive genes recovered; endpoint labels show recovery at 300 predictions per dataset.",
            ),
            "predictor_combination_expanded_frontier": (
                "Predictor-combination essential candidates",
                "Coverage-versus-APS comparison of the best single predictors and best combinations; the full summary is printed in the preceding table and saved as TSV.",
            ),
        })
        for key, (title, caption) in plot_descriptions.items():
            path = combined_outputs.get("plots", {}).get(key)
            if path:
                plot_items.append((title, caption, pathlib.Path(path)))

        for title, caption, path in plot_items:
            if not path.is_file():
                continue
            image = plt.imread(path)
            fig, _ = new_page()
            fig.text(
                0.06,
                0.975,
                title,
                fontsize=12,
                fontweight="bold",
                color="#17324D",
                va="top",
                ha="left",
            )
            fig.text(
                0.06,
                0.94,
                caption,
                fontsize=9.4,
                color="#22303C",
                va="top",
                ha="left",
            )
            image_ax = fig.add_axes([0.04, 0.04, 0.92, 0.82])
            image_ax.imshow(image)
            image_ax.axis("off")
            save_page(fig)

    return report_path


def finalize_run_bundle(
    out_dir,
    *,
    out_root,
    config_path,
    dataset_outputs,
    tool_ids,
    metric_rows,
    joined_frames,
    tool_labels,
    fdr_threshold,
    effect_threshold,
    predictor_top_fraction,
    logger_info,
):
    layout = _init_run_layout(out_dir)
    metric_tables = write_metric_tables(
        metric_rows,
        layout["per_experiment_tables_dir"],
        logger=logger_info,
    )
    combined_outputs = write_cross_dataset_summaries(
        metric_rows,
        layout["combined_tables_dir"],
        layout["combined_plots_dir"],
        joined_frames=joined_frames,
        fdr_threshold=fdr_threshold,
        effect_threshold=effect_threshold,
        predictor_top_fraction=predictor_top_fraction,
        tool_labels=tool_labels,
        logger=logger_info,
    )
    readme_path = write_run_readme(
        out_dir,
        config_path=config_path,
        dataset_outputs=dataset_outputs,
        tool_ids=tool_ids,
        metric_tables=metric_tables,
        combined_outputs=combined_outputs,
        fdr_threshold=fdr_threshold,
        effect_threshold=effect_threshold,
        predictor_top_fraction=predictor_top_fraction,
    )
    report_path = write_run_pdf_report(
        out_dir,
        config_path=config_path,
        dataset_outputs=dataset_outputs,
        tool_ids=tool_ids,
        metric_tables=metric_tables,
        combined_outputs=combined_outputs,
        fdr_threshold=fdr_threshold,
        effect_threshold=effect_threshold,
        predictor_top_fraction=predictor_top_fraction,
    )
    summary = {
        "config": str(config_path),
        "out_root": str(out_root),
        "out_dir": str(out_dir),
        "run_dir_name": out_dir.name,
        "dataset_ids": [item["dataset_id"] for item in dataset_outputs],
        "tool_ids": tool_ids,
        "readme": str(readme_path),
        "report_pdf": str(report_path),
        "metric_tables": metric_tables,
        "cross_dataset_outputs": combined_outputs,
        "datasets": dataset_outputs,
    }
    summary_path = out_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    logger_info(f"Wrote summary: {summary_path}")
    return {
        "metric_tables": metric_tables,
        "combined_outputs": combined_outputs,
        "readme_path": readme_path,
        "report_path": report_path,
        "summary": summary,
        "summary_path": summary_path,
    }
