"""Cross-dataset metric tables and summary plots."""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import NullLocator, PercentFormatter

import funmirbench.evaluate as ev


CONTROL_PLOT_TOOL_IDS = set()
SUMMARY_STATS = ["count", "mean", "median", "std", "min", "max"]


def _is_real_plot_tool(tool_id):
    return ev._is_publication_tool(tool_id)


def _publication_metrics_df(metrics_df):
    display_df = metrics_df[metrics_df["tool_id"].map(_is_real_plot_tool)].copy()
    return display_df if not display_df.empty else metrics_df.copy()


def write_metric_tables(metric_rows, tables_dir, *, logger=None):
    tables_dir.mkdir(parents=True, exist_ok=True)
    metrics_df = pd.DataFrame(metric_rows)
    if metrics_df.empty:
        raise ValueError("No metric rows were produced.")
    metrics_df = _publication_metrics_df(metrics_df)

    id_cols = ["dataset_id", "mirna", "cell_line", "perturbation", "geo_accession"]
    metrics_df[id_cols] = metrics_df[id_cols].fillna("NA")
    out_paths = {}
    for metric_name, filename in [
        ("coverage", "coverage_per_experiment.tsv"),
        ("positive_coverage", "positive_coverage_per_experiment.tsv"),
        ("aps", "aps_per_experiment.tsv"),
        ("pr_auc", "pr_auc_per_experiment.tsv"),
        ("spearman", "spearman_per_experiment.tsv"),
        ("auroc", "auroc_per_experiment.tsv"),
    ]:
        wide = metrics_df.pivot_table(
            index=id_cols, columns="tool_id", values=metric_name, aggfunc="first",
        ).reset_index()
        wide.columns.name = None
        out_path = tables_dir / filename
        wide.to_csv(out_path, sep="\t", index=False)
        out_paths[metric_name] = str(out_path)
        ev._emit_log(logger, f"  Wrote {metric_name} table: {out_path}")
    return out_paths


def _summarize_by_tool(metrics_df, metric_names, *, group_cols=("tool_id",)):
    summary = metrics_df.groupby(list(group_cols))[metric_names].agg(SUMMARY_STATS)
    summary.columns = [f"{metric_name}_{stat_name}" for metric_name, stat_name in summary.columns]
    return summary.reset_index()


def write_continuous_metric_table(continuous_rows, tables_dir, *, logger=None):
    """Write per dataset x predictor continuous metrics, including class-NA rows."""
    tables_dir.mkdir(parents=True, exist_ok=True)
    out_path = tables_dir / "continuous_metrics_per_experiment.tsv"
    pd.DataFrame(continuous_rows).to_csv(out_path, sep="\t", index=False)
    ev._emit_log(logger, f"  Wrote continuous metric table: {out_path}")
    return str(out_path)


def write_reference_baseline_tables(reference_rows, per_experiment_dir, combined_dir, *, logger=None):
    """Write reference-baseline rows and their summary, kept out of the publication tables."""
    per_experiment_dir.mkdir(parents=True, exist_ok=True)
    combined_dir.mkdir(parents=True, exist_ok=True)
    reference_df = pd.DataFrame(reference_rows)
    per_experiment_path = per_experiment_dir / "reference_baselines_per_experiment.tsv"
    reference_df.to_csv(per_experiment_path, sep="\t", index=False)
    summary = _summarize_by_tool(
        reference_df,
        ["coverage", "positive_coverage", "aps", "pr_auc", "spearman", "auroc"],
    )
    summary.insert(1, "is_reference", True)
    summary_path = combined_dir / "cross_dataset_reference_summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)
    ev._emit_log(logger, f"  Wrote reference baseline tables: {per_experiment_path}, {summary_path}")
    return {"per_experiment": str(per_experiment_path), "summary": str(summary_path)}


def _plot_cross_dataset_metric_distribution(metrics_df, *, metric_name, out_path):
    metrics_df = _publication_metrics_df(metrics_df)
    if metrics_df.empty:
        return False
    fig, ax = plt.subplots(figsize=(8.6, 5.6))
    tool_ids = list(metrics_df["tool_id"].drop_duplicates())
    positions = np.arange(len(tool_ids), dtype=float)
    ev._style_axes(ax, grid_axis="y")
    data = []
    for tool_id in tool_ids:
        values = metrics_df.loc[metrics_df["tool_id"] == tool_id, metric_name].dropna().astype(float).tolist()
        data.append(values)
    box = ax.boxplot(
        data,
        positions=positions,
        widths=0.55,
        patch_artist=True,
        showfliers=False,
    )
    for patch, tool_id in zip(box["boxes"], tool_ids):
        color = ev._tool_color(tool_id)
        patch.set_facecolor(color)
        patch.set_alpha(0.35)
        patch.set_edgecolor(color)
    for median in box["medians"]:
        median.set_color("#22303C")
        median.set_linewidth(1.4)
    for whisker in box["whiskers"]:
        whisker.set_color("#7A8798")
    for cap in box["caps"]:
        cap.set_color("#7A8798")

    for tool_index, values in enumerate(data):
        if not values:
            continue
        jitter = np.linspace(-0.09, 0.09, num=len(values)) if len(values) > 1 else np.array([0.0])
        ax.scatter(
            np.full(len(values), positions[tool_index]) + jitter,
            values,
            s=30,
            alpha=0.75,
            color=ev._tool_color(tool_ids[tool_index]),
            edgecolors="white",
            linewidths=0.45,
            zorder=3,
        )

    all_values = [value for values in data for value in values if np.isfinite(value)]
    if all_values:
        ymin, ymax = _metric_data_limits(all_values, metric_name)
    else:
        ymin, ymax = ev._metric_plot_limits(metric_name)
    ax.set_ylim(ymin, ymax)
    ax.set_ylabel(metric_name.upper(), fontsize=ev.PLOT_AXIS_LABEL_SIZE)
    ax.set_xticks(positions)
    ax.set_xticklabels(
        [ev._wrap_axis_label(ev._tool_label(tool_id), width=13) for tool_id in tool_ids],
        rotation=0,
        ha="center",
    )
    ev._add_figure_heading(
        fig,
        title=f"Cross-dataset {metric_name.upper()} distribution",
        subtitle="Each box summarizes per-dataset metric values.",
    )
    ev._save_figure(fig, out_path)
    return True


def _metric_data_limits(values, metric_name):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return ev._metric_plot_limits(metric_name)
    min_value = float(values.min())
    max_value = float(values.max())
    span = max(max_value - min_value, 0.04)
    pad = span * 0.22
    lower = min_value - pad
    upper = max_value + pad
    if metric_name in {"coverage", "positive_coverage", "aps", "pr_auc", "auroc"}:
        lower = max(0.0, lower)
        upper = min(1.02, upper)
    elif metric_name == "spearman":
        lower = max(-1.02, lower)
        upper = min(1.02, upper)
    if upper <= lower:
        upper = lower + 0.05
    return lower, upper


def _rank_distribution_metadata(rank_type):
    if rank_type == "local":
        return {
            "title": "Positive vs background local rank distributions",
            "axis_label": "Local rank within dataset",
            "subtitle": (
                "Dense ranks are dataset-local; GT positives use predictor colors."
            ),
        }
    if rank_type == "global":
        return {
            "title": "Positive vs background global rank distributions",
            "axis_label": "Global rank across predictor file",
            "subtitle": (
                "Dense ranks use each full predictor file; GT positives use predictor colors."
            ),
        }
    raise ValueError(f"Unsupported rank distribution type: {rank_type}")


def _rank_class_plot_data(joined_frames, *, fdr_threshold, effect_threshold, rank_type):
    combined = pd.concat(joined_frames, ignore_index=True)
    rank_specs = ev._rank_distribution_specs(combined, rank_types=(rank_type,))
    rank_specs = [
        (tool_id, column, column_type)
        for tool_id, column, column_type in rank_specs
        if _is_real_plot_tool(tool_id)
    ]
    if not rank_specs:
        return None

    keep_cols = ["logFC", "FDR", *[column for _, column, _ in rank_specs]]
    for optional in ("dataset_id", "perturbation", *ev.FDR_AUXILIARY_COLUMNS):
        if optional in combined.columns:
            keep_cols.append(optional)
    keep_cols = list(dict.fromkeys(keep_cols))
    work = ev._filter_usable_gt_rows(combined[keep_cols], fdr_threshold=fdr_threshold)
    if work.empty:
        return None
    work = ev._annotate_ground_truth(work)
    work["is_positive"] = ev._positive_mask(
        work,
        fdr_threshold=fdr_threshold,
        effect_threshold=effect_threshold,
    ).astype(int)

    tool_ids = []
    positive_data = []
    background_data = []
    for tool_id, column, column_type in rank_specs:
        if column_type == "score":
            values = ev._rank_scale_scores(work[column])
        else:
            values = work[column].astype(float)
        pos_values = values.loc[work["is_positive"] == 1].dropna().astype(float).tolist()
        bg_values = values.loc[work["is_positive"] == 0].dropna().astype(float).tolist()
        tool_ids.append(tool_id)
        positive_data.append(pos_values)
        background_data.append(bg_values)
    if not any(positive_data) or not any(background_data):
        return None
    return {
        "tool_ids": tool_ids,
        "positive_data": positive_data,
        "background_data": background_data,
    }


def _plot_rank_class_distributions(
    joined_frames, *, out_path, fdr_threshold, effect_threshold, rank_type
):
    plot_data = _rank_class_plot_data(
        joined_frames,
        fdr_threshold=fdr_threshold,
        effect_threshold=effect_threshold,
        rank_type=rank_type,
    )
    if plot_data is None:
        return False
    tool_ids = plot_data["tool_ids"]
    positive_data = plot_data["positive_data"]
    background_data = plot_data["background_data"]

    fig, ax = plt.subplots(figsize=(max(8.4, len(tool_ids) * 2.2), 6.0))
    ev._style_axes(ax, grid_axis="y")
    positions = np.arange(len(tool_ids), dtype=float) * 1.35
    negative_positions = positions - 0.18
    positive_positions = positions + 0.18
    colors = [ev._tool_color(tool_id) for tool_id in tool_ids]
    valid_bg_positions = [pos for pos, values in zip(negative_positions, background_data) if values]
    valid_bg_data = [values for values in background_data if values]
    valid_pos_positions = [pos for pos, values in zip(positive_positions, positive_data) if values]
    valid_pos_data = [values for values in positive_data if values]

    if valid_bg_data:
        bg_violin = ax.violinplot(
            valid_bg_data,
            positions=valid_bg_positions,
            widths=0.30,
            showmeans=False,
            showextrema=False,
            showmedians=True,
        )
        for body in bg_violin["bodies"]:
            body.set_facecolor(ev.NEGATIVE_COLOR)
            body.set_edgecolor(ev.NEGATIVE_COLOR)
            body.set_alpha(0.45)
        bg_violin["cmedians"].set_color("#22303C")
        bg_violin["cmedians"].set_linewidth(1.5)

    if valid_pos_data:
        pos_violin = ax.violinplot(
            valid_pos_data,
            positions=valid_pos_positions,
            widths=0.30,
            showmeans=False,
            showextrema=False,
            showmedians=True,
        )
        for body, color in zip(pos_violin["bodies"], [c for c, values in zip(colors, positive_data) if values]):
            body.set_facecolor(color)
            body.set_edgecolor(color)
            body.set_alpha(0.40)
        pos_violin["cmedians"].set_color("#22303C")
        pos_violin["cmedians"].set_linewidth(1.5)

    meta = _rank_distribution_metadata(rank_type)
    ax.set_ylim(0, 1.02)
    ax.set_ylabel(meta["axis_label"], fontsize=ev.PLOT_AXIS_LABEL_SIZE)
    ax.set_xticks(positions)
    ax.set_xticklabels(
        [ev._wrap_axis_label(ev._tool_label(tool_id), width=12) for tool_id in tool_ids],
        rotation=0,
        ha="center",
    )
    ev._add_figure_heading(
        fig,
        title=meta["title"],
        subtitle=meta["subtitle"],
    )
    ax.legend(
        handles=[
            Patch(facecolor=ev.NEGATIVE_COLOR, edgecolor=ev.NEGATIVE_COLOR, alpha=0.45, label="Background genes"),
            Patch(
                facecolor="#C8D6EA",
                edgecolor="#6E89A8",
                alpha=0.50,
                label="GT positives (predictor color)",
            ),
        ],
        frameon=False,
        fontsize=ev.PLOT_LEGEND_SIZE,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.22),
        ncol=2,
        borderaxespad=0.0,
    )
    fig.subplots_adjust(top=0.84, bottom=0.30, left=0.10, right=0.96)
    ev._save_figure(fig, out_path)
    return True


def _plot_rank_class_count_distributions(
    joined_frames, *, out_path, fdr_threshold, effect_threshold, rank_type
):
    plot_data = _rank_class_plot_data(
        joined_frames,
        fdr_threshold=fdr_threshold,
        effect_threshold=effect_threshold,
        rank_type=rank_type,
    )
    if plot_data is None:
        return False
    tool_ids = plot_data["tool_ids"]
    positive_data = plot_data["positive_data"]
    background_data = plot_data["background_data"]

    column_count = min(3, max(1, len(tool_ids)))
    row_count = int(np.ceil(len(tool_ids) / column_count))
    fig, axes = plt.subplots(
        row_count,
        column_count,
        figsize=(max(8.4, column_count * 3.1), row_count * 3.8 + 1.9),
        sharex=True,
        sharey=True,
    )
    axes = np.atleast_1d(axes).ravel()
    bins = np.linspace(0, 1, 21)
    max_count = 1
    for values in [*positive_data, *background_data]:
        if not values:
            continue
        counts, _ = np.histogram(values, bins=bins)
        max_count = max(max_count, int(counts.max()))

    meta = _rank_distribution_metadata(rank_type)
    for ax, tool_id, pos_values, bg_values in zip(
        axes,
        tool_ids,
        positive_data,
        background_data,
    ):
        ev._style_axes(ax, grid_axis="y")
        if bg_values:
            ax.hist(
                bg_values,
                bins=bins,
                color=ev.NEGATIVE_COLOR,
                alpha=0.55,
            )
        if pos_values:
            ax.hist(
                pos_values,
                bins=bins,
                histtype="step",
                color=ev._tool_color(tool_id),
                linewidth=2.2,
            )
        ax.set_title(ev._tool_label(tool_id), fontsize=12, loc="left", pad=8)
        ax.set_xlim(0, 1)
        ax.set_yscale("log")
        ax.yaxis.set_minor_locator(NullLocator())
        ax.set_ylim(0.8, max_count * 1.6)
        ax.tick_params(
            axis="both",
            which="major",
            bottom=True,
            left=True,
            labelbottom=True,
            labelleft=True,
        )

    for ax in axes[len(tool_ids) :]:
        ax.set_visible(False)

    ev._add_figure_heading(
        fig,
        title=meta["title"].replace("distributions", "counts"),
        subtitle="Genes are binned by normalized rank; y-axis shows log-scaled counts.",
        x=0.08,
        subtitle_y=0.91,
    )
    fig.supxlabel(meta["axis_label"], fontsize=ev.PLOT_AXIS_LABEL_SIZE, y=0.070)
    fig.supylabel("Gene count (log scale)", fontsize=ev.PLOT_AXIS_LABEL_SIZE, x=0.035)
    fig.subplots_adjust(top=0.70, bottom=0.24, hspace=0.48, wspace=0.24)
    legend_handles = [
        Patch(
            facecolor=ev.NEGATIVE_COLOR,
            edgecolor=ev.NEGATIVE_COLOR,
            alpha=0.55,
            label="Background genes",
        ),
        Line2D(
            [0],
            [0],
            color=ev.POSITIVE_COLOR,
            linewidth=2.2,
            label="GT positives (predictor color)",
        ),
    ]
    fig.legend(
        handles=legend_handles,
        frameon=False,
        fontsize=ev.PLOT_LEGEND_SIZE,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.82),
        ncol=2,
    )
    ev._save_figure(fig, out_path)
    return True


def _rank_values_for_spec(frame, column, column_type):
    if column_type == "score":
        return ev._rank_scale_scores(frame[column])
    return frame[column].astype(float)


def _recovery_percent_decimals(max_value):
    max_value = float(max_value or 0.0)
    if max_value < 0.01:
        return 2
    if max_value < 0.10:
        return 1
    return 0


def _recovery_fraction_y_max(max_endpoint):
    target = max(float(max_endpoint or 0.0) * 1.35, 0.005)
    for limit in (0.005, 0.01, 0.02, 0.05, 0.10, 0.20, 0.50, 1.02):
        if target <= limit:
            return limit
    return 1.02


def _add_recovery_endpoint_labels(ax, label_items, *, x_value, y_min, y_max, percent_decimals=0):
    if not label_items:
        return
    y_range = max(float(y_max) - float(y_min), 1e-9)
    min_gap = y_range * 0.035
    sorted_items = sorted(label_items, key=lambda item: item["y"])
    adjusted = []
    next_y = float(y_min)
    for item in sorted_items:
        y_text = min(max(float(item["y"]), float(y_min)), float(y_max))
        y_text = max(y_text, next_y)
        adjusted.append((item, y_text))
        next_y = y_text + min_gap
    if adjusted and adjusted[-1][1] > y_max:
        overflow = adjusted[-1][1] - float(y_max)
        adjusted = [(item, max(float(y_min), y_text - overflow)) for item, y_text in adjusted]

    label_x = float(x_value) + max(float(x_value) * 0.018, 1.5)
    for item, y_text in adjusted:
        label = f"{item['label']} ({item['value']:.{int(percent_decimals)}%})"
        ax.annotate(
            label,
            (x_value, item["y"]),
            xytext=(label_x, y_text),
            textcoords="data",
            ha="left",
            va="center",
            fontsize=ev.PLOT_LEGEND_SIZE,
            color=item["color"],
            annotation_clip=False,
            clip_on=False,
            arrowprops={
                "arrowstyle": "-",
                "color": item["color"],
                "linewidth": 0.75,
                "alpha": 0.8,
                "shrinkA": 2,
                "shrinkB": 2,
            },
        )


def _plot_positive_recovery_by_prediction_count(
    joined_frames,
    *,
    out_path,
    fdr_threshold,
    effect_threshold,
    max_predictions=300,
    normalized=False,
    endpoint_labels=False,
    excluded_tool_ids=None,
):
    combined = pd.concat(joined_frames, ignore_index=True)
    rank_specs = ev._rank_distribution_specs(combined, rank_types=("local", "global", "score"))
    if excluded_tool_ids:
        excluded_tool_ids = {str(tool_id) for tool_id in excluded_tool_ids}
        rank_specs = [
            (tool_id, column, column_type)
            for tool_id, column, column_type in rank_specs
            if str(tool_id) not in excluded_tool_ids
        ]
    if not rank_specs:
        return False

    dataset_col = "dataset_id" if "dataset_id" in combined.columns else None
    keep_cols = ["gene_id", "logFC", "FDR", *[column for _, column, _ in rank_specs]]
    for optional in ("dataset_id", "perturbation", *ev.FDR_AUXILIARY_COLUMNS):
        if optional in combined.columns:
            keep_cols.append(optional)
    keep_cols = list(dict.fromkeys(keep_cols))
    work = ev._filter_usable_gt_rows(combined[keep_cols], fdr_threshold=fdr_threshold)
    if work.empty:
        return False

    work = ev._annotate_ground_truth(work)
    work["is_positive"] = ev._positive_mask(
        work,
        fdr_threshold=fdr_threshold,
        effect_threshold=effect_threshold,
    ).astype(int)
    groups = (
        list(work.groupby(dataset_col, sort=True))
        if dataset_col
        else [("all_datasets", work)]
    )

    max_predictions = int(max_predictions)
    if max_predictions <= 0:
        return False
    x_values = np.arange(1, max_predictions + 1, dtype=int)
    plotted = False
    endpoint_items = []

    fig, ax = plt.subplots(figsize=(9.2 if endpoint_labels else 8.8, 5.8))
    ev._style_axes(ax, grid_axis="both")
    for tool_id, column, column_type in rank_specs:
        dataset_curves = []
        for _, dataset_frame in groups:
            positive_total = int(dataset_frame["is_positive"].sum())
            if normalized and positive_total <= 0:
                continue
            values = _rank_values_for_spec(dataset_frame, column, column_type)
            scored = dataset_frame.loc[values.notna(), ["gene_id", "is_positive"]].copy()
            if scored.empty:
                continue
            scored["rank_value"] = values.loc[scored.index].astype(float)
            scored = scored.sort_values(
                ["rank_value", "gene_id"],
                ascending=[False, True],
                kind="mergesort",
            )
            hits = scored["is_positive"].to_numpy(dtype=int)
            if hits.size == 0:
                continue
            cumulative = np.cumsum(hits).astype(float)
            if normalized:
                cumulative = cumulative / float(positive_total)
            curve = np.repeat(float(cumulative[-1]), max_predictions)
            observed = min(max_predictions, cumulative.size)
            curve[:observed] = cumulative[:observed]
            dataset_curves.append(curve)

        if not dataset_curves:
            continue
        mean_curve = np.vstack(dataset_curves).mean(axis=0)
        if not np.isfinite(mean_curve).any():
            continue
        plotted = True
        color = ev._tool_color(tool_id)
        ax.plot(
            x_values,
            mean_curve,
            linewidth=2.8,
            color=color,
            label=ev._tool_label(tool_id),
        )
        if endpoint_labels:
            endpoint_items.append(
                {
                    "label": ev._tool_label(tool_id),
                    "value": float(mean_curve[-1]),
                    "y": float(mean_curve[-1]),
                    "color": color,
                }
            )

    if not plotted:
        plt.close(fig)
        return False

    ax.set_xlim(1, max_predictions)
    if normalized:
        max_endpoint = max((item["y"] for item in endpoint_items), default=0.0)
        y_max = (
            1.02
            if not excluded_tool_ids
            else _recovery_fraction_y_max(max_endpoint)
        )
        ax.set_ylim(0, y_max)
        percent_decimals = _recovery_percent_decimals(y_max)
        ax.yaxis.set_major_formatter(PercentFormatter(xmax=1.0, decimals=percent_decimals))
    else:
        ax.set_ylim(bottom=0)
        percent_decimals = 0
    ax.set_xlabel("Predicted targets per dataset", fontsize=ev.PLOT_AXIS_LABEL_SIZE)
    ax.set_ylabel(
        "Mean GT-positive recovery fraction" if normalized else "Mean GT positives recovered",
        fontsize=ev.PLOT_AXIS_LABEL_SIZE,
    )
    if normalized:
        title = "GT-positive recovery fraction by prediction count"
        subtitle = (
            "Each curve shows the mean fraction of GT positives recovered; "
            f"endpoint labels show recovery at {max_predictions} predictions."
        )
    else:
        title = "GT-positive recovery by prediction count"
        subtitle = (
            "Each curve shows cumulative mean GT-positive genes recovered as top-ranked predictions are admitted per dataset."
        )
    ev._add_figure_heading(fig, title=title, subtitle=subtitle)
    if endpoint_labels:
        _add_recovery_endpoint_labels(
            ax,
            endpoint_items,
            x_value=max_predictions,
            y_min=ax.get_ylim()[0],
            y_max=ax.get_ylim()[1],
            percent_decimals=percent_decimals,
        )
    else:
        ax.legend(frameon=False, fontsize=ev.PLOT_LEGEND_SIZE, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    ev._save_figure(fig, out_path)
    return True


def write_cross_dataset_summaries(
    metric_rows,
    tables_dir,
    plots_dir,
    *,
    joined_frames=None,
    fdr_threshold=0.05,
    effect_threshold=1.0,
    predictor_top_fraction=0.10,
    tool_labels=None,
    logger=None,
):
    del predictor_top_fraction
    tables_dir.mkdir(parents=True, exist_ok=True)
    plots_dir.mkdir(parents=True, exist_ok=True)
    ev._set_tool_labels(tool_labels)
    metric_plots_dir = plots_dir / "metrics"
    rank_plots_dir = plots_dir / "ranks"
    for path in (metric_plots_dir, rank_plots_dir):
        path.mkdir(parents=True, exist_ok=True)
    metrics_df = pd.DataFrame(metric_rows)
    if metrics_df.empty:
        return {
            "tables": {},
            "plots": {},
        }
    metrics_df = _publication_metrics_df(metrics_df)
    real_tool_ids = [
        tool_id
        for tool_id in metrics_df["tool_id"].drop_duplicates().tolist()
        if _is_real_plot_tool(tool_id)
    ]
    ev._set_tool_colors(real_tool_ids)

    metric_names = ["coverage", "positive_coverage", "aps", "pr_auc", "spearman", "auroc"]
    summary = _summarize_by_tool(metrics_df, metric_names)
    summary_path = tables_dir / "cross_dataset_predictor_summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)
    ev._emit_log(logger, f"  Wrote cross-dataset summary table: {summary_path}")
    distribution_paths = {}
    for metric_name in metric_names:
        metric_path = metric_plots_dir / f"cross_dataset_{metric_name}_distribution.png"
        wrote_distribution = _plot_cross_dataset_metric_distribution(
            metrics_df,
            metric_name=metric_name,
            out_path=metric_path,
        )
        if wrote_distribution:
            distribution_paths[f"cross_dataset_{metric_name}_distribution"] = str(metric_path)
            ev._emit_log(logger, f"  Wrote cross-dataset {metric_name} distribution: {metric_path}")

    rank_distribution_paths = {}
    if joined_frames:
        for rank_type, plot_key, filename, count_plot_key, count_filename in (
            (
                "local",
                "positive_background_local_rank_distributions",
                "positive_background_local_rank_distributions.png",
                "positive_background_local_rank_counts",
                "positive_background_local_rank_counts.png",
            ),
            (
                "global",
                "positive_background_global_rank_distributions",
                "positive_background_global_rank_distributions.png",
                "positive_background_global_rank_counts",
                "positive_background_global_rank_counts.png",
            ),
        ):
            rank_distribution_path = rank_plots_dir / filename
            wrote_rank_distribution = _plot_rank_class_distributions(
                joined_frames,
                out_path=rank_distribution_path,
                fdr_threshold=fdr_threshold,
                effect_threshold=effect_threshold,
                rank_type=rank_type,
            )
            if wrote_rank_distribution:
                rank_distribution_paths[plot_key] = str(rank_distribution_path)
                ev._emit_log(logger, f"  Wrote {rank_type} rank distribution plot: {rank_distribution_path}")
            rank_count_path = rank_plots_dir / count_filename
            wrote_rank_count = _plot_rank_class_count_distributions(
                joined_frames,
                out_path=rank_count_path,
                fdr_threshold=fdr_threshold,
                effect_threshold=effect_threshold,
                rank_type=rank_type,
            )
            if wrote_rank_count:
                rank_distribution_paths[count_plot_key] = str(rank_count_path)
                ev._emit_log(logger, f"  Wrote {rank_type} rank count plot: {rank_count_path}")
        normalized_recovery_path = rank_plots_dir / "positive_recovery_fraction_by_prediction_count.png"
        wrote_normalized_recovery = _plot_positive_recovery_by_prediction_count(
            joined_frames,
            out_path=normalized_recovery_path,
            fdr_threshold=fdr_threshold,
            effect_threshold=effect_threshold,
            normalized=True,
            endpoint_labels=True,
            excluded_tool_ids=CONTROL_PLOT_TOOL_IDS,
        )
        if wrote_normalized_recovery:
            rank_distribution_paths["positive_recovery_fraction_by_prediction_count"] = str(
                normalized_recovery_path
            )
            ev._emit_log(
                logger,
                f"  Wrote normalized positive recovery plot: {normalized_recovery_path}",
            )

    return {
        "tables": {
            "cross_dataset_predictor_summary": str(summary_path),
        },
        "plots": {
            **distribution_paths,
            **rank_distribution_paths,
        },
    }
