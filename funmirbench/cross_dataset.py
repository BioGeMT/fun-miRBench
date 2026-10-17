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
MATCHED_POPULATION = "matched"
MATCHED_SUMMARY_METRICS = ["aps", "pr_auc", "auroc", "spearman", "pearson"]
MATCHED_COHORT_STATS = ["mean", "median", "std", "min", "max"]
MATCHED_METRIC_COLUMNS = [
    "dataset_id",
    "tool_id",
    "tool_type",
    "population",
    "matched_tool_ids",
    "common_rows",
    "common_positives",
    *MATCHED_SUMMARY_METRICS,
    "metrics_na_reason",
]
MATCHED_ROSTER_COLUMNS = ["tool_id", "tool_type", "configured", "included", "exclusion_reason"]
PAIRED_REFERENCE_METRICS = ["aps", "pr_auc", "auroc", "spearman", "pearson"]
PAIRED_REFERENCE_COLUMNS = [
    "tool_id",
    "reference_tool_id",
    "population",
    "summary_type",
    "metric",
    "comparisons",
    "tool_mean",
    "reference_mean",
    "mean_difference_reference_minus_tool",
    "median_difference_reference_minus_tool",
]
NO_METRIC_ROWS_REASON = "no metric rows in any evaluated dataset"


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


def _paired_reference_rows(metric_rows, reference_df):
    """Pair each predictor with its ``matched_to_<tool_id>`` references, metric by metric.

    For every metric, only datasets where both the predictor and the reference have a
    finite value are compared, so both means cover the same comparisons.
    """
    if not metric_rows or reference_df.empty or "matched_tool_id" not in reference_df.columns:
        return []
    tool_df = _publication_metrics_df(pd.DataFrame(metric_rows))
    matched = reference_df[reference_df["matched_tool_id"].fillna("").astype(str) != ""]
    rows = []
    for (tool_id, reference_id), reference_group in matched.groupby(["matched_tool_id", "tool_id"], sort=True):
        tool_metrics = tool_df[tool_df["tool_id"].astype(str) == str(tool_id)].set_index("dataset_id")
        reference_metrics = reference_group.set_index("dataset_id")
        for metric_name in PAIRED_REFERENCE_METRICS:
            if metric_name not in tool_metrics.columns or metric_name not in reference_metrics.columns:
                continue
            pairs = pd.concat(
                {"tool": tool_metrics[metric_name], "reference": reference_metrics[metric_name]},
                axis=1,
            ).astype(float)
            pairs = pairs[np.isfinite(pairs).all(axis=1)]
            difference = pairs["reference"] - pairs["tool"]
            rows.append(
                {
                    "tool_id": tool_id,
                    "reference_tool_id": reference_id,
                    "population": f"matched_to_{tool_id}",
                    "summary_type": "paired_common_cohort",
                    "metric": metric_name,
                    "comparisons": int(len(pairs)),
                    "tool_mean": float(pairs["tool"].mean()) if len(pairs) else float("nan"),
                    "reference_mean": float(pairs["reference"].mean()) if len(pairs) else float("nan"),
                    "mean_difference_reference_minus_tool": (
                        float(difference.mean()) if len(pairs) else float("nan")
                    ),
                    "median_difference_reference_minus_tool": (
                        float(difference.median()) if len(pairs) else float("nan")
                    ),
                }
            )
    return rows


def write_reference_baseline_tables(
    reference_rows, per_experiment_dir, combined_dir, *, metric_rows=None, logger=None,
):
    """Write reference-baseline rows, their descriptive summary, and paired comparisons.

    Reference baselines are kept out of the publication metric tables and summaries. The
    per-reference summary (``summary_type=descriptive_unpaired``) is grouped by population
    and averages each reference over its own datasets. Tool-vs-reference comparisons belong
    in the paired table, which uses a metric-specific common cohort for each pair.
    """
    per_experiment_dir.mkdir(parents=True, exist_ok=True)
    combined_dir.mkdir(parents=True, exist_ok=True)
    reference_df = pd.DataFrame(reference_rows)
    per_experiment_path = per_experiment_dir / "reference_baselines_per_experiment.tsv"
    reference_df.to_csv(per_experiment_path, sep="\t", index=False)
    summary = _summarize_by_tool(
        reference_df,
        ["coverage", "positive_coverage", "aps", "pr_auc", "spearman", "auroc"],
        group_cols=("population", "matched_tool_id", "tool_id"),
    )
    summary = summary[
        ["tool_id", "population", "matched_tool_id"]
        + [column for column in summary.columns if column not in {"tool_id", "population", "matched_tool_id"}]
    ]
    summary.insert(1, "is_reference", True)
    summary.insert(4, "summary_type", "descriptive_unpaired")
    summary_path = combined_dir / "cross_dataset_reference_summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)
    paired_path = combined_dir / "cross_dataset_tool_vs_reference_paired.tsv"
    pd.DataFrame(
        _paired_reference_rows(metric_rows, reference_df),
        columns=PAIRED_REFERENCE_COLUMNS,
    ).to_csv(paired_path, sep="\t", index=False)
    ev._emit_log(
        logger,
        f"  Wrote reference baseline tables: {per_experiment_path}, {summary_path}, {paired_path}",
    )
    return {
        "per_experiment": str(per_experiment_path),
        "summary": str(summary_path),
        "paired_summary": str(paired_path),
    }


def _matched_tool_roster(metrics_df, *, configured_tool_ids=None, reference_baselines=False):
    """Return the configured, included, and excluded tools of the matched population."""
    evaluated = set(metrics_df["tool_id"].astype(str))
    if configured_tool_ids is None:
        configured_tool_ids = metrics_df["tool_id"].astype(str).drop_duplicates().tolist()
    configured_tool_ids = [str(tool_id) for tool_id in configured_tool_ids]
    rows = []
    for tool_id in configured_tool_ids:
        included = tool_id in evaluated
        rows.append(
            {
                "tool_id": tool_id,
                "tool_type": "predictor",
                "configured": True,
                "included": included,
                "exclusion_reason": "" if included else NO_METRIC_ROWS_REASON,
            }
        )
    for tool_id in sorted(evaluated - set(configured_tool_ids)):
        rows.append(
            {
                "tool_id": tool_id,
                "tool_type": "predictor",
                "configured": False,
                "included": False,
                "exclusion_reason": "not in the configured predictor roster",
            }
        )
    if reference_baselines:
        for tool_id in sorted(ev.REFERENCE_TOOL_IDS):
            rows.append(
                {
                    "tool_id": tool_id,
                    "tool_type": "reference",
                    "configured": True,
                    "included": True,
                    "exclusion_reason": "",
                }
            )
    return pd.DataFrame(rows, columns=MATCHED_ROSTER_COLUMNS)


def _matched_score_metrics(common, scores, reason):
    """Return matched metrics for one score vector; undefined values stay NaN.

    Predictors define the common genes, so only a reference can miss scores there; it is
    then reported as unavailable instead of being scored on fewer genes.
    """
    metrics = {metric_name: float("nan") for metric_name in MATCHED_SUMMARY_METRICS}
    missing_genes = int(scores.isna().sum())
    if missing_genes:
        return metrics, f"reference_missing_on_{missing_genes}_of_{len(scores)}_common_genes"
    frame = common[["expected_effect"]].assign(matched_score=scores.astype(float))
    if len(frame):
        metrics["pearson"], metrics["spearman"] = ev._continuous_metrics(frame, score_col="matched_score")
    if not reason:
        metrics["pr_auc"], metrics["aps"] = ev._compute_pr_metrics(common["is_positive"], frame["matched_score"])
        metrics["auroc"] = ev._compute_auroc(common["is_positive"], frame["matched_score"])
    return metrics, reason


def _matched_metric_rows(
    metrics_df, joined_frames, *, tool_ids, fdr_threshold, effect_threshold, reference_seed=None,
):
    """Score every included predictor (and reference) on one shared population per dataset.

    The population is the usable genes scored by all included predictors, restricted to
    datasets with a metric row for every included predictor. References are scored on the
    same genes when ``reference_seed`` is given.
    """
    if not tool_ids:
        return []
    score_cols = [f"{ev.SCORE_PREFIX}{tool_id}" for tool_id in tool_ids]
    tools_by_dataset = {
        str(dataset_id): set(group.astype(str))
        for dataset_id, group in metrics_df.groupby("dataset_id")["tool_id"]
    }
    rows = []
    for joined in joined_frames or []:
        if joined.empty or "dataset_id" not in joined.columns:
            continue
        dataset_id = str(joined["dataset_id"].iloc[0])
        if not set(tool_ids) <= tools_by_dataset.get(dataset_id, set()):
            continue
        if any(score_col not in joined.columns for score_col in score_cols):
            continue
        keep_cols = ["gene_id", "logFC", "FDR", *score_cols]
        for optional in ("dataset_id", "perturbation", *ev.FDR_AUXILIARY_COLUMNS):
            if optional in joined.columns:
                keep_cols.append(optional)
        work = ev._filter_usable_gt_rows(joined[keep_cols], fdr_threshold=fdr_threshold)
        work = ev._annotate_ground_truth(work)
        work["is_positive"] = ev._positive_mask(
            work,
            fdr_threshold=fdr_threshold,
            effect_threshold=effect_threshold,
        ).astype(int)
        common = work.dropna(subset=score_cols)
        positives = int(common["is_positive"].sum())
        negatives = int(len(common) - positives)
        if common.empty:
            reason = "No genes are scored by every matched predictor."
        elif positives == 0:
            reason = "No positives remain in the matched population."
        elif negatives == 0:
            reason = "No negatives remain in the matched population."
        else:
            reason = ""
        entries = [
            (tool_id, "predictor", common[score_col])
            for tool_id, score_col in zip(tool_ids, score_cols)
        ]
        if reference_seed is not None:
            reference_scores = ev._reference_scores(joined, dataset_id=dataset_id, seed=reference_seed)
            for tool_id in sorted(ev.REFERENCE_TOOL_IDS):
                score_col = f"{ev.SCORE_PREFIX}{tool_id}"
                if score_col in reference_scores.columns:
                    scores = reference_scores.loc[common.index, score_col]
                else:
                    scores = pd.Series(float("nan"), index=common.index)
                entries.append((tool_id, "reference", scores))
        for tool_id, tool_type, scores in entries:
            metrics, metrics_na_reason = _matched_score_metrics(common, scores, reason)
            rows.append(
                {
                    "dataset_id": dataset_id,
                    "tool_id": tool_id,
                    "tool_type": tool_type,
                    "population": MATCHED_POPULATION,
                    "matched_tool_ids": ",".join(tool_ids),
                    "common_rows": int(len(common)),
                    "common_positives": positives,
                    **metrics,
                    "metrics_na_reason": metrics_na_reason,
                }
            )
    return rows


def _matched_cohort_summary(matched_df, *, cohort_tool_ids):
    """Summarize each metric over the datasets where every included predictor is finite.

    The cohort is metric-specific but identical for all tools, so every tool is averaged
    over the same comparisons. Undefined values are never filled; a reference without a
    finite value on every cohort dataset gets NaN statistics for that metric.
    """
    tool_types = matched_df.drop_duplicates("tool_id").set_index("tool_id")["tool_type"]
    summary = pd.DataFrame({"tool_id": tool_types.index, "tool_type": tool_types.to_numpy()})
    for metric_name in MATCHED_SUMMARY_METRICS:
        wide = matched_df.pivot(index="dataset_id", columns="tool_id", values=metric_name).astype(float)
        cohort = wide.loc[np.isfinite(wide[list(cohort_tool_ids)]).all(axis=1)]
        summary[f"{metric_name}_comparisons"] = int(len(cohort))
        for stat_name in MATCHED_COHORT_STATS:
            summary[f"{metric_name}_{stat_name}"] = [
                float(getattr(cohort[tool_id], stat_name)(skipna=False)) if len(cohort) else float("nan")
                for tool_id in summary["tool_id"]
            ]
    return summary


def _write_matched_summary(
    metrics_df,
    joined_frames,
    tables_dir,
    *,
    fdr_threshold,
    effect_threshold,
    configured_tool_ids=None,
    reference_seed=None,
    logger=None,
):
    roster = _matched_tool_roster(
        metrics_df,
        configured_tool_ids=configured_tool_ids,
        reference_baselines=reference_seed is not None,
    )
    predictors = roster[roster["tool_type"] == "predictor"]
    configured = predictors.loc[predictors["configured"], "tool_id"].tolist()
    included = predictors.loc[predictors["included"], "tool_id"].tolist()
    excluded = predictors.loc[~predictors["included"], "tool_id"].tolist()
    matched_df = pd.DataFrame(
        _matched_metric_rows(
            metrics_df,
            joined_frames,
            tool_ids=included,
            fdr_threshold=fdr_threshold,
            effect_threshold=effect_threshold,
            reference_seed=reference_seed,
        ),
        columns=MATCHED_METRIC_COLUMNS,
    )
    per_dataset_path = tables_dir / "matched_metrics_per_experiment.tsv"
    matched_df.to_csv(per_dataset_path, sep="\t", index=False)
    if matched_df.empty:
        summary = pd.DataFrame(columns=["tool_id", "tool_type"])
    else:
        summary = _matched_cohort_summary(matched_df, cohort_tool_ids=included)
    summary.insert(2, "is_reference", summary["tool_type"] == "reference")
    summary.insert(3, "population", MATCHED_POPULATION)
    summary.insert(4, "configured_tool_ids", ",".join(configured))
    summary.insert(5, "matched_tool_ids", ",".join(included))
    summary.insert(6, "excluded_tool_ids", ",".join(excluded))
    summary.insert(7, "matched_dataset_count", int(matched_df["dataset_id"].nunique()))
    summary_path = tables_dir / "cross_dataset_predictor_summary_matched.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)
    roster_path = tables_dir / "matched_tool_roster.tsv"
    roster.to_csv(roster_path, sep="\t", index=False)
    ev._emit_log(logger, f"  Wrote matched cross-dataset summary table: {summary_path}")
    if excluded:
        ev._emit_log(
            logger,
            "  Matched summary excludes configured predictor(s): "
            + ", ".join(f"{tool_id} ({NO_METRIC_ROWS_REASON})" for tool_id in excluded),
        )
    return {
        "cross_dataset_predictor_summary_matched": str(summary_path),
        "matched_metrics_per_experiment": str(per_dataset_path),
        "matched_tool_roster": str(roster_path),
    }


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
    tool_ids=None,
    reference_baselines=False,
    reference_seed=0,
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
    # The summary above averages each predictor over its own scored genes and evaluated
    # datasets; the matched summary compares predictors on one shared population.
    matched_tables = {}
    if joined_frames:
        matched_tables = _write_matched_summary(
            metrics_df,
            joined_frames,
            tables_dir,
            fdr_threshold=fdr_threshold,
            effect_threshold=effect_threshold,
            configured_tool_ids=tool_ids,
            reference_seed=reference_seed if reference_baselines else None,
            logger=logger,
        )

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
            **matched_tables,
        },
        "plots": {
            **distribution_paths,
            **rank_distribution_paths,
        },
    }
