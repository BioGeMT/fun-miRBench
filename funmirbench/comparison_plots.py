"""Common comparison plot helpers layered on top of evaluation outputs.

These helpers keep the main metric evaluation unchanged. They repair or augment
figure products when an ultra-sparse predictor would otherwise make common-scored
comparison plots disappear.
"""

from __future__ import annotations

import pathlib

from funmirbench.evaluate import (
    _is_publication_tool,
    SCORE_PREFIX,
    _plot_predictor_gsea_curves,
    _plot_predictor_pr_curves,
    _plot_predictor_roc_curves,
    _prepare_common_scored_frame,
)


DEFAULT_MIN_COMMON_COVERAGE = 0.10


def _score_col_for_tool(tool_id: str) -> str:
    return f"{SCORE_PREFIX}{tool_id}"


def _eligible_tools(metric_rows, *, min_coverage: float) -> list[str]:
    eligible = []
    for row in metric_rows:
        tool_id = str(row.get("tool_id", ""))
        if not tool_id or not _is_publication_tool(tool_id):
            continue
        if row.get("class_metrics_na_reason"):
            # Continuous-only rows (no positives or negatives) cannot enter PR/ROC plots.
            continue
        try:
            coverage = float(row.get("coverage"))
        except (TypeError, ValueError):
            continue
        if coverage >= float(min_coverage):
            eligible.append(tool_id)
    return sorted(dict.fromkeys(eligible))


def write_common_comparison_plots(
    joined,
    *,
    evaluation,
    dataset_metric_rows,
    plots_dir,
    dataset_id,
    fdr_threshold,
    effect_threshold,
    perturbation=None,
    min_common_coverage: float = DEFAULT_MIN_COMMON_COVERAGE,
    logger=None,
):
    """Write common-scored PR/ROC/GSEA plots for coverage-eligible predictors.

    The standard evaluator uses every selected predictor for common-set plots.
    That is statistically strict, but it can remove useful common comparison
    figures when a very sparse method is present. This step keeps all predictor
    metrics and all own-scored plots intact, but regenerates common comparison
    figures using predictors with at least ``min_common_coverage`` within this
    dataset.
    """
    eligible_tools = [
        tool_id
        for tool_id in _eligible_tools(dataset_metric_rows, min_coverage=min_common_coverage)
        if _score_col_for_tool(tool_id) in joined.columns
    ]
    if len(eligible_tools) < 2:
        if logger is not None:
            logger(
                f"  Common plots skipped for {dataset_id}: "
                f"only {len(eligible_tools)} predictor(s) reached {min_common_coverage:.0%} coverage."
            )
        return []

    score_cols = [_score_col_for_tool(tool_id) for tool_id in eligible_tools]
    try:
        common = _prepare_common_scored_frame(
            joined,
            score_cols=score_cols,
            fdr_threshold=fdr_threshold,
            effect_threshold=effect_threshold,
            perturbation=perturbation,
        )
    except ValueError as exc:
        if logger is not None:
            logger(f"  Common plots skipped for {dataset_id}: {exc}")
        return []

    comparisons = []
    for tool_id, score_col in zip(eligible_tools, score_cols):
        comparisons.append(
            {
                "tool_id": tool_id,
                "y_true": common["is_positive"].astype(int).to_numpy(),
                "y_score": common[score_col].astype(float).to_numpy(),
                "gene_id": common["gene_id"].astype(str).to_numpy() if "gene_id" in common.columns else None,
            }
        )

    comparisons_dir = pathlib.Path(plots_dir) / "comparisons"
    comparisons_dir.mkdir(parents=True, exist_ok=True)
    outputs = []

    pr_path = comparisons_dir / "precision_recall_common.png"
    _plot_predictor_pr_curves(comparisons, dataset_id=dataset_id, out_path=pr_path)
    evaluation["plots"]["predictor_pr_curves"] = str(pr_path)
    outputs.append(pr_path)

    roc_path = comparisons_dir / "roc_common.png"
    _plot_predictor_roc_curves(comparisons, dataset_id=dataset_id, out_path=roc_path)
    evaluation["plots"]["predictor_roc_curves"] = str(roc_path)
    outputs.append(roc_path)

    gsea_path = comparisons_dir / "gsea_common.png"
    _plot_predictor_gsea_curves(comparisons, dataset_id=dataset_id, out_path=gsea_path)
    evaluation["plots"]["predictor_gsea_curves"] = str(gsea_path)
    outputs.append(gsea_path)

    if logger is not None:
        logger(
            f"  Wrote common plots for {dataset_id} using "
            f"{len(eligible_tools)} predictors with >= {min_common_coverage:.0%} coverage: "
            + ", ".join(eligible_tools)
        )
    return outputs
