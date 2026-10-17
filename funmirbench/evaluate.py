"""Evaluate joined GT/prediction tables: metrics, plots, reports."""

from funmirbench.evaluate_common import *
from funmirbench.evaluate_plots import *
from funmirbench.heatmap_plots import *
from funmirbench.evaluate_reports import *
from funmirbench.de_table import FDR_DERIVED_COLUMNS, add_fdr_derived_columns
from funmirbench.validate_experiments import (
    DEFAULT_ZERO_POSITIVE_POLICY,
    resolve_zero_positive_policy,
)


def _continuous_metric_row(identity, coverage_info, *, pearson, spearman):
    reason = coverage_info["class_metrics_na_reason"]
    return {
        **identity,
        "rows_total": coverage_info["rows_total"],
        "rows_scored": coverage_info["rows_scored"],
        "coverage": coverage_info["coverage"],
        "positives_total": coverage_info["positives_total"],
        "positives_scored": coverage_info["positives_scored"],
        "negatives_scored": coverage_info["rows_scored"] - coverage_info["positives_scored"],
        "rows_excluded_fdr_na": coverage_info["rows_excluded_fdr_na"],
        "spearman": spearman,
        "pearson": pearson,
        "class_metrics_status": "not_evaluable" if reason else "evaluated",
        "class_metrics_na_reason": reason,
    }


REFERENCE_ALL_GENES_POPULATION = "all_usable_genes"


def _evaluate_reference_baselines(
    joined, *, identity, fdr_threshold, effect_threshold, perturbation, seed,
    matched_gene_ids=None, logger=None,
):
    """Score reference baselines with the predictors' GT labels on two kinds of population.

    ``all_usable_genes`` rows use every usable gene with a reference score and are
    descriptive only: sparse predictors are scored on fewer genes. ``matched_to_<tool_id>``
    rows use exactly the usable genes that predictor scored, so predictor-vs-reference
    differences can be read on identical genes. A reference with a missing score on any of
    those genes is reported as unavailable (NaN metrics plus counts), never on fewer genes.
    """
    if "gene_id" not in joined.columns:
        return []
    base_cols = [
        column
        for column in joined.columns
        if not column.startswith((SCORE_PREFIX, GLOBAL_RANK_PREFIX, LOCAL_RANK_PREFIX))
    ]
    reference_scores = _reference_scores(joined, dataset_id=identity["dataset_id"], seed=seed)
    gene_ids = joined["gene_id"].astype(str)
    populations = [(REFERENCE_ALL_GENES_POPULATION, "", reference_scores, None)]
    for tool_id, tool_gene_ids in (matched_gene_ids or {}).items():
        in_population = gene_ids.isin(set(tool_gene_ids))
        matched_scores = reference_scores.copy()
        matched_scores.loc[~in_population, :] = float("nan")
        populations.append((f"matched_to_{tool_id}", tool_id, matched_scores, in_population))
    rows = []
    for population, matched_tool_id, population_scores, in_population in populations:
        rows.extend(
            _reference_population_rows(
                pd.concat([joined[base_cols], population_scores], axis=1),
                score_cols=list(population_scores.columns),
                identity=identity,
                population=population,
                matched_tool_id=matched_tool_id,
                in_population=in_population,
                fdr_threshold=fdr_threshold,
                effect_threshold=effect_threshold,
                perturbation=perturbation,
                seed=seed,
                logger=logger if not matched_tool_id else None,
            )
        )
    return rows


def _unavailable_reference_row(identity, *, tool_id, population, matched_tool_id, seed, matched_genes, missing_genes):
    return {
        **identity,
        "tool_id": tool_id,
        "is_reference": True,
        "population": population,
        "matched_tool_id": matched_tool_id,
        "matched_genes": matched_genes,
        "reference_missing_genes": missing_genes,
        "reference_description": REFERENCE_TOOL_DESCRIPTIONS[tool_id],
        "reference_seed": int(seed),
        **{
            column: float("nan")
            for column in (
                "rows_total", "rows_scored", "rows_missing_score", "coverage",
                "positives_total", "positives_scored", "positive_coverage",
                "aps", "pr_auc", "auroc", "spearman", "pearson",
            )
        },
        "class_metrics_na_reason": "",
        "reference_unavailable_reason": (
            f"reference_missing_on_{missing_genes}_of_{matched_genes}_predictor_genes"
        ),
    }


def _reference_population_rows(
    frame, *, score_cols, identity, population, matched_tool_id, in_population,
    fdr_threshold, effect_threshold, perturbation, seed, logger=None,
):
    """Score each reference on one population; ``in_population`` marks matched genes."""
    rows = []
    for score_col in score_cols:
        tool_id = _tool_id_from_score_col(score_col)
        matched_genes = missing_genes = float("nan")
        if in_population is not None:
            matched_genes = int(in_population.sum())
            missing_genes = int(frame.loc[in_population, score_col].isna().sum())
            if missing_genes:
                rows.append(
                    _unavailable_reference_row(
                        identity,
                        tool_id=tool_id,
                        population=population,
                        matched_tool_id=matched_tool_id,
                        seed=seed,
                        matched_genes=matched_genes,
                        missing_genes=missing_genes,
                    )
                )
                continue
        try:
            scored, coverage_info = _prepare_scored_frame(
                frame, score_col=score_col,
                fdr_threshold=fdr_threshold, effect_threshold=effect_threshold,
                perturbation=perturbation,
                require_both_classes=False,
            )
        except ValueError as exc:
            _emit_log(logger, f"    Reference: {tool_id} [{population}] | skipped: {exc}")
            continue
        pearson, spearman = _continuous_metrics(scored, score_col=score_col)
        aps = pr_auc = auroc = float("nan")
        if not coverage_info["class_metrics_na_reason"]:
            pr_auc, aps = _compute_pr_metrics(scored["is_positive"], scored[score_col])
            auroc = _compute_auroc(scored["is_positive"], scored[score_col])
        rows.append({
            **identity,
            "tool_id": tool_id,
            "is_reference": True,
            "population": population,
            "matched_tool_id": matched_tool_id,
            "matched_genes": matched_genes,
            "reference_missing_genes": missing_genes,
            "reference_description": REFERENCE_TOOL_DESCRIPTIONS[tool_id],
            "reference_seed": int(seed),
            "rows_total": coverage_info["rows_total"],
            "rows_scored": coverage_info["rows_scored"],
            "rows_missing_score": coverage_info["rows_missing_score"],
            "coverage": coverage_info["coverage"],
            "positives_total": coverage_info["positives_total"],
            "positives_scored": coverage_info["positives_scored"],
            "positive_coverage": coverage_info["positive_coverage"],
            "aps": aps, "pr_auc": pr_auc, "auroc": auroc,
            "spearman": spearman, "pearson": pearson,
            "class_metrics_na_reason": coverage_info["class_metrics_na_reason"],
            "reference_unavailable_reason": "",
        })
        _emit_log(
            logger,
            f"    Reference: {tool_id} [{population}] | APS={aps:.3f} | AUROC={auroc:.3f} "
            f"| Spearman={spearman:.3f}",
        )
    return rows


def evaluate_joined_dataframe(
    joined, *, plots_dir, reports_dir,
    fdr_threshold, effect_threshold, predictor_top_fraction,
    dataset_id=None, mirna=None, cell_line=None,
    perturbation=None, geo_accession=None,
    de_table_path=None, joined_tsv=None,
    predictor_output_paths=None, dataset_label=None,
    tool_labels=None,
    write_top_prediction_cdfs=False,
    fdr_na_policy=None,
    zero_positive_policy=DEFAULT_ZERO_POSITIVE_POLICY,
    reference_baselines=True,
    reference_seed=0,
    logger=None,
):
    zero_positive_policy = resolve_zero_positive_policy(zero_positive_policy)
    plots_dir.mkdir(parents=True, exist_ok=True)
    reports_dir.mkdir(parents=True, exist_ok=True)
    _set_tool_labels(tool_labels)

    score_cols = sorted(c for c in joined.columns if c.startswith(SCORE_PREFIX))
    if not score_cols:
        raise ValueError("No score_<tool_id> columns found in joined dataframe.")
    if fdr_na_policy is not None and "FDR" in joined.columns:
        # Re-derive benchmark_FDR in place so every downstream consumer of `joined`
        # applies the same policy; the DE loader already does this for benchmark runs.
        derived = add_fdr_derived_columns(joined, fdr_na_policy=fdr_na_policy)
        for column in FDR_DERIVED_COLUMNS:
            joined[column] = derived[column]

    dataset_id = dataset_id or (
        str(joined["dataset_id"].iloc[0]) if "dataset_id" in joined.columns else "NA"
    )
    mirna = mirna or (
        str(joined["mirna"].iloc[0]) if "mirna" in joined.columns else None
    )
    _set_dataset_label(dataset_id, dataset_label)
    dataset_plots_dir = plots_dir
    predictor_plots_dir = dataset_plots_dir / "predictors"
    comparison_plots_dir = dataset_plots_dir / "comparisons"
    heatmap_plots_dir = dataset_plots_dir / "heatmaps"
    for path in (dataset_plots_dir, predictor_plots_dir, comparison_plots_dir, heatmap_plots_dir):
        path.mkdir(parents=True, exist_ok=True)
    tool_ids = [_tool_id_from_score_col(sc) for sc in score_cols]
    publication_tool_ids = _publication_tool_ids(tool_ids)
    _set_tool_colors(publication_tool_ids or tool_ids)
    global_rank_cols = []
    local_rank_cols = []
    for score_col, tool_id in zip(score_cols, tool_ids):
        global_rank_col = _rank_col_for_tool(tool_id)
        local_rank_col = _rank_col_for_tool(tool_id, prefix=LOCAL_RANK_PREFIX)
        if local_rank_col not in joined.columns:
            joined[local_rank_col] = _rank_scale_scores(joined[score_col])
        if global_rank_col not in joined.columns:
            joined[global_rank_col] = joined[local_rank_col]
        global_rank_cols.append(global_rank_col)
        local_rank_cols.append(local_rank_col)

    metric_rows = []
    skipped_tool_rows = []
    continuous_metric_rows = []
    # Usable scored genes of each predictor with a metric row, for matched references.
    matched_gene_ids = {}
    dataset_plots = {}
    predictor_correlation_tsv = None
    comparisons = []
    coverage_by_tool = {}
    scored_frames_by_tool = {}
    evaluated_score_cols = []
    evaluated_tool_ids = []
    evaluated_local_rank_cols = []

    _emit_log(logger, f"    Evaluation start: {dataset_id} | tools={tool_ids}")

    for score_col, tool_id, local_rank_col in zip(score_cols, tool_ids, local_rank_cols):
        _emit_log(logger, f"    Tool: {tool_id} | preparing scored pairs")
        tool_plots_dir = predictor_plots_dir / tool_id
        tool_plots_dir.mkdir(parents=True, exist_ok=True)
        identity = {
            "dataset_id": dataset_id, "mirna": mirna, "cell_line": cell_line,
            "perturbation": perturbation, "geo_accession": geo_accession,
            "tool_id": tool_id,
        }
        try:
            scored, coverage_info = _prepare_scored_frame(
                joined, score_col=score_col,
                fdr_threshold=fdr_threshold, effect_threshold=effect_threshold,
                perturbation=perturbation,
                require_both_classes=False,
            )
        except ValueError as exc:
            message = str(exc)
            valid_rows = _usable_gt_row_mask(joined, fdr_threshold=fdr_threshold)
            scored_rows = int(joined.loc[valid_rows, score_col].notna().sum())
            total_valid_rows = int(valid_rows.sum())
            sparse_scoring = scored_rows < total_valid_rows
            if not (message.startswith("No scored rows remain") and sparse_scoring):
                raise
            _emit_log(logger, f"    Tool: {tool_id} | skipped: {exc}")
            positives_total = 0
            positives_scored = 0
            try:
                skip_cols = ["logFC", "FDR", score_col]
                for optional in FDR_AUXILIARY_COLUMNS:
                    if optional in joined.columns:
                        skip_cols.append(optional)
                skip_frame = joined.loc[valid_rows, skip_cols].copy()
                skip_frame = _annotate_ground_truth(skip_frame, perturbation=perturbation)
                positive_mask = _positive_mask(
                    skip_frame,
                    fdr_threshold=fdr_threshold,
                    effect_threshold=effect_threshold,
                )
                positives_total = int(positive_mask.sum())
                positives_scored = int(positive_mask.loc[skip_frame[score_col].notna()].sum())
            except (KeyError, ValueError, TypeError):
                pass
            skipped_tool_rows.append(
                {
                    "dataset_id": dataset_id,
                    "mirna": mirna,
                    "cell_line": cell_line,
                    "perturbation": perturbation,
                    "geo_accession": geo_accession,
                    "tool_id": tool_id,
                    "skip_reason": message,
                    "rows_total": total_valid_rows,
                    "rows_scored": scored_rows,
                    "rows_missing_score": max(total_valid_rows - scored_rows, 0),
                    "positives_total": positives_total,
                    "positives_scored": positives_scored,
                }
            )
            continue
        class_metrics_na_reason = coverage_info["class_metrics_na_reason"]
        if class_metrics_na_reason:
            # Continuous metrics need no GT classes, so they are always recorded.
            pearson, spearman = _continuous_metrics(scored, score_col=score_col)
            continuous_metric_rows.append(
                _continuous_metric_row(identity, coverage_info, pearson=pearson, spearman=spearman)
            )
            _emit_log(
                logger,
                f"    Tool: {tool_id} | class metrics NA: {class_metrics_na_reason} "
                f"| Spearman={spearman:.3f}",
            )
            if zero_positive_policy == "continuous_only":
                metric_rows.append({
                    **identity,
                    "rows_total": coverage_info["rows_total"],
                    "rows_scored": coverage_info["rows_scored"],
                    "rows_missing_score": coverage_info["rows_missing_score"],
                    "coverage": coverage_info["coverage"],
                    "positives_total": coverage_info["positives_total"],
                    "positives_scored": coverage_info["positives_scored"],
                    "positive_coverage": coverage_info["positive_coverage"],
                    **{metric_name: float("nan") for metric_name in CLASS_METRIC_NAMES},
                    "spearman": spearman, "pearson": pearson,
                    "rows_excluded_fdr_na": coverage_info["rows_excluded_fdr_na"],
                    "rows_excluded_fdr_na_scored": coverage_info["rows_excluded_fdr_na_scored"],
                    "class_metrics_na_reason": class_metrics_na_reason,
                })
                matched_gene_ids[tool_id] = scored["gene_id"].astype(str).tolist()
            else:
                skipped_tool_rows.append({
                    **identity,
                    "skip_reason": class_metrics_na_reason,
                    "rows_total": coverage_info["rows_total"],
                    "rows_scored": coverage_info["rows_scored"],
                    "rows_missing_score": coverage_info["rows_missing_score"],
                    "positives_total": coverage_info["positives_total"],
                    "positives_scored": coverage_info["positives_scored"],
                    "spearman": spearman,
                    "pearson": pearson,
                })
            continue
        scatter_png = tool_plots_dir / "score_vs_expected_effect.png"
        gsea_png = tool_plots_dir / "gsea_enrichment.png"
        pr_curve_png = tool_plots_dir / "precision_recall_curve.png"
        roc_curve_png = tool_plots_dir / "roc_curve.png"
        pearson, spearman = _plot_scatter_with_correlation(
            scored,
            score_col=score_col,
            dataset_id=dataset_id,
            tool_id=tool_id,
            positives_total=coverage_info["positives_total"],
            out_path=scatter_png,
        )
        enrichment_score = _plot_gsea_enrichment(
            scored,
            score_col=score_col,
            dataset_id=dataset_id,
            tool_id=tool_id,
            positives_total=coverage_info["positives_total"],
            out_path=gsea_png,
        )
        pr_auc, aps = _compute_pr_metrics(scored["is_positive"], scored[score_col])
        auroc = _compute_auroc(scored["is_positive"], scored[score_col])
        _plot_single_predictor_pr_curve(
            {
                "tool_id": tool_id,
                "y_true": scored["is_positive"],
                "y_score": scored[score_col],
                "coverage": coverage_info["coverage"],
                "positives_total": coverage_info["positives_total"],
            },
            dataset_id=dataset_id,
            out_path=pr_curve_png,
        )
        _plot_single_predictor_roc_curve(
            {
                "tool_id": tool_id,
                "y_true": scored["is_positive"],
                "y_score": scored[score_col],
                "coverage": coverage_info["coverage"],
                "positives_total": coverage_info["positives_total"],
            },
            dataset_id=dataset_id,
            out_path=roc_curve_png,
        )

        metrics = {
            "rows_total": float(coverage_info["rows_total"]),
            "rows_used": float(len(scored)),
            "positives": float(scored["is_positive"].sum()),
            "negatives": float(len(scored) - int(scored["is_positive"].sum())),
            "rows_missing_score": float(coverage_info["rows_missing_score"]),
            "coverage": float(coverage_info["coverage"]),
            "positives_total": float(coverage_info["positives_total"]),
            "positives_scored": float(coverage_info["positives_scored"]),
            "positive_coverage": float(coverage_info["positive_coverage"]),
            "pearson": pearson, "spearman": spearman,
            "aps": aps, "pr_auc": pr_auc, "auroc": auroc,
        }
        if _is_publication_tool(tool_id):
            report_md = reports_dir / f"{dataset_id}__{tool_id}_evaluation_report.md"
            report_pdf = reports_dir / f"{dataset_id}__{tool_id}_evaluation_report.pdf"
            _write_tool_report(
                dataset_id=dataset_id, mirna=mirna, cell_line=cell_line,
                perturbation=perturbation, geo_accession=geo_accession,
                de_table_path=de_table_path, joined_tsv=joined_tsv,
                tool_id=tool_id,
                predictor_output_path=(predictor_output_paths or {}).get(tool_id),
                metrics=metrics, markdown_path=report_md, pdf_path=report_pdf,
                coverage_info=coverage_info,
                scatter_png=scatter_png,
                pr_curve_png=pr_curve_png,
                roc_curve_png=roc_curve_png,
                gsea_png=gsea_png,
                fdr_threshold=fdr_threshold,
                effect_threshold=effect_threshold,
            )
        _emit_log(
            logger,
            (
                f"    Tool: {tool_id} | coverage={coverage_info['coverage']:.1%} "
                f"| positive_cov={coverage_info['positive_coverage']:.1%} "
                f"| rows={coverage_info['rows_scored']}/{coverage_info['rows_total']} "
                f"| APS={aps:.3f} | AUROC={auroc:.3f} | ES={enrichment_score:.3f}"
            ),
        )
        if _is_publication_tool(tool_id):
            _emit_log(logger, f"    Tool: {tool_id} | wrote scatter/report")
        else:
            _emit_log(logger, f"    Tool: {tool_id} | wrote diagnostic plots")

        metric_rows.append({
            **identity,
            "rows_total": coverage_info["rows_total"],
            "rows_scored": coverage_info["rows_scored"],
            "rows_missing_score": coverage_info["rows_missing_score"],
            "coverage": coverage_info["coverage"],
            "positives_total": coverage_info["positives_total"],
            "positives_scored": coverage_info["positives_scored"],
            "positive_coverage": coverage_info["positive_coverage"],
            "aps": aps, "spearman": spearman, "pearson": pearson,
            "auroc": auroc, "pr_auc": pr_auc, "gsea_es": enrichment_score,
            "rows_excluded_fdr_na": coverage_info["rows_excluded_fdr_na"],
            "rows_excluded_fdr_na_scored": coverage_info["rows_excluded_fdr_na_scored"],
            "class_metrics_na_reason": "",
        })
        continuous_metric_rows.append(
            _continuous_metric_row(identity, coverage_info, pearson=pearson, spearman=spearman)
        )
        matched_gene_ids[tool_id] = scored["gene_id"].astype(str).tolist()
        comparisons.append({
            "tool_id": tool_id,
            "y_true": scored["is_positive"],
            "y_score": scored[score_col],
            "coverage": coverage_info["coverage"],
        })
        coverage_by_tool[tool_id] = coverage_info["coverage"]
        scored_frames_by_tool[tool_id] = scored
        evaluated_score_cols.append(score_col)
        evaluated_tool_ids.append(tool_id)
        evaluated_local_rank_cols.append(local_rank_col)
        dataset_plots[f"{tool_id}_scatter"] = str(scatter_png)
        dataset_plots[f"{tool_id}_gsea_enrichment"] = str(gsea_png)
        dataset_plots[f"{tool_id}_pr_curve"] = str(pr_curve_png)
        dataset_plots[f"{tool_id}_roc_curve"] = str(roc_curve_png)

    reference_metric_rows = []
    if reference_baselines:
        reference_metric_rows = _evaluate_reference_baselines(
            joined,
            identity={
                "dataset_id": dataset_id, "mirna": mirna, "cell_line": cell_line,
                "perturbation": perturbation, "geo_accession": geo_accession,
            },
            fdr_threshold=fdr_threshold,
            effect_threshold=effect_threshold,
            perturbation=perturbation,
            seed=reference_seed,
            matched_gene_ids=matched_gene_ids,
            logger=logger,
        )

    if not evaluated_score_cols:
        if metric_rows:
            _emit_log(
                logger,
                f"    Dataset: {dataset_id} | class metrics NA for every predictor; "
                "recorded continuous metrics only",
            )
            return {
                "metric_rows": metric_rows,
                "skipped_tool_rows": skipped_tool_rows,
                "plots": dataset_plots,
                "predictor_correlation_tsv": None,
                "continuous_metric_rows": continuous_metric_rows,
                "reference_metric_rows": reference_metric_rows,
            }
        reason = "No selected predictor had evaluable scored rows after DE/GT filtering."
        _emit_log(logger, f"    Dataset: {dataset_id} | skipped: {reason}")
        return {
            "metric_rows": [],
            "skipped_tool_rows": skipped_tool_rows,
            "plots": dataset_plots,
            "predictor_correlation_tsv": None,
            "dataset_skip_reason": reason,
            "continuous_metric_rows": continuous_metric_rows,
            "reference_metric_rows": reference_metric_rows,
        }

    plot_entries = [
        (score_col, tool_id, rank_col, comparison)
        for score_col, tool_id, rank_col, comparison in zip(
            evaluated_score_cols,
            evaluated_tool_ids,
            evaluated_local_rank_cols,
            comparisons,
        )
        if _is_publication_tool(tool_id)
    ]
    if not plot_entries:
        plot_entries = list(
            zip(evaluated_score_cols, evaluated_tool_ids, evaluated_local_rank_cols, comparisons)
        )
    plot_score_cols = [item[0] for item in plot_entries]
    plot_tool_ids = [item[1] for item in plot_entries]
    plot_local_rank_cols = [item[2] for item in plot_entries]
    plot_comparisons = [item[3] for item in plot_entries]

    heatmap_png = heatmap_plots_dir / "algorithms_vs_genes.png"
    _plot_algorithms_vs_genes_heatmap(
        joined,
        score_cols=plot_score_cols,
        rank_cols=plot_local_rank_cols,
        tool_ids=plot_tool_ids,
        dataset_id=dataset_id, out_path=heatmap_png,
        fdr_threshold=fdr_threshold, effect_threshold=effect_threshold,
        perturbation=perturbation,
    )
    dataset_plots["algorithms_vs_genes_heatmap"] = str(heatmap_png)
    _emit_log(logger, f"    Dataset: {dataset_id} | wrote gene-level heatmap")

    top_positive_heatmap_png = heatmap_plots_dir / "top_10pct_positive_genes.png"
    wrote_top_positive_heatmap = _plot_top_positive_heatmap(
        joined,
        rank_cols=plot_local_rank_cols,
        tool_ids=plot_tool_ids,
        dataset_id=dataset_id,
        out_path=top_positive_heatmap_png,
        fdr_threshold=fdr_threshold,
        effect_threshold=effect_threshold,
        positive_fraction=0.10,
        perturbation=perturbation,
    )
    if wrote_top_positive_heatmap:
        dataset_plots["top_10pct_positive_heatmap"] = str(top_positive_heatmap_png)
        _emit_log(logger, f"    Dataset: {dataset_id} | wrote top-positive heatmap")

    if len(plot_score_cols) >= 2:
        comparison_pr_png = comparison_plots_dir / "precision_recall_common.png"
        comparison_pr_all_png = comparison_plots_dir / "precision_recall_all_scored.png"
        comparison_roc_png = comparison_plots_dir / "roc_common.png"
        comparison_roc_all_png = comparison_plots_dir / "roc_all_scored.png"
        comparison_gsea_png = comparison_plots_dir / "gsea_common.png"
        comparison_cdf_png = comparison_plots_dir / f"top_{TOP_PREDICTION_CDF_N}_effect_cdfs.png"
        try:
            common_pr = _prepare_common_scored_frame(
                joined,
                score_cols=plot_score_cols,
                fdr_threshold=fdr_threshold,
                effect_threshold=effect_threshold,
                perturbation=perturbation,
            )
            common_comparisons = [
                {
                    "tool_id": tool_id,
                    "y_true": common_pr["is_positive"],
                    "y_score": common_pr[score_col],
                    "coverage": coverage_by_tool.get(tool_id, float("nan")),
                }
                for score_col, tool_id in zip(plot_score_cols, plot_tool_ids)
            ]
            _plot_predictor_pr_curves(common_comparisons, dataset_id=dataset_id, out_path=comparison_pr_png)
            _plot_predictor_roc_curves(common_comparisons, dataset_id=dataset_id, out_path=comparison_roc_png)
            _plot_predictor_gsea_curves(common_comparisons, dataset_id=dataset_id, out_path=comparison_gsea_png)
            dataset_plots["predictor_pr_curves"] = str(comparison_pr_png)
            dataset_plots["predictor_roc_curves"] = str(comparison_roc_png)
            dataset_plots["predictor_gsea_curves"] = str(comparison_gsea_png)
        except ValueError as exc:
            _emit_log(logger, f"    Dataset: {dataset_id} | skipped common comparison plots: {exc}")

        try:
            own_scored_comparisons = []
            for score_col, tool_id, comparison in zip(plot_score_cols, plot_tool_ids, plot_comparisons):
                scored = scored_frames_by_tool[tool_id]
                own_scored_comparisons.append(
                    {
                        "tool_id": tool_id,
                        "y_true": scored["is_positive"],
                        "y_score": scored[score_col],
                        "coverage": comparison["coverage"],
                    }
                )
            _plot_predictor_pr_curves_own_scored(
                own_scored_comparisons,
                dataset_id=dataset_id,
                out_path=comparison_pr_all_png,
            )
            _plot_predictor_roc_curves_own_scored(
                own_scored_comparisons,
                dataset_id=dataset_id,
                out_path=comparison_roc_all_png,
            )
            dataset_plots["predictor_pr_curves_all_scored"] = str(comparison_pr_all_png)
            dataset_plots["predictor_roc_curves_all_scored"] = str(comparison_roc_all_png)
            if write_top_prediction_cdfs:
                _plot_top_prediction_effect_cdfs(
                    joined,
                    score_cols=plot_score_cols,
                    tool_ids=plot_tool_ids,
                    dataset_id=dataset_id,
                    out_path=comparison_cdf_png,
                    top_n=TOP_PREDICTION_CDF_N,
                    perturbation=perturbation,
                )
                dataset_plots["top_prediction_effect_cdfs"] = str(comparison_cdf_png)
            _emit_log(logger, f"    Dataset: {dataset_id} | wrote PR/ROC/GSEA comparison plots")
        except ValueError as exc:
            _emit_log(logger, f"    Dataset: {dataset_id} | skipped own-scored comparison plots: {exc}")

    if len(plot_local_rank_cols) >= 2:
        corr_tsv = heatmap_plots_dir / "predictor_rank_correlation.tsv"
        corr = _build_predictor_correlation_matrix(
            joined,
            rank_cols=plot_local_rank_cols,
            tool_ids=plot_tool_ids,
            top_fraction=predictor_top_fraction,
        )
        corr.to_csv(corr_tsv, sep="\t")
        predictor_correlation_tsv = str(corr_tsv)
        _emit_log(logger, f"    Dataset: {dataset_id} | wrote predictor correlation table")

    _emit_log(logger, f"    Evaluation complete: {dataset_id}")
    return {
        "metric_rows": metric_rows,
        "skipped_tool_rows": skipped_tool_rows,
        "plots": dataset_plots,
        "predictor_correlation_tsv": predictor_correlation_tsv,
        "continuous_metric_rows": continuous_metric_rows,
        "reference_metric_rows": reference_metric_rows,
    }
