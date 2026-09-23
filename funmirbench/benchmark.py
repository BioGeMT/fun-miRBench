"""Run the fun-miRBenvh benchmark from a single YAML config."""

from __future__ import annotations

import argparse
import json
import logging
import pathlib
import shutil

import yaml

import funmirbench.evaluate as evaluate_module
from funmirbench.benchmark_config import (
    build_run_dir_name,
    filter_df,
    load_experiments,
    load_predictions,
    validate_predictor_output_files,
)
from funmirbench.benchmark_reports import (
    _init_run_layout,
    write_run_readme,
)
from funmirbench.common_predictions import (
    write_combined_common_prediction_summary,
    write_common_prediction_summary,
)
from funmirbench.comparison_plots import write_common_comparison_plots
from funmirbench.cross_dataset import write_cross_dataset_summaries, write_metric_tables
from funmirbench.dataset_reports import write_predictor_reports
from funmirbench.evaluate import (
    REPORT_PAGE_SIZE,
    evaluate_joined_dataframe,
)
from funmirbench.join import build_joined, load_predictor_score_cache
from funmirbench.logger import parse_log_level, setup_logging
from funmirbench.predictor_combinations import write_predictor_combination_outputs
from funmirbench.protein_coding import (
    DEFAULT_CACHE_REL_PATH,
    DEFAULT_GTF_REL_PATH,
    ENSEMBL_RELEASE,
    load_protein_coding_gene_ids,
)
from funmirbench.run_report import write_run_pdf_report
from funmirbench.validate_experiments import (
    format_validation_failure,
    log_validation_summary,
    validate_experiments,
)


logger = logging.getLogger(__name__)


def parse_args():
    parser = argparse.ArgumentParser(description="Run the fun-miRBenvh benchmark.")
    parser.add_argument("--config", type=pathlib.Path, required=True)
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
    )
    return parser.parse_args()


def _optional_float(value, default=None):
    if value is None:
        return default
    return float(value)


def _optional_bool(value, default=False):
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in {"1", "true", "yes", "y", "on"}:
            return True
        if normalized in {"0", "false", "no", "n", "off"}:
            return False
    raise ValueError(f"Expected a boolean value, got {value!r}")


def _positive_int(value, *, name):
    parsed = int(value)
    if parsed < 1:
        raise ValueError(f"{name} must be >= 1, got {value!r}")
    return parsed


def clear_dataset_outputs(dataset_id, plots_dir, reports_dir):
    dataset_plots_dir = plots_dir / dataset_id
    if dataset_plots_dir.exists():
        shutil.rmtree(dataset_plots_dir)

    for stale_report in reports_dir.glob(f"{dataset_id}__*"):
        if stale_report.is_file():
            stale_report.unlink()


def _finalize_run_bundle(
    out_dir,
    *,
    out_root,
    config_path,
    config_snapshot_path,
    dataset_outputs,
    tool_ids,
    metric_rows,
    joined_frames,
    common_prediction_summaries,
    tool_labels,
    fdr_threshold,
    effect_threshold,
    predictor_top_fraction,
    protein_coding_filter,
    skipped_datasets,
):
    layout = _init_run_layout(out_dir)
    if metric_rows:
        metric_tables = write_metric_tables(
            metric_rows,
            layout["per_experiment_tables_dir"],
            logger=logger.info,
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
            logger=logger.info,
        )
        common_summary_path = write_combined_common_prediction_summary(
            common_prediction_summaries,
            layout["combined_tables_dir"],
        )
        combined_outputs.setdefault("tables", {})["common_prediction_summary"] = str(common_summary_path)
        combination_outputs = write_predictor_combination_outputs(
            joined_frames,
            layout["combined_tables_dir"],
            layout["combined_plots_dir"],
            tool_ids=tool_ids,
            fdr_threshold=fdr_threshold,
            effect_threshold=effect_threshold,
            predictor_top_fraction=predictor_top_fraction,
            logger=logger.info,
        )
        combined_outputs.setdefault("tables", {}).update(combination_outputs.get("tables", {}))
        combined_outputs.setdefault("plots", {}).update(combination_outputs.get("plots", {}))
    else:
        metric_tables = {}
        combined_outputs = {"tables": {}, "plots": {}}
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
        protein_coding_filter=protein_coding_filter,
        skipped_datasets=skipped_datasets,
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
        skipped_datasets=skipped_datasets,
    )
    summary = {
        "config": str(config_path),
        "config_snapshot": str(config_snapshot_path),
        "out_root": str(out_root),
        "out_dir": str(out_dir),
        "run_dir_name": out_dir.name,
        "dataset_ids": [item["dataset_id"] for item in dataset_outputs],
        "tool_ids": tool_ids,
        "readme": str(readme_path),
        "report_pdf": str(report_path),
        "metric_tables": metric_tables,
        "cross_dataset_outputs": combined_outputs,
        "protein_coding_filter": protein_coding_filter,
        "datasets": dataset_outputs,
        "skipped_datasets": skipped_datasets,
    }
    summary_path = out_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    logger.info(f"Wrote summary: {summary_path}")
    return {
        "metric_tables": metric_tables,
        "combined_outputs": combined_outputs,
        "readme_path": readme_path,
        "report_path": report_path,
        "summary": summary,
        "summary_path": summary_path,
    }


def run_benchmark(config_path):
    config_path = config_path.expanduser().resolve()
    logger.info(f"Config: {config_path}")
    logger.info("Loading benchmark config...")
    with config_path.open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    root = config_path.parent
    experiments_tsv = root / config["experiments_tsv"]
    experiment_filters = config.get("experiments")
    eval_cfg = config.get("evaluation", {})
    raw_fdr_threshold = eval_cfg.get("fdr_threshold", 0.05)
    fdr_threshold = None if raw_fdr_threshold is None else float(raw_fdr_threshold)
    effect_threshold = float(eval_cfg.get("effect_threshold", 1.0))
    predictor_top_fraction = float(eval_cfg.get("predictor_top_fraction", 0.10))
    write_top_prediction_cdfs = _optional_bool(eval_cfg.get("write_top_prediction_cdfs"), True)
    report_min_common_coverage = float(
        eval_cfg.get("report_min_common_coverage", eval_cfg.get("publication_min_common_coverage", 0.10))
    )
    protein_coding_only = _optional_bool(eval_cfg.get("protein_coding_only"), True)
    predictor_load_workers = _positive_int(
        eval_cfg.get("predictor_load_workers", 1),
        name="evaluation.predictor_load_workers",
    )

    logger.info("Loading predictors...")
    predictions = load_predictions(
        root / config["predictions_tsv"],
        config.get("predictors"),
    )
    if not predictions:
        raise ValueError("Predictor selection resolved to no predictors.")

    validate_predictor_output_files(predictions, root)

    logger.info("Validating selected experiments...")
    validation_summary = validate_experiments(
        experiments_tsv,
        root=root,
        filters=experiment_filters,
        fdr_threshold=fdr_threshold,
        effect_threshold=effect_threshold,
    )
    log_validation_summary(validation_summary)
    if not validation_summary.ok:
        raise ValueError(format_validation_failure(validation_summary))

    logger.info("Loading experiments...")
    experiments = load_experiments(
        experiments_tsv,
        root,
        experiment_filters,
    )
    if validation_summary.excluded_dataset_ids:
        experiments = [
            experiment
            for experiment in experiments
            if experiment.id not in validation_summary.excluded_dataset_ids
        ]
    if not experiments:
        raise ValueError("Experiment selection resolved to no datasets.")

    protein_coding_gene_ids = None
    protein_coding_filter = {
        "enabled": protein_coding_only,
        "ensembl_release": ENSEMBL_RELEASE,
        "gtf_path": str(eval_cfg.get("protein_coding_gtf") or DEFAULT_GTF_REL_PATH),
        "gene_cache": str(eval_cfg.get("protein_coding_gene_cache") or DEFAULT_CACHE_REL_PATH),
        "gene_count": None,
    }
    if protein_coding_only:
        logger.info("Protein-coding-only evaluation filter is enabled.")
        protein_coding_gene_ids = load_protein_coding_gene_ids(
            root=root,
            gtf_path=eval_cfg.get("protein_coding_gtf"),
            cache_path=eval_cfg.get("protein_coding_gene_cache"),
        )
        protein_coding_filter["gene_count"] = len(protein_coding_gene_ids)
        logger.info("Protein-coding gene set contains %d genes.", len(protein_coding_gene_ids))

    evaluate_module.FIGURE_DPI = int(eval_cfg.get("figure_dpi", eval_cfg.get("publication_figure_dpi", 450)))
    out_root = (root / config.get("out_dir", "results")).resolve()
    out_root.mkdir(parents=True, exist_ok=True)
    run_dir_name = build_run_dir_name()
    out_dir = out_root / run_dir_name
    if out_dir.exists():
        raise FileExistsError(
            f"Run output directory already exists: {out_dir}. "
            "Start the benchmark again to create a new timestamped run directory."
        )
    logger.info(f"Results root: {out_root}")
    logger.info(f"Run output dir: {out_dir}")
    main_layout = _init_run_layout(out_dir)
    config_snapshot_path = out_dir / "benchmark_config.yaml"
    shutil.copy2(config_path, config_snapshot_path)
    logger.info(f"Wrote benchmark config snapshot: {config_snapshot_path}")

    tool_ids = list(predictions)
    tool_labels = {
        str(tool_id): str(meta.get("official_name") or tool_id)
        for tool_id, meta in predictions.items()
    }
    metric_rows = []
    dataset_outputs = []
    joined_frames = []
    common_prediction_summaries = []
    skipped_datasets = []
    logger.info(f"Experiments: {len(experiments)}")
    logger.info(f"Predictors:  {tool_ids}")
    logger.info("Loading predictor score files once for this run with %d worker(s)...", predictor_load_workers)
    predictor_cache = load_predictor_score_cache(
        tool_ids,
        predictions,
        root,
        max_workers=predictor_load_workers,
        logger=logger.info,
    )

    for meta in experiments:
        logger.info(f"Dataset: {meta.id} | {meta.miRNA} | {meta.cell_line}")
        dataset_dir = main_layout["datasets_dir"] / meta.id
        if dataset_dir.exists():
            shutil.rmtree(dataset_dir)
        (dataset_dir / "plots").mkdir(parents=True, exist_ok=True)
        (dataset_dir / "reports").mkdir(parents=True, exist_ok=True)
        logger.info(f"  Joining predictions for {meta.id}...")
        joined, predictor_output_paths = build_joined(
            meta,
            tool_ids,
            predictions,
            root,
            protein_coding_gene_ids=protein_coding_gene_ids,
            predictor_cache=predictor_cache,
            logger=logger.info,
        )
        joined_path = dataset_dir / "joined.tsv"
        joined.to_csv(joined_path, sep="\t", index=False)
        logger.info(f"  Wrote joined table: {joined_path}")

        score_columns = [
            f"score_{tool_id}"
            for tool_id in tool_ids
            if f"score_{tool_id}" in joined.columns
        ]
        if not score_columns or not joined[score_columns].notna().any().any():
            reason = (
                f"No selected predictor has predictions for miRNA {meta.miRNA}."
            )
            logger.info("  Skipping %s: %s", meta.id, reason)
            skipped_datasets.append(
                {
                    "dataset_id": meta.id,
                    "mirna": meta.miRNA,
                    "cell_line": meta.cell_line,
                    "perturbation": meta.perturbation,
                    "reason": reason,
                    "joined_tsv": str(joined_path),
                }
            )
            continue

        logger.info(f"  Evaluating metrics and plots for {meta.id}...")
        evaluation = evaluate_joined_dataframe(
            joined,
            plots_dir=dataset_dir / "plots",
            reports_dir=dataset_dir / "reports",
            fdr_threshold=fdr_threshold,
            effect_threshold=effect_threshold,
            predictor_top_fraction=predictor_top_fraction,
            dataset_id=meta.id,
            mirna=meta.miRNA,
            cell_line=meta.cell_line,
            perturbation=meta.perturbation,
            geo_accession=meta.geo_accession,
            de_table_path=str(meta.full_path),
            joined_tsv=joined_path,
            predictor_output_paths=predictor_output_paths,
            tool_labels=tool_labels,
            write_top_prediction_cdfs=write_top_prediction_cdfs,
            logger=logger.info,
        )
        write_common_comparison_plots(
            joined,
            evaluation=evaluation,
            dataset_metric_rows=evaluation["metric_rows"],
            plots_dir=dataset_dir / "plots",
            dataset_id=meta.id,
            fdr_threshold=fdr_threshold,
            effect_threshold=effect_threshold,
            perturbation=meta.perturbation,
            min_common_coverage=report_min_common_coverage,
            logger=logger.info,
        )
        common_summary_path, common_prediction_summary = write_common_prediction_summary(
            joined,
            dataset_dir / "reports",
            dataset_id=meta.id,
            tool_ids=tool_ids,
            report_min_common_coverage=report_min_common_coverage,
        )
        common_prediction_summaries.append(common_prediction_summary)
        write_predictor_reports(
            reports_dir=dataset_dir / "reports",
            plots_dir=dataset_dir / "plots",
            dataset_id=meta.id,
            mirna=meta.miRNA,
            cell_line=meta.cell_line,
            perturbation=meta.perturbation,
            geo_accession=meta.geo_accession,
            de_table_path=str(meta.full_path),
            predictor_output_paths=predictor_output_paths,
            metric_rows=evaluation["metric_rows"],
            skipped_tool_rows=evaluation.get("skipped_tool_rows", []),
            tool_labels=tool_labels,
            fdr_threshold=fdr_threshold,
            effect_threshold=effect_threshold,
            common_prediction_summary=common_prediction_summary,
        )
        joined_frames.append(joined.copy())
        metric_rows.extend(evaluation["metric_rows"])
        dataset_outputs.append(
            {
                "dataset_id": meta.id,
                "mirna": meta.miRNA,
                "cell_line": meta.cell_line,
                "perturbation": meta.perturbation,
                "geo_accession": meta.geo_accession,
                "de_table_path": str(meta.full_path),
                "joined_tsv": str(joined_path),
                "dataset_dir": str(dataset_dir),
                "predictor_output_paths": predictor_output_paths,
                "plots": evaluation["plots"],
                "common_prediction_summary_tsv": str(common_summary_path),
            }
        )
        logger.info(f"  Finished {meta.id}")

    logger.info("Writing metric tables...")
    logger.info("Writing cross-dataset summaries...")
    _finalize_run_bundle(
        out_dir,
        out_root=out_root,
        config_path=config_path,
        config_snapshot_path=config_snapshot_path,
        dataset_outputs=dataset_outputs,
        tool_ids=tool_ids,
        metric_rows=metric_rows,
        joined_frames=joined_frames,
        common_prediction_summaries=common_prediction_summaries,
        tool_labels=tool_labels,
        fdr_threshold=fdr_threshold,
        effect_threshold=effect_threshold,
        predictor_top_fraction=predictor_top_fraction,
        protein_coding_filter=protein_coding_filter,
        skipped_datasets=skipped_datasets,
    )
    return out_dir


def main():
    args = parse_args()
    setup_logging(parse_log_level(args.log_level))
    out_dir = run_benchmark(args.config)
    logger.info(f"Done. Results in {out_dir}")


if __name__ == "__main__":
    main()
