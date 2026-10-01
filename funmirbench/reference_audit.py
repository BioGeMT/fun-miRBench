"""Audit reference-label composition against expression and provenance signals."""

from __future__ import annotations

import argparse
import json
import logging
import pathlib
from dataclasses import dataclass
from typing import Iterable

import numpy as np
import pandas as pd
import yaml

from funmirbench.benchmark_config import filter_df
from funmirbench.de_table import find_gene_id_column, read_de_table
from funmirbench.experiment_store import sync_zenodo_experiments
from funmirbench.logger import parse_log_level, setup_logging
from funmirbench.protein_coding import (
    DEFAULT_CACHE_REL_PATH,
    DEFAULT_GTF_REL_PATH,
    ENSEMBL_RELEASE,
    load_protein_coding_gene_ids,
)

logger = logging.getLogger(__name__)

PERTURBATION_EFFECT_SIGN = {
    "overexpression": -1.0,
    "oe": -1.0,
    "knockout": 1.0,
    "ko": 1.0,
    "knockdown": 1.0,
    "kd": 1.0,
}

BASELINE_EXPRESSION_ALIASES = (
    "normalized_control_mean",
    "raw_control_mean",
    "control_mean",
    "controlMean",
    "mean_control",
    "control_baseMean",
    "baseMean_control",
    "baseMean",
)
TREATED_EXPRESSION_ALIASES = (
    "normalized_treated_mean",
    "raw_treated_mean",
    "treated_mean",
    "treatedMean",
    "mean_treated",
    "condition_mean",
    "treated_baseMean",
    "baseMean_treated",
)
UNCERTAINTY_COLUMN_ALIASES = (
    "lfcSE",
    "logFC_SE",
    "log2FoldChangeSE",
    "stat",
    "dispersion",
)


@dataclass(frozen=True)
class ReferenceAuditConfig:
    config_path: pathlib.Path
    out_dir: pathlib.Path
    fdr_threshold: float | None
    effect_threshold: float
    protein_coding_only: bool
    protein_coding_gtf: str | None
    protein_coding_gene_cache: str | None
    sync_experiments: bool
    expression_bins: int
    write_row_table: bool


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


def normalize_perturbation(value: object) -> str:
    text = str(value or "").strip().lower()
    if text not in PERTURBATION_EFFECT_SIGN:
        raise ValueError(
            f"Unsupported perturbation {value!r}; expected one of "
            "Overexpression, Knockout, or Knockdown."
        )
    return text


def perturbation_effect_sign(value: object) -> float:
    return PERTURBATION_EFFECT_SIGN[normalize_perturbation(value)]


def first_present_column(df: pd.DataFrame, candidates: Iterable[str]) -> str | None:
    columns_by_key = {str(column).lower(): str(column) for column in df.columns}
    for candidate in candidates:
        column = columns_by_key.get(candidate.lower())
        if column is not None:
            return column
    return None


def prepare_de_frame(
    df: pd.DataFrame,
    *,
    dataset_id: str,
    mirna: str,
    perturbation: str,
    fdr_threshold: float | None,
    effect_threshold: float,
    protein_coding_gene_ids: set[str] | None = None,
) -> tuple[pd.DataFrame, dict]:
    """Return an analysis-ready DE frame and a compact row accounting summary."""
    work = df.copy()
    gene_col = find_gene_id_column(work)
    if gene_col == "__index__":
        work["gene_id"] = work.index.astype(str)
    elif gene_col != "gene_id":
        work = work.rename(columns={gene_col: "gene_id"})
    work["gene_id"] = work["gene_id"].astype(str).str.replace(r"\.[0-9]+$", "", regex=True)

    rows_raw = len(work)
    if protein_coding_gene_ids is not None:
        work = work.loc[work["gene_id"].isin(protein_coding_gene_ids)].copy()
    rows_after_gene_universe = len(work)

    sign = perturbation_effect_sign(perturbation)
    work["dataset_id"] = dataset_id
    work["mirna"] = mirna
    work["perturbation"] = perturbation
    work["logFC"] = pd.to_numeric(work.get("logFC"), errors="coerce")
    if "FDR" in work.columns:
        work["FDR"] = pd.to_numeric(work["FDR"], errors="coerce")
    else:
        work["FDR"] = np.nan
    work["expected_effect"] = work["logFC"] * sign

    usable = work["logFC"].notna()
    if fdr_threshold is not None:
        usable = usable & work["FDR"].notna() & work["FDR"].between(0.0, 1.0)
    work["is_assessable"] = usable

    positive = work["is_assessable"] & (work["expected_effect"] > float(effect_threshold))
    if fdr_threshold is not None:
        positive = positive & (work["FDR"] < float(fdr_threshold))
    work["is_current_positive"] = positive

    baseline_col = first_present_column(work, BASELINE_EXPRESSION_ALIASES)
    treated_col = first_present_column(work, TREATED_EXPRESSION_ALIASES)
    if baseline_col:
        work["baseline_expression"] = pd.to_numeric(work[baseline_col], errors="coerce")
    else:
        work["baseline_expression"] = np.nan
    if treated_col:
        work["treated_expression"] = pd.to_numeric(work[treated_col], errors="coerce")
    else:
        work["treated_expression"] = np.nan

    summary = {
        "dataset_id": dataset_id,
        "mirna": mirna,
        "perturbation": perturbation,
        "rows_raw": rows_raw,
        "rows_after_gene_universe": rows_after_gene_universe,
        "rows_assessable": int(work["is_assessable"].sum()),
        "rows_missing_logFC": int(work["logFC"].isna().sum()),
        "rows_missing_FDR": int(work["FDR"].isna().sum()),
        "positives": int(work["is_current_positive"].sum()),
        "positive_rate": _safe_rate(work["is_current_positive"].sum(), work["is_assessable"].sum()),
        "baseline_expression_column": baseline_col or "",
        "treated_expression_column": treated_col or "",
        "available_uncertainty_columns": ",".join(
            column for column in UNCERTAINTY_COLUMN_ALIASES if column in work.columns
        ),
    }
    return work, summary


def _safe_rate(numerator, denominator) -> float:
    denominator = int(denominator)
    if denominator == 0:
        return float("nan")
    return float(numerator) / float(denominator)


def assign_expression_bins(values: pd.Series, *, bins: int) -> pd.Series:
    """Bin baseline expression with zeros separated from positive quantiles."""
    numeric = pd.to_numeric(values, errors="coerce")
    labels = pd.Series("missing", index=values.index, dtype="object")
    labels.loc[numeric.eq(0.0)] = "zero"

    positive = numeric[numeric.gt(0.0)]
    if positive.empty:
        return labels
    unique_count = positive.nunique(dropna=True)
    bin_count = max(1, min(int(bins), int(unique_count)))
    if bin_count == 1:
        labels.loc[positive.index] = "positive"
        return labels

    quantiles = pd.qcut(positive, q=bin_count, duplicates="drop")
    categories = list(quantiles.cat.categories)
    name_by_interval = {
        category: f"q{index + 1}_positive"
        for index, category in enumerate(categories)
    }
    labels.loc[positive.index] = quantiles.map(name_by_interval).astype(str)
    return labels


def summarize_expression_bins(frame: pd.DataFrame, *, bins: int) -> pd.DataFrame:
    if frame["baseline_expression"].notna().sum() == 0:
        return pd.DataFrame()
    work = frame.loc[frame["is_assessable"]].copy()
    work["baseline_bin"] = assign_expression_bins(work["baseline_expression"], bins=bins)
    grouped = []
    for (dataset_id, baseline_bin), group in work.groupby(["dataset_id", "baseline_bin"], dropna=False):
        grouped.append(_summarize_group(group, dataset_id=dataset_id, group_name=baseline_bin))
    return pd.DataFrame(grouped)


def summarize_zero_categories(frame: pd.DataFrame) -> pd.DataFrame:
    if frame["baseline_expression"].notna().sum() == 0:
        return pd.DataFrame()
    work = frame.loc[frame["is_assessable"]].copy()
    if work.empty:
        return pd.DataFrame()

    baseline_zero = work["baseline_expression"].eq(0.0)
    if work["treated_expression"].notna().sum() == 0:
        work["zero_category"] = np.where(baseline_zero, "baseline_zero", "baseline_positive_or_missing")
    else:
        treated_zero = work["treated_expression"].eq(0.0)
        work["zero_category"] = np.select(
            [
                baseline_zero & treated_zero,
                baseline_zero & ~treated_zero,
                ~baseline_zero & treated_zero,
            ],
            [
                "baseline_zero_treated_zero",
                "baseline_zero_treated_positive",
                "baseline_positive_treated_zero",
            ],
            default="baseline_positive_treated_positive_or_missing",
        )
    grouped = []
    for (dataset_id, category), group in work.groupby(["dataset_id", "zero_category"], dropna=False):
        grouped.append(_summarize_group(group, dataset_id=dataset_id, group_name=category))
    return pd.DataFrame(grouped).rename(columns={"group": "zero_category"})


def summarize_gene_response_propensity(frame: pd.DataFrame) -> pd.DataFrame:
    work = frame.loc[frame["is_assessable"]].copy()
    if work.empty:
        return pd.DataFrame()
    grouped = work.groupby("gene_id", dropna=False)
    out = grouped.agg(
        assessable_experiments=("dataset_id", "nunique"),
        assessable_rows=("dataset_id", "size"),
        positive_rows=("is_current_positive", "sum"),
        mean_expected_effect=("expected_effect", "mean"),
        median_expected_effect=("expected_effect", "median"),
        mean_logFC=("logFC", "mean"),
        median_logFC=("logFC", "median"),
    ).reset_index()
    out["positive_fraction"] = out["positive_rows"] / out["assessable_rows"]
    return out.sort_values(
        ["positive_fraction", "assessable_rows", "mean_expected_effect"],
        ascending=[False, False, False],
        kind="mergesort",
    )


def _summarize_group(group: pd.DataFrame, *, dataset_id: str, group_name: str) -> dict:
    logfc = pd.to_numeric(group["logFC"], errors="coerce")
    expected = pd.to_numeric(group["expected_effect"], errors="coerce")
    q75, q25 = np.nanpercentile(logfc.dropna(), [75, 25]) if logfc.notna().any() else (np.nan, np.nan)
    return {
        "dataset_id": dataset_id,
        "group": str(group_name),
        "rows": len(group),
        "positives": int(group["is_current_positive"].sum()),
        "positive_rate": _safe_rate(group["is_current_positive"].sum(), len(group)),
        "median_logFC": float(logfc.median()) if logfc.notna().any() else float("nan"),
        "iqr_logFC": float(q75 - q25) if logfc.notna().any() else float("nan"),
        "variance_logFC": float(logfc.var(ddof=1)) if logfc.notna().sum() > 1 else float("nan"),
        "median_expected_effect": float(expected.median()) if expected.notna().any() else float("nan"),
    }


def load_benchmark_config(config_path: pathlib.Path) -> dict:
    with config_path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def selected_experiment_paths(tsv_path: pathlib.Path, filters: dict | None) -> list[str]:
    df = pd.read_csv(tsv_path, sep="\t")
    if filters:
        df = filter_df(df, filters)
    return [str(value) for value in df["de_table_path"].dropna().tolist()]


def build_audit_config(args: argparse.Namespace) -> ReferenceAuditConfig:
    config_path = args.config.expanduser().resolve()
    benchmark_config = load_benchmark_config(config_path)
    eval_cfg = benchmark_config.get("evaluation", {})
    raw_fdr_threshold = eval_cfg.get("fdr_threshold", 0.05)
    return ReferenceAuditConfig(
        config_path=config_path,
        out_dir=args.out_dir.expanduser().resolve(),
        fdr_threshold=None if raw_fdr_threshold is None else float(raw_fdr_threshold),
        effect_threshold=float(eval_cfg.get("effect_threshold", 1.0)),
        protein_coding_only=(
            False if args.all_genes else _optional_bool(eval_cfg.get("protein_coding_only"), True)
        ),
        protein_coding_gtf=eval_cfg.get("protein_coding_gtf") or str(DEFAULT_GTF_REL_PATH),
        protein_coding_gene_cache=eval_cfg.get("protein_coding_gene_cache") or str(DEFAULT_CACHE_REL_PATH),
        sync_experiments=args.sync,
        expression_bins=args.expression_bins,
        write_row_table=args.write_row_table,
    )


def run_reference_audit(config: ReferenceAuditConfig) -> pathlib.Path:
    benchmark_config = load_benchmark_config(config.config_path)
    root = config.config_path.parent
    experiments_tsv = root / benchmark_config["experiments_tsv"]
    experiment_filters = benchmark_config.get("experiments")

    if config.sync_experiments:
        logger.info("Syncing selected experiment DE tables...")
        sync_zenodo_experiments(
            selected_experiment_paths(experiments_tsv, experiment_filters),
            repo=root,
        )

    experiments_df = pd.read_csv(experiments_tsv, sep="\t")
    if experiment_filters:
        experiments_df = filter_df(experiments_df, experiment_filters)
    if experiments_df.empty:
        raise ValueError("Experiment selection resolved to no rows.")

    protein_coding_gene_ids = None
    if config.protein_coding_only:
        logger.info("Loading protein-coding gene universe...")
        protein_coding_gene_ids = load_protein_coding_gene_ids(
            root=root,
            gtf_path=config.protein_coding_gtf,
            cache_path=config.protein_coding_gene_cache,
        )

    tables_dir = config.out_dir / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)
    config.out_dir.mkdir(parents=True, exist_ok=True)

    prepared_frames = []
    summaries = []
    missing = []
    for row in experiments_df.to_dict(orient="records"):
        dataset_id = str(row["id"])
        de_path = root / str(row["de_table_path"])
        if not de_path.exists():
            missing.append({
                "dataset_id": dataset_id,
                "path": str(de_path),
                "issue": "missing_de_table",
            })
            logger.warning("Missing DE table for %s: %s", dataset_id, de_path)
            continue

        de = read_de_table(de_path)
        prepared, summary = prepare_de_frame(
            de,
            dataset_id=dataset_id,
            mirna=str(row.get("mirna_name", "")),
            perturbation=str(row.get("experiment_type", "")),
            fdr_threshold=config.fdr_threshold,
            effect_threshold=config.effect_threshold,
            protein_coding_gene_ids=protein_coding_gene_ids,
        )
        summaries.append(summary)
        prepared_frames.append(prepared)
        if not summary["baseline_expression_column"]:
            missing.append({
                "dataset_id": dataset_id,
                "path": str(de_path),
                "issue": "missing_baseline_expression_column",
            })
        if not summary["available_uncertainty_columns"]:
            missing.append({
                "dataset_id": dataset_id,
                "path": str(de_path),
                "issue": "missing_uncertainty_columns",
            })

    if prepared_frames:
        combined = pd.concat(prepared_frames, ignore_index=True)
    else:
        combined = pd.DataFrame()

    pd.DataFrame(summaries).to_csv(tables_dir / "per_experiment_summary.tsv", sep="\t", index=False)
    pd.DataFrame(missing).to_csv(tables_dir / "missing_diagnostics.tsv", sep="\t", index=False)

    if not combined.empty:
        expression_bins = summarize_expression_bins(combined, bins=config.expression_bins)
        expression_bins.to_csv(tables_dir / "expression_bins.tsv", sep="\t", index=False)

        zero_categories = summarize_zero_categories(combined)
        zero_categories.to_csv(tables_dir / "zero_count_categories.tsv", sep="\t", index=False)

        gene_propensity = summarize_gene_response_propensity(combined)
        gene_propensity.to_csv(tables_dir / "gene_response_propensity.tsv", sep="\t", index=False)

        if config.write_row_table:
            combined.to_csv(tables_dir / "audit_rows.tsv.gz", sep="\t", index=False, compression="gzip")

    manifest = {
        "config": str(config.config_path),
        "out_dir": str(config.out_dir),
        "fdr_threshold": config.fdr_threshold,
        "effect_threshold": config.effect_threshold,
        "protein_coding_only": config.protein_coding_only,
        "ensembl_release": ENSEMBL_RELEASE if config.protein_coding_only else None,
        "experiment_count": int(len(experiments_df)),
        "processed_experiment_count": int(len(prepared_frames)),
        "outputs": {
            "per_experiment_summary": str(tables_dir / "per_experiment_summary.tsv"),
            "missing_diagnostics": str(tables_dir / "missing_diagnostics.tsv"),
            "expression_bins": str(tables_dir / "expression_bins.tsv"),
            "zero_count_categories": str(tables_dir / "zero_count_categories.tsv"),
            "gene_response_propensity": str(tables_dir / "gene_response_propensity.tsv"),
        },
    }
    (config.out_dir / "audit_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n",
        encoding="utf-8",
    )
    write_readme(config.out_dir, manifest)
    return config.out_dir


def write_readme(out_dir: pathlib.Path, manifest: dict) -> None:
    lines = [
        "# Reference Bias Audit",
        "",
        "This folder contains diagnostic tables for the current FuNmiRBench reference construction.",
        "It does not change benchmark labels or predictor evaluation.",
        "",
        "## Inputs",
        "",
        f"- Config: `{manifest['config']}`",
        f"- FDR threshold: `{manifest['fdr_threshold']}`",
        f"- Effect threshold: `{manifest['effect_threshold']}`",
        f"- Protein-coding only: `{manifest['protein_coding_only']}`",
        "",
        "## Tables",
        "",
        "- `tables/per_experiment_summary.tsv`: row counts, positives, and available diagnostic columns.",
        "- `tables/expression_bins.tsv`: positive rate and logFC summaries by baseline-expression bin.",
        "- `tables/zero_count_categories.tsv`: zero-baseline/treated categories when expression columns exist.",
        "- `tables/gene_response_propensity.tsv`: repeated gene-level response tendency across selected experiments.",
        "- `tables/missing_diagnostics.tsv`: missing input or missing diagnostic columns to resolve before deeper analysis.",
        "",
    ]
    out_dir.joinpath("README.md").write_text("\n".join(lines), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Audit FuNmiRBench reference labels against expression/provenance diagnostics."
    )
    parser.add_argument("--config", type=pathlib.Path, required=True)
    parser.add_argument(
        "--out-dir",
        type=pathlib.Path,
        default=pathlib.Path("results/reference_bias_audit"),
    )
    parser.add_argument(
        "--sync",
        action="store_true",
        help="Sync selected curated DE tables from Zenodo before auditing.",
    )
    parser.add_argument(
        "--all-genes",
        action="store_true",
        help="Disable the config's protein-coding-only gene universe filter.",
    )
    parser.add_argument("--expression-bins", type=int, default=5)
    parser.add_argument(
        "--write-row-table",
        action="store_true",
        help="Also write the per-row audit table as tables/audit_rows.tsv.gz.",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    setup_logging(parse_log_level(args.log_level))
    out_dir = run_reference_audit(build_audit_config(args))
    logger.info("Done. Audit outputs in %s", out_dir)


if __name__ == "__main__":
    main()
