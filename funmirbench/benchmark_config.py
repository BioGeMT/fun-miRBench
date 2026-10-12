"""Configuration and metadata helpers for fun-miRBench benchmark runs."""

from __future__ import annotations

import datetime as dt
import os
import pathlib

import pandas as pd

from funmirbench import DatasetMeta


def clean_optional_string(value, default=""):
    if value is None or pd.isna(value):
        return default
    text = str(value).strip()
    return text if text else default


def build_run_dir_name(
    *,
    experiments,
    tool_ids,
    eval_cfg,
    tags=None,
    run_date=None,
):
    """Return the timestamp-based run directory name."""
    del experiments, tool_ids, eval_cfg, tags
    run_timestamp = run_date or dt.datetime.now()
    if isinstance(run_timestamp, dt.datetime):
        return run_timestamp.strftime("%Y%m%d_%H%M%S")
    return run_timestamp.strftime("%Y%m%d")


def filter_df(df, filters):
    """AND across columns, OR within each column's value list."""
    for col, values in filters.items():
        if col not in df.columns:
            raise ValueError(
                f"Filter column {col!r} was not found. Available columns: {sorted(df.columns)}"
            )
        if not isinstance(values, list):
            values = [values]
        df = df[df[col].isin(values)]
    return df


def load_experiments(tsv_path, root, filters):
    df = pd.read_csv(tsv_path, sep="\t")
    if filters:
        df = filter_df(df, filters)

    metas = []
    for _, row in df.iterrows():
        geo = clean_optional_string(row.get("geo_accession"), default="")
        metas.append(
            DatasetMeta(
                id=str(row["id"]),
                miRNA=str(row["mirna_name"]),
                cell_line=clean_optional_string(row.get("tested_cell_line")),
                tissue=clean_optional_string(row.get("tissue")),
                perturbation=clean_optional_string(row.get("experiment_type")),
                organism=clean_optional_string(row.get("organism")),
                geo_accession=geo,
                data_path=str(row["de_table_path"]),
                root=root,
            )
        )
    return metas


def load_predictions(tsv_path, filters):
    df = pd.read_csv(tsv_path, sep="\t")
    if filters:
        df = filter_df(df, filters)
    if df["tool_id"].duplicated().any():
        raise ValueError("Duplicate tool_id values found after predictor filtering.")
    return {row["tool_id"]: row.to_dict() for _, row in df.iterrows()}


def resolve_predictor_output_path(predictor_output_path, root):
    path = pathlib.Path(os.path.expandvars(str(predictor_output_path))).expanduser()
    if not path.is_absolute():
        path = pathlib.Path(root) / path
    return path


def validate_predictor_output_files(predictions, root):
    missing = []
    for tool_id, meta in predictions.items():
        configured_path = meta.get("predictor_output_path")
        if configured_path is None or pd.isna(configured_path) or str(configured_path).strip() == "":
            missing.append((tool_id, configured_path))
            continue
        resolved_path = resolve_predictor_output_path(configured_path, root)
        if not resolved_path.is_file():
            missing.append((tool_id, configured_path))

    if not missing:
        return

    lines = ["Missing predictor output files:", ""]
    lines.extend(f"{tool_id}: {configured_path}" for tool_id, configured_path in missing)
    lines.extend(
        [
            "",
            (
                "Generate or download the selected standardized predictor files, "
                "or remove unavailable predictors from benchmark.yaml."
            ),
        ]
    )
    raise FileNotFoundError("\n".join(lines))
