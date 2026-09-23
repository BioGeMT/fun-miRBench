"""Configuration and metadata helpers for FuNmiRBench benchmark runs."""

from __future__ import annotations

import datetime as dt
import os
import pathlib
import urllib.parse

import pandas as pd

from funmirbench import DatasetMeta


def clean_optional_string(value, default=""):
    if value is None or pd.isna(value):
        return default
    text = str(value).strip()
    return text if text else default


def build_run_dir_name(*, run_date=None):
    """Return the date-based run directory name."""
    run_date = run_date or dt.date.today()
    if isinstance(run_date, dt.datetime):
        run_date = run_date.date()
    return run_date.strftime("%Y%m%d")


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
        parsed = urllib.parse.urlparse(str(row.get("gse_url", "") or ""))
        geo = urllib.parse.parse_qs(parsed.query).get("acc", [None])[0]
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


def selected_experiment_paths(tsv_path, filters) -> list[str]:
    df = pd.read_csv(tsv_path, sep="\t")
    if filters:
        df = filter_df(df, filters)
    return [str(value) for value in df["de_table_path"].tolist()]


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
