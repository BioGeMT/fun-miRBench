"""Configuration and metadata helpers for fun-miRBenvh benchmark runs."""

from __future__ import annotations

import datetime as dt
import difflib
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


def build_run_dir_name(*, run_time=None):
    """Return the timestamp-based run directory name."""
    run_time = run_time or dt.datetime.now()
    if isinstance(run_time, dt.date) and not isinstance(run_time, dt.datetime):
        run_time = dt.datetime.combine(run_time, dt.time())
    return run_time.strftime("%Y%m%d_%H%M%S")


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

    if filters and "tool_id" in filters:
        requested = filters["tool_id"]
        if not isinstance(requested, list):
            requested = [requested]

        available = [str(value) for value in df["tool_id"].dropna().tolist()]
        unknown = [str(value) for value in requested if str(value) not in available]
        if unknown:
            lines = ["Unknown predictor tool_id value(s) in benchmark config:", ""]
            for tool_id in unknown:
                lines.append(f"- {tool_id}")
                matches = difflib.get_close_matches(tool_id, available, n=1, cutoff=0.5)
                if matches:
                    lines.append(f"  Did you mean: {matches[0]}?")
            lines.extend(
                [
                    "",
                    "Registered predictor IDs:",
                    *[f"- {tool_id}" for tool_id in available],
                    "",
                    f"Predictor registry: {tsv_path}",
                ]
            )
            raise ValueError("\n".join(lines))

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
