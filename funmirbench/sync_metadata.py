"""Sync candidate experiment metadata rows into the experiment registry."""

from __future__ import annotations

import argparse
import logging
import pathlib
import re
import sys

import pandas as pd

from funmirbench.logger import parse_log_level, setup_logging


logger = logging.getLogger(__name__)

REGISTRY_PATH = pathlib.Path("metadata") / "mirna_experiment_info.tsv"
CANDIDATE_PATTERN = "candidate_metadata.tsv"
REGISTRY_KEY = "id"


def repo_root() -> pathlib.Path:
    return pathlib.Path(__file__).resolve().parents[1]


def read_tsv(path: pathlib.Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def read_candidate_metadata_with_source(path: pathlib.Path) -> pd.DataFrame:
    df = read_tsv(path)
    df["_source_path"] = str(path)
    df["_source_mtime"] = path.stat().st_mtime
    return df


def collect_input_paths(inputs: list[pathlib.Path], repo: pathlib.Path) -> list[pathlib.Path]:
    if not inputs:
        base = repo / "pipelines" / "experiments" / "runs"
        if not base.exists():
            return []
        return sorted(base.glob(f"*/{CANDIDATE_PATTERN}"))

    paths = []
    seen: set[str] = set()
    for item in inputs:
        resolved = item.expanduser().resolve()
        candidates = sorted(resolved.rglob(CANDIDATE_PATTERN)) if resolved.is_dir() else [resolved]
        for path in candidates:
            if str(path) not in seen:
                paths.append(path)
                seen.add(str(path))
    return paths


def _normalize_pubmed_id(value) -> str:
    text = str(value or "").strip()
    if not text or text.upper() == "NA" or text.rstrip("/").upper().endswith("/NA"):
        return ""
    match = re.search(r"(\d{5,})", text)
    return match.group(1) if match else text


def _prepare_canonical_columns(incoming: pd.DataFrame) -> pd.DataFrame:
    incoming = incoming.copy()

    if "pubmed_id" not in incoming.columns and "article_pubmed_id" in incoming.columns:
        incoming["pubmed_id"] = incoming["article_pubmed_id"].map(_normalize_pubmed_id)
    elif "pubmed_id" in incoming.columns:
        incoming["pubmed_id"] = incoming["pubmed_id"].map(_normalize_pubmed_id)

    if "geo_accession" not in incoming.columns:
        incoming["geo_accession"] = ""
    if "gse_url" in incoming.columns:
        missing_geo = incoming["geo_accession"].astype(str).str.strip().eq("")
        extracted = incoming.loc[missing_geo, "gse_url"].astype(str).str.extract(
            r"[?&]acc=(GSE\d+)",
            expand=False,
        )
        incoming.loc[missing_geo, "geo_accession"] = extracted.fillna("")

    return incoming


def merge_registry(existing: pd.DataFrame, incoming: pd.DataFrame) -> pd.DataFrame:
    if REGISTRY_KEY not in existing.columns:
        raise ValueError(f"Registry is missing key column {REGISTRY_KEY!r}.")
    if REGISTRY_KEY not in incoming.columns:
        raise ValueError(f"Incoming rows are missing key column {REGISTRY_KEY!r}.")

    for column in existing.columns:
        if column not in incoming.columns:
            incoming[column] = ""
    incoming = incoming[existing.columns].copy()

    existing_by_key = existing.set_index(REGISTRY_KEY)
    for column in ("control_samples", "condition_samples"):
        if column not in existing.columns:
            continue
        for index, key in incoming[REGISTRY_KEY].items():
            value = incoming.at[index, column]
            if (pd.isna(value) or str(value).strip() == "") and key in existing_by_key.index:
                incoming.at[index, column] = existing_by_key.at[key, column]

    incoming_keys = incoming[REGISTRY_KEY].tolist()
    if len(set(incoming_keys)) != len(incoming_keys):
        raise ValueError(f"Incoming metadata contains duplicate {REGISTRY_KEY} values: {incoming_keys}")

    existing = existing[~existing[REGISTRY_KEY].isin(incoming[REGISTRY_KEY])].copy()
    return pd.concat([existing, incoming], ignore_index=True)


def sync_metadata(
    *,
    inputs: list[pathlib.Path],
    repo: pathlib.Path | None = None,
    registry_path: pathlib.Path | None = None,
) -> dict:
    repo = (repo or repo_root()).resolve()
    registry = (registry_path or (repo / REGISTRY_PATH)).resolve()
    input_paths = collect_input_paths(inputs, repo)
    if not input_paths:
        raise ValueError("No candidate metadata TSV files found.")

    existing = read_tsv(registry)
    incoming_frames = [read_candidate_metadata_with_source(path) for path in input_paths]
    incoming = pd.concat(incoming_frames, ignore_index=True)

    if REGISTRY_KEY in incoming.columns:
        duplicate_keys = incoming[REGISTRY_KEY][incoming[REGISTRY_KEY].duplicated(keep=False)].unique().tolist()
        if duplicate_keys:
            logger.warning(
                "Found duplicate incoming %s values; keeping the most recent candidate_metadata.tsv for each: %s",
                REGISTRY_KEY,
                duplicate_keys,
            )
            incoming = (
                incoming.sort_values("_source_mtime")
                .drop_duplicates(subset=[REGISTRY_KEY], keep="last")
                .reset_index(drop=True)
            )

    incoming = incoming.drop(columns=["_source_path", "_source_mtime"], errors="ignore")
    incoming = _prepare_canonical_columns(incoming)
    merged = merge_registry(existing, incoming)
    merged.to_csv(registry, sep="\t", index=False)
    return {
        "registry": str(registry),
        "rows_before": int(existing.shape[0]),
        "rows_added_or_updated": int(incoming.shape[0]),
        "rows_after": int(merged.shape[0]),
        "inputs": [str(path) for path in input_paths],
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Sync candidate experiment metadata rows into the registry.")
    parser.add_argument("--input", type=pathlib.Path, action="append", default=[])
    parser.add_argument("--registry", type=pathlib.Path)
    parser.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"])
    return parser.parse_args(sys.argv[1:])


def main() -> None:
    args = parse_args()
    setup_logging(parse_log_level(args.log_level))
    result = sync_metadata(inputs=args.input, registry_path=args.registry)
    logger.info("Registry: %s", result["registry"])
    logger.info("Rows before: %s", result["rows_before"])
    logger.info("Rows added or updated: %s", result["rows_added_or_updated"])
    logger.info("Rows after: %s", result["rows_after"])


if __name__ == "__main__":
    main()
