"""Cache the final per-gene tables and the family-expanded TargetScan scores from their Zenodo record."""

from __future__ import annotations

import gzip
import logging
import pathlib
import shutil
import tempfile
import zipfile

import pandas as pd

from funmirbench.benchmark_config import resolve_predictor_output_path
from funmirbench.predictor_store import _download_archive, repo_root
from funmirbench.zenodo_store import FINAL_DATA_RECORD, fetch_zenodo_file_registry

logger = logging.getLogger(__name__)

FINAL_DATA_ARCHIVE = "final_data.zip"
FINAL_DATA_DIR = pathlib.Path("results/final_data")
FAMILY_EXPANDED_TOOL = "targetscan_family_expanded"
FAMILY_EXPANDED_ARCHIVE = "targetscan_family_expanded_standardized.tsv.gz"


def _registry(registry: dict[str, dict] | None, timeout: int) -> dict[str, dict]:
    return registry if registry is not None else fetch_zenodo_file_registry(record=FINAL_DATA_RECORD, timeout=timeout)


def sync_final_data(
    *,
    repo: pathlib.Path | None = None,
    registry: dict[str, dict] | None = None,
    timeout: int = 120,
    force: bool = False,
) -> pathlib.Path:
    """Ensure results/final_data is available locally; an existing folder is kept unless force is set."""
    destination = (repo or repo_root()).resolve() / FINAL_DATA_DIR
    if destination.exists() and not force:
        logger.info("Final tables already present: %s", destination)
        return destination
    destination.parent.mkdir(parents=True, exist_ok=True)
    archive = _download_archive(_registry(registry, timeout)[FINAL_DATA_ARCHIVE], directory=destination.parent,
                                timeout=timeout)
    try:
        with tempfile.TemporaryDirectory(dir=destination.parent, prefix=".final_data.") as unpacked:
            with zipfile.ZipFile(archive) as handle:
                handle.extractall(unpacked)
            if destination.exists():
                shutil.rmtree(destination)
            (pathlib.Path(unpacked) / "final_data").replace(destination)
    finally:
        archive.unlink(missing_ok=True)
    return destination


def sync_family_expanded_targetscan(
    *,
    repo: pathlib.Path | None = None,
    registry: dict[str, dict] | None = None,
    timeout: int = 120,
    force: bool = False,
) -> pathlib.Path:
    """Ensure the family-expanded TargetScan scores registered in metadata/predictions_info.tsv are available."""
    repo = (repo or repo_root()).resolve()
    tools = pd.read_csv(repo / "metadata" / "predictions_info.tsv", sep="\t").set_index("tool_id")
    destination = resolve_predictor_output_path(tools.loc[FAMILY_EXPANDED_TOOL, "predictor_output_path"], repo)
    if destination.is_file() and not force:
        logger.info("Family-expanded TargetScan scores already present: %s", destination)
        return destination
    destination.parent.mkdir(parents=True, exist_ok=True)
    archive = _download_archive(_registry(registry, timeout)[FAMILY_EXPANDED_ARCHIVE], directory=destination.parent,
                                timeout=timeout)
    partial = destination.with_name(f".{destination.name}.tmp")
    try:
        with gzip.open(archive, "rb") as source, partial.open("wb") as target:
            shutil.copyfileobj(source, target)
        partial.replace(destination)
    finally:
        archive.unlink(missing_ok=True)
        partial.unlink(missing_ok=True)
    return destination
