"""Shared helpers for standardized predictor pipeline entrypoints.

Keep predictor-specific parsing and biological mapping logic inside each predictor's
``utils.py``. This module is intentionally small: it standardizes repository path
resolution, logging setup, progress messages, and final standardized TSV writing.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Any, Iterable


LOG_LEVEL_CHOICES = ["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]
STANDARDIZED_COLUMNS = ["Ensembl_ID", "Gene_Name", "miRNA_ID", "miRNA_Name", "Score"]
MIRBASE_MATURE_RELATIVE_PATH = Path("mirbase") / "mature.fa"
ENSEMBL_GTF_RELATIVE_PATH = Path("ensembl") / "Homo_sapiens.GRCh38.115.gtf.gz"


def repo_root() -> Path:
    """Return the repository root from this shared module location."""
    return Path(__file__).resolve().parents[2]


ROOT = repo_root()
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from funmirbench.logger import parse_log_level, setup_logging  # noqa: E402


def predictor_dir(tool_id: str, *, root: Path | None = None) -> Path:
    """Return ``pipelines/standardized_predictors/<tool_id>``."""
    return (root or ROOT) / "pipelines" / "standardized_predictors" / tool_id


def common_resources_dir(*, root: Path | None = None) -> Path:
    """Return the shared downloaded-annotation cache directory."""
    return (root or ROOT) / "data" / "common_resources"


def resolve_cli_path(path: str | Path, root: Path | None = None) -> Path:
    """Resolve CLI paths relative to the repository root unless absolute."""
    value = Path(path)
    if value.is_absolute():
        return value
    return (root or ROOT) / value


def display_path(path: str | Path, *, root: Path | None = None) -> Path:
    """Return a repository-relative path for portable log messages."""
    value = Path(path)
    if not value.is_absolute():
        return value

    try:
        return value.resolve().relative_to((root or ROOT).resolve())
    except ValueError:
        return value


def add_log_level_arg(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=LOG_LEVEL_CHOICES,
        help="Logging level. Default: INFO",
    )


def add_standard_pipeline_args(
    parser: argparse.ArgumentParser,
    *,
    tool_id: str,
    root: Path,
    include_resources_dir: bool,
) -> None:
    """Add harmonized path and logging arguments for a predictor pipeline."""
    pipeline_dir = predictor_dir(tool_id, root=root)
    predictions_dir = root / "data" / "predictions" / tool_id
    parser.add_argument(
        "--common-resources-dir",
        type=Path,
        default=common_resources_dir(root=root),
        help="Shared miRBase and Ensembl resource directory",
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=pipeline_dir / "data",
        help="Directory for downloaded predictor data",
    )
    if include_resources_dir:
        parser.add_argument(
            "--resources-dir",
            type=Path,
            default=pipeline_dir / "resources",
            help="Directory for predictor-specific supporting resources",
        )
    parser.add_argument(
        "--standardized-output-file",
        type=Path,
        default=predictions_dir / f"{tool_id}_standardized.tsv",
        help="Output standardized TSV file",
    )
    parser.add_argument(
        "--log-file",
        type=Path,
        default=predictions_dir / f"{tool_id}_pipeline.log",
        help="Log file path",
    )
    add_log_level_arg(parser)


def configure_file_logging(log_file: Path, log_level: str) -> None:
    """Configure console logging plus one predictor-specific file handler."""
    log_file.parent.mkdir(parents=True, exist_ok=True)
    setup_logging(parse_log_level(log_level))
    root_logger = logging.getLogger()

    # Avoid duplicate file handlers when multiple pipelines run in one process.
    log_file = log_file.resolve()
    for handler in root_logger.handlers:
        if isinstance(handler, logging.FileHandler) and Path(handler.baseFilename) == log_file:
            return

    file_handler = logging.FileHandler(log_file, mode="w", encoding="utf-8")
    file_handler.setLevel(root_logger.level)
    file_handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(name)s: %(message)s"))
    root_logger.addHandler(file_handler)


def log_step(logger: logging.Logger, step_number: int, total_steps: int, message: str) -> None:
    logger.info("Step %d/%d: %s", step_number, total_steps, message)


def validate_standardized_table(
    df: Any,
    *,
    required_columns: Iterable[str] = STANDARDIZED_COLUMNS,
) -> None:
    """Validate the common standardized predictor output schema.

    The check is intentionally schema-only. Some predictors may legitimately have
    blank gene names or MIMAT IDs before future curation passes, so row-level QC
    should stay predictor-specific instead of blocking shared output writing.
    """
    required = list(required_columns)
    missing = [column for column in required if column not in df.columns]
    if missing:
        raise ValueError(f"Standardized predictor output is missing required columns: {missing}")


def write_standardized_table(
    df: Any,
    output: Path,
    *,
    logger: logging.Logger | None = None,
    columns: Iterable[str] = STANDARDIZED_COLUMNS,
) -> None:
    """Validate and write a standardized predictor TSV."""
    output.parent.mkdir(parents=True, exist_ok=True)
    columns = list(columns)
    validate_standardized_table(df, required_columns=columns)
    df.loc[:, columns].to_csv(output, sep="\t", index=False)
    if logger is not None:
        logger.info("Output written to: %s", display_path(output))
        logger.info("Rows written: %d", len(df))
