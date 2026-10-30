"""Download the published final tables and standardized predictors of the fun-miRBench benchmark."""

from __future__ import annotations

import argparse
import logging
import pathlib

from funmirbench.final_data_store import sync_family_expanded_targetscan, sync_final_data
from funmirbench.predictor_store import sync_all_zenodo_predictors
from funmirbench.logger import parse_log_level, setup_logging
from funmirbench.zenodo_store import FINAL_DATA_RECORD, ZENODO_RECORD


logger = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            f"Download the final tables and family-expanded TargetScan scores from Zenodo record {FINAL_DATA_RECORD} "
            f"and the other standardized predictors from Zenodo record {ZENODO_RECORD}."
        )
    )
    parser.add_argument("--repo", type=pathlib.Path, default=None)
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument("--force", action="store_true")
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    setup_logging(parse_log_level(args.log_level))
    logger.info("Preparing published fun-miRBench data from Zenodo records %s and %s.", FINAL_DATA_RECORD, ZENODO_RECORD)

    final_data = sync_final_data(repo=args.repo, timeout=args.timeout, force=args.force)
    predictors = sync_all_zenodo_predictors(
        repo=args.repo,
        timeout=args.timeout,
        force=args.force,
    )

    expanded = sync_family_expanded_targetscan(repo=args.repo, timeout=args.timeout, force=args.force)

    logger.info("Prepared %s and %d standardized predictor tables.", final_data, len(predictors) + 1)
    logger.info("Family-expanded TargetScan scores: %s", expanded)


if __name__ == "__main__":
    main()
