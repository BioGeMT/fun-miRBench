"""Download the complete published fun-miRBenvh benchmark dataset."""

from __future__ import annotations

import argparse
import logging
import pathlib

from funmirbench.experiment_store import sync_all_zenodo_experiments
from funmirbench.predictor_store import sync_all_zenodo_predictors
from funmirbench.logger import parse_log_level, setup_logging
from funmirbench.zenodo_store import ZENODO_RECORD


logger = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Download all published fun-miRBenvh experiments and standardized "
            f"predictors from Zenodo record {ZENODO_RECORD}."
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
    logger.info("Preparing published fun-miRBenvh data from Zenodo record %s.", ZENODO_RECORD)

    experiments = sync_all_zenodo_experiments(
        repo=args.repo,
        timeout=args.timeout,
        force=args.force,
    )
    predictors = sync_all_zenodo_predictors(
        repo=args.repo,
        timeout=args.timeout,
        force=args.force,
    )

    logger.info(
        "Prepared %d experiment tables and %d standardized predictor tables.",
        len(experiments),
        len(predictors),
    )


if __name__ == "__main__":
    main()
