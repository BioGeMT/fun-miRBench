"""Download the complete published FuNmiRBench benchmark dataset."""

from __future__ import annotations

import argparse
import pathlib

from funmirbench.experiment_store import sync_all_zenodo_experiments
from funmirbench.predictor_store import sync_all_zenodo_predictors
from funmirbench.zenodo_store import ZENODO_RECORD


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Download all published FuNmiRBench experiments and standardized "
            f"predictors from Zenodo record {ZENODO_RECORD}."
        )
    )
    parser.add_argument("--repo", type=pathlib.Path, default=None)
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

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

    print(
        f"Prepared {len(experiments)} experiment tables and "
        f"{len(predictors)} standardized predictor tables "
        f"from Zenodo record {ZENODO_RECORD}."
    )


if __name__ == "__main__":
    main()
