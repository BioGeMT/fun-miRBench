#!/usr/bin/env python3
"""CLI entrypoint for the TargetScan standardization pipeline."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import (  # noqa: E402
    ENSEMBL_GTF_RELATIVE_PATH,
    MIRBASE_MATURE_RELATIVE_PATH,
    add_standard_pipeline_args,
    configure_file_logging,
    display_path,
    log_step,
    repo_root,
    resolve_cli_path,
)
from utils import (  # noqa: E402
    compute_final_statistics,
    parse_mirbase_mature,
    step1_download_targetscan_files,
    step2_build_representative_transcript_index,
    step3_download_ensembl115_gtf,
    step4_build_and_cache_ensembl115_tables,
    step5_qc_targetscan_vs_ensembl_transcripts,
    step6_download_mirbase_mature,
    step_build_human_mirna_annotations,
    step_build_human_mirna_family_members,
    step_write_standardized_predictions,
)


logger = logging.getLogger("pipeline")
FAMILY_EXPANDED_STEM = "targetscan_family_expanded"


def parse_args(root: Path, argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Standardize TargetScan predictions for fun-miRBench.")
    add_standard_pipeline_args(
        parser,
        tool_id="targetscan",
        root=root,
        include_resources_dir=True,
    )
    parser.add_argument(
        "--expand-family-members",
        action="store_true",
        help=(
            "Opt-in family-score projection: write each family score for every member miRNA "
            "of the family instead of only the TargetScan representative miRNA. Members get "
            "the family-level score, not an individual-miRNA score. Adds provenance columns "
            f"and, unless output paths are given, writes {FAMILY_EXPANDED_STEM}_standardized.tsv "
            "so the default standardized file is not overwritten."
        ),
    )
    args = parser.parse_args(argv)
    if args.expand_family_members:
        output_dir = Path(parser.get_default("standardized_output_file")).parent
        if args.standardized_output_file == parser.get_default("standardized_output_file"):
            args.standardized_output_file = output_dir / f"{FAMILY_EXPANDED_STEM}_standardized.tsv"
        if args.log_file == parser.get_default("log_file"):
            args.log_file = output_dir / f"{FAMILY_EXPANDED_STEM}_pipeline.log"
    return args


def main() -> None:
    root = repo_root()
    args = parse_args(root)
    args.common_resources_dir = resolve_cli_path(args.common_resources_dir, root)
    args.data_dir = resolve_cli_path(args.data_dir, root)
    args.resources_dir = resolve_cli_path(args.resources_dir, root)
    args.standardized_output_file = resolve_cli_path(args.standardized_output_file, root)
    args.log_file = resolve_cli_path(args.log_file, root)

    configure_file_logging(args.log_file, args.log_level)
    logger.info("Starting TargetScan standardization pipeline")

    total_steps = 8

    log_step(logger, 1, total_steps, "Download and unzip TargetScan inputs")
    files = step1_download_targetscan_files(
        args.data_dir,
        args.resources_dir,
        force=False,
    )

    log_step(logger, 2, total_steps, "Build representative-transcript index")
    tx_index = step2_build_representative_transcript_index(files["Gene_info.txt"], species_id="9606")

    log_step(logger, 3, total_steps, "Download Ensembl v115 GTF")
    ensembl_gtf = step3_download_ensembl115_gtf(
        args.common_resources_dir / ENSEMBL_GTF_RELATIVE_PATH,
        force=False,
    )

    log_step(logger, 4, total_steps, "Build/cache Ensembl v115 tables")
    ensembl_tables = step4_build_and_cache_ensembl115_tables(
        ensembl_gtf,
        cache_dir=args.resources_dir,
        force_rebuild=False,
    )

    log_step(logger, 5, total_steps, "Run TargetScan-vs-Ensembl transcript QC")
    step5_qc_targetscan_vs_ensembl_transcripts(
        tx_index=tx_index,
        ensembl_tables=ensembl_tables,
    )

    log_step(logger, 6, total_steps, "Download and parse miRBase mature annotations")
    mirbase_fa = step6_download_mirbase_mature(
        args.common_resources_dir / MIRBASE_MATURE_RELATIVE_PATH,
        force=False,
    )
    mirbase_acc2name = parse_mirbase_mature(mirbase_fa)

    log_step(logger, 7, total_steps, "Build human miRNA annotations")
    mirna_annotations = step_build_human_mirna_annotations(
        files["miR_Family_Info.txt"],
        mirbase_acc2name=mirbase_acc2name,
        species_id="9606",
    )

    family_members = None
    if args.expand_family_members:
        family_members = step_build_human_mirna_family_members(
            files["miR_Family_Info.txt"],
            species_id="9606",
        )

    log_step(logger, 8, total_steps, "Write standardized TargetScan predictions")
    step_write_standardized_predictions(
        files["Summary_Counts.all_predictions.txt"],
        tx_index=tx_index,
        ensembl_tables=ensembl_tables,
        mirna_annotations=mirna_annotations,
        output_path=args.standardized_output_file,
        species_id="9606",
        family_members=family_members,
    )

    compute_final_statistics(args.standardized_output_file)
    logger.info("Relative output path: %s", display_path(args.standardized_output_file))


if __name__ == "__main__":
    main()
