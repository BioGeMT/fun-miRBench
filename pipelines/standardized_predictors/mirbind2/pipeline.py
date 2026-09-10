#!/usr/bin/env python3

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import (  # noqa: E402
    MIRBASE_MATURE_RELATIVE_PATH,
    add_standard_pipeline_args,
    configure_file_logging,
    log_step,
    repo_root,
    resolve_cli_path,
    write_standardized_table,
)
from utils import (  # noqa: E402
    build_output_table,
    create_mirna_name_to_mimat_mapping,
    download_file,
    load_predictions,
    map_mirna_names_to_mimat,
    resolve_path_relative_to_root,
)

logger = logging.getLogger("pipeline")


def parse_args(root: Path) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Standardize miRBind2 predictions for FuNmiRBench.")
    add_standard_pipeline_args(
        parser,
        tool_id="mirbind2",
        root=root,
        include_resources_dir=False,
    )
    return parser.parse_args()


def main() -> None:
    root = repo_root()
    args = parse_args(root)
    args.common_resources_dir = resolve_cli_path(args.common_resources_dir, root)
    args.data_dir = resolve_cli_path(args.data_dir, root)
    args.standardized_output_file = resolve_cli_path(args.standardized_output_file, root)
    args.log_file = resolve_cli_path(args.log_file, root)

    configure_file_logging(args.log_file, args.log_level)
    logger.info("Starting miRBind2 standardization pipeline")

    raw_ensembl_column = "Gene_ID"
    raw_gene_name_column = "Gene_Symbol"
    raw_mirna_name_column = "miRNA_Name"
    raw_prediction_column = "miRBind2_3UTR_prediction"
    raw_columns = [
        "Transcript_ID",
        raw_ensembl_column,
        raw_gene_name_column,
        raw_mirna_name_column,
        raw_prediction_column,
    ]

    final_columns = ["Ensembl_ID", "Gene_Name", "miRNA_ID", "miRNA_Name", "Score"]

    mirbind2_predictions_url = "https://zenodo.org/records/20609975/files/3utrs_mirbind2_predictions.tsv.gz?download=1"
    total_steps = 6

    log_step(logger, 1, total_steps, "Download or resolve shared miRBase and raw miRBind2 predictions")
    mirbase_url = "https://mirbase.org/download_version_files/22.1/mature.fa"
    mirbase_path = download_file(
        mirbase_url,
        args.common_resources_dir / MIRBASE_MATURE_RELATIVE_PATH,
        resource_label="miRBase mature.fa resource",
    )
    raw_predictions_path = download_file(
        mirbind2_predictions_url,
        args.data_dir / "3utrs_mirbind2_predictions.tsv.gz",
        timeout=360,
        resource_label="miRBind2 raw prediction file",
    )

    log_step(logger, 2, total_steps, "Load raw miRBind2 predictions")
    pred_df = load_predictions(
        raw_predictions_path,
        required_columns=raw_columns,
        ensembl_column=raw_ensembl_column,
        gene_name_column=raw_gene_name_column,
        mirna_name_column=raw_mirna_name_column,
        score_column=raw_prediction_column,
    )

    log_step(logger, 3, total_steps, "Create miRNA name-to-MIMAT mapping from miRBase mature.fa")
    mirna_name_to_mimat = create_mirna_name_to_mimat_mapping(mirbase_path)

    log_step(logger, 4, total_steps, "Annotate raw miRNA names to MIMAT IDs")
    pred_df = map_mirna_names_to_mimat(
        pred_df,
        mirna_name_to_mimat,
        mirna_name_column=raw_mirna_name_column,
        mimat_column="miRNA_ID",
    )

    log_step(logger, 5, total_steps, "Build standardized output columns")
    final_df = build_output_table(
        pred_df,
        raw_ensembl_column=raw_ensembl_column,
        raw_gene_name_column=raw_gene_name_column,
        raw_mirna_name_column=raw_mirna_name_column,
        raw_prediction_column=raw_prediction_column,
        ensembl_column="Ensembl_ID",
        gene_name_column="Gene_Name",
        mimat_column="miRNA_ID",
        mirna_name_column="miRNA_Name",
        score_column="Score",
        final_columns=final_columns,
    )

    log_step(logger, 6, total_steps, "Validate and write standardized output table")
    write_standardized_table(final_df, args.standardized_output_file, logger=logger, columns=final_columns)
    logger.info("Relative output path: %s", resolve_path_relative_to_root(args.standardized_output_file))


if __name__ == "__main__":
    main()
