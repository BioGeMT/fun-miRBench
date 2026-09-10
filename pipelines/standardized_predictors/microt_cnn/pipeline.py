#!/usr/bin/env python3
"""CLI entrypoint for the microT-CNN standardization pipeline."""

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
    log_step,
    repo_root,
    resolve_cli_path,
    write_standardized_table,
)
from utils import (  # noqa: E402
    build_ensembl_tx_to_gene_from_gtf,
    build_output_table,
    collapse_transcript_rows_to_genes,
    create_mirna_name_to_mimat_mapping,
    download_file,
    load_ensembl_tx_to_gene,
    load_prediction_files,
    map_mirna_names_to_mimat,
    map_transcripts_to_genes,
    resolve_path_relative_to_root,
)

logger = logging.getLogger("pipeline")


def parse_args(root: Path) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Standardize microT-CNN predictions for FuNmiRBench.")
    add_standard_pipeline_args(
        parser,
        tool_id="microt_cnn",
        root=root,
        include_resources_dir=True,
    )
    return parser.parse_args()


def main() -> None:
    root = repo_root()
    args = parse_args(root)
    args.common_resources_dir = resolve_cli_path(args.common_resources_dir, root)
    args.data_dir = resolve_cli_path(args.data_dir, root)
    args.resources_dir = resolve_cli_path(args.resources_dir, root)
    args.standardized_output_file = resolve_cli_path(args.standardized_output_file, root)
    args.log_file = resolve_cli_path(args.log_file, root)

    configure_file_logging(args.log_file, args.log_level)
    logger.info("Starting microT-CNN standardization pipeline")
    total_steps = 7

    log_step(logger, 1, total_steps, "Resolve raw microT-CNN predictions and external resources")
    mirbase_url = "https://mirbase.org/download_version_files/22.1/mature.fa"
    mirbase_path = download_file(
        mirbase_url,
        args.common_resources_dir / MIRBASE_MATURE_RELATIVE_PATH,
        resource_label="miRBase mature.fa resource",
    )

    microt_cnn_predictions_url = "10.5281/zenodo.20313523"
    raw_predictions_path = download_file(
        microt_cnn_predictions_url,
        args.data_dir / "microT_CNN_prediction_result_human_all_scores_gene_level.tsv.gz",
        timeout=360,
        resource_label="microT-CNN all-score prediction file",
    )
    raw_gene_id_column = "Ensembl_ID"
    raw_gene_name_column = "Gene_Name"
    raw_mirna_column = "miRNA_Name"
    raw_prediction_column = "Score"

    log_step(logger, 2, total_steps, "Load raw microT-CNN predictions")
    pred_df = load_prediction_files(
        raw_predictions_path,
        raw_gene_id_column,
        raw_gene_name_column,
        raw_mirna_column,
        raw_prediction_column,
    )

    log_step(logger, 3, total_steps, "Create miRNA and transcript mapping resources")
    mirna_name_to_mimat_map = create_mirna_name_to_mimat_mapping(mirbase_path)
    tx2gene_file = args.resources_dir / "ensembl115_tx2gene.tsv.gz"
    if tx2gene_file.exists():
        logger.info("Loading Ensembl transcript-to-gene mapping: %s", resolve_path_relative_to_root(tx2gene_file))
        tx_to_gene = load_ensembl_tx_to_gene(tx2gene_file)
    else:
        ensembl_gtf_url = "https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz"
        ensembl_gtf_path = download_file(
            ensembl_gtf_url,
            args.common_resources_dir / ENSEMBL_GTF_RELATIVE_PATH,
            timeout=360,
            resource_label="Ensembl v115 GTF resource",
        )
        tx_to_gene = build_ensembl_tx_to_gene_from_gtf(ensembl_gtf_path, tx2gene_file)

    mimat_column = "miRNA_ID"
    ensembl_id_column = "Ensembl_ID"
    gene_name_column = "Gene_Name"
    mirna_name_column = "miRNA_Name"
    score_column = "Score"
    final_columns = [ensembl_id_column, gene_name_column, mimat_column, mirna_name_column, score_column]

    log_step(logger, 4, total_steps, "Map prediction miRNA names to MIMAT IDs")
    pred_df = map_mirna_names_to_mimat(
        pred_df,
        mirna_name_to_mimat_map,
        raw_mirna_column,
        mirna_name_column,
        mimat_column,
    )

    log_step(logger, 5, total_steps, "Map microT-CNN Ensembl transcript IDs to Ensembl gene IDs")
    pred_df = map_transcripts_to_genes(
        pred_df,
        tx_to_gene,
        raw_gene_id_column,
        ensembl_id_column,
    )

    log_step(logger, 6, total_steps, "Build standardized output columns")
    final_df = build_output_table(
        pred_df,
        raw_prediction_column,
        score_column,
        final_columns,
    )
    final_df = collapse_transcript_rows_to_genes(
        final_df,
        ensembl_id_column,
        gene_name_column,
        mirna_name_column,
        mimat_column,
        score_column,
    )

    log_step(logger, 7, total_steps, "Validate and write standardized output table")
    write_standardized_table(final_df, args.standardized_output_file, logger=logger, columns=final_columns)
    logger.info("Relative output path: %s", resolve_path_relative_to_root(args.standardized_output_file))


if __name__ == "__main__":
    main()
