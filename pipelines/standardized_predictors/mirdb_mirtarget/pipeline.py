#!/usr/bin/env python3
"""CLI entrypoint for the miRDB standardization pipeline."""

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
    create_ncbi_gene_id_to_ensembl_mapping,
    create_refseq_to_ensembl_mapping,
    download_file,
    fill_unmapped_rows_with_refseq_to_ensembl,
    load_prediction_files,
    map_mirna_names_to_mimat,
    map_ncbi_gene_id_to_ensembl,
    resolve_path_relative_to_root,
)

logger = logging.getLogger("pipeline")


def parse_args(root: Path) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Standardize miRDB predictions for FuNmiRBench.")
    add_standard_pipeline_args(
        parser,
        tool_id="mirdb_mirtarget",
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
    logger.info("Starting miRDB standardization pipeline")
    total_steps = 9

    log_step(logger, 1, total_steps, "Resolve raw miRDB predictions and external miRBase/BioMart resources")
    mirbase_url = "https://mirbase.org/download_version_files/22.1/mature.fa"
    mirbase_path = download_file(
        mirbase_url,
        args.common_resources_dir / MIRBASE_MATURE_RELATIVE_PATH,
        resource_label="miRBase mature.fa resource",
    )

    biomart_query = """<?xml version="1.0" encoding="UTF-8"?>
        <!DOCTYPE Query>
        <Query  virtualSchemaName = "default" formatter = "TSV" header = "1" uniqueRows = "1" count = "" datasetConfigVersion = "0.6" >
            <Dataset name = "hsapiens_gene_ensembl" interface = "default" >
                <Attribute name = "ensembl_gene_id" />
                <Attribute name = "external_gene_name" />
                <Attribute name = "refseq_mrna" />
                <Attribute name = "entrezgene_id" />
            </Dataset>
        </Query>"""
    biomart_url = "https://Sep2025.archive.ensembl.org/biomart/martservice"
    biomart_path = download_file(
        biomart_url,
        args.resources_dir / "hsapiens_ncbi_gene_id_refseq_to_ensembl.tsv",
        params={"query": biomart_query},
        timeout=360,
        resource_label="BioMart NCBI Gene ID/RefSeq-to-Ensembl mapping table",
    )

    mirdb_predictions_url = "https://mirdb.org/download/miRDB_v6.0_prediction_result_human_all_scores.txt.gz"
    raw_predictions_path = download_file(
        mirdb_predictions_url,
        args.data_dir / "miRDB_v6.0_prediction_result_human_all_scores.txt.gz",
        timeout=360,
        resource_label="miRDB v6.0 all-score prediction file",
    )

    raw_mirna_column = "miRNA"
    raw_transcript_column = "refseq_id"
    raw_prediction_column = "prediction"
    raw_ncbi_gene_id_column = "ncbi_gene_id"
    final_columns = ["Ensembl_ID", "Gene_Name", "miRNA_ID", "miRNA_Name", "Score"]

    log_step(logger, 2, total_steps, "Load raw miRDB predictions")
    pred_df = load_prediction_files(
        raw_predictions_path,
        raw_mirna_column,
        raw_transcript_column,
        raw_prediction_column,
        raw_ncbi_gene_id_column,
    )

    log_step(logger, 3, total_steps, "Create miRNA name-to-MIMAT mapping from miRBase mature.fa")
    mirna_name_to_mimat_map = create_mirna_name_to_mimat_mapping(mirbase_path)

    log_step(logger, 4, total_steps, "Create NCBI Gene ID-to-Ensembl gene mapping from BioMart table")
    biomart_ensembl_id_column = "Gene stable ID"
    biomart_gene_name_column = "Gene name"
    biomart_refseq_column = "RefSeq mRNA ID"
    biomart_ncbi_gene_id_column = "NCBI gene (formerly Entrezgene) ID"
    ncbi_gene_id_to_ensembl_map = create_ncbi_gene_id_to_ensembl_mapping(
        biomart_path,
        biomart_ensembl_id_column,
        biomart_gene_name_column,
        biomart_ncbi_gene_id_column,
    )

    log_step(logger, 5, total_steps, "Create RefSeq-to-Ensembl fallback gene mapping from BioMart table")
    refseq_to_ensembl_map = create_refseq_to_ensembl_mapping(
        biomart_path,
        biomart_ensembl_id_column,
        biomart_gene_name_column,
        biomart_refseq_column,
    )

    log_step(logger, 6, total_steps, "Map prediction miRNA names to MIMAT IDs")
    pred_df = map_mirna_names_to_mimat(
        pred_df,
        mirna_name_to_mimat_map,
        raw_mirna_column,
        "miRNA_Name",
        "miRNA_ID",
    )

    log_step(logger, 7, total_steps, "Map prediction NCBI Gene IDs to Ensembl gene IDs and gene names")
    pred_df = map_ncbi_gene_id_to_ensembl(
        pred_df,
        ncbi_gene_id_to_ensembl_map,
        raw_ncbi_gene_id_column,
        "Ensembl_ID",
        "Gene_Name",
        drop_unmapped=False,
    )

    log_step(logger, 8, total_steps, "Map remaining prediction RefSeq transcript IDs to Ensembl gene IDs and gene names")
    pred_df = fill_unmapped_rows_with_refseq_to_ensembl(
        pred_df,
        refseq_to_ensembl_map,
        raw_transcript_column,
        "Ensembl_ID",
        "Gene_Name",
    )

    log_step(logger, 9, total_steps, "Build, validate, and write standardized output table")
    final_df = build_output_table(
        pred_df,
        raw_prediction_column,
        "Score",
        final_columns,
        "Ensembl_ID",
        "miRNA_ID",
    )
    write_standardized_table(final_df, args.standardized_output_file, logger=logger, columns=final_columns)
    logger.info("Relative output path: %s", resolve_path_relative_to_root(args.standardized_output_file))


if __name__ == "__main__":
    main()
