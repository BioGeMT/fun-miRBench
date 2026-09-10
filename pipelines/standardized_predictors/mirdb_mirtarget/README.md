# miRDB MirTarget v4.0

This directory contains the standardization pipeline for miRDB gene-level predictions.

## Files

- `pipeline.py`: CLI entrypoint for the pipeline.
- `utils.py`: shared helpers for logging, downloads, cleaning, mapping, and output construction.
- `data/predictions/mirdb_mirtarget/mirdb_mirtarget_pipeline.log`: log file written by the default run, relative to the repository root.

The pipeline downloads and reuses the raw miRDB predictions file at:

```text
pipelines/standardized_predictors/mirdb_mirtarget/data/miRDB_v6.0_prediction_result_human_all_scores.txt.gz
```

This file is downloaded from:

```text
https://mirdb.org/download/miRDB_v6.0_prediction_result_human_all_scores.txt.gz
```

If this cache file already exists, the pipeline reuses it. Otherwise, it downloads the miRDB v6.0 all-score prediction file and logs whether the resource was reused, downloaded, or failed.

Local prediction-file overrides are not supported. The miRDB download and cache path above are the canonical source for reproducible standardization.

The pipeline downloads and reuses miRBase from the shared annotation cache:

```text
data/common_resources/mirbase/mature.fa
```

The predictor-specific BioMart mapping remains under:

```text
pipelines/standardized_predictors/mirdb_mirtarget/resources/hsapiens_ncbi_gene_id_refseq_to_ensembl.tsv
```

If either cache file already exists, the pipeline reuses it. Otherwise, it downloads miRBase
`mature.fa` version 22.1 and the BioMart NCBI Gene ID/RefSeq-to-Ensembl mapping table, logging
whether each resource was reused, downloaded, or failed.

The raw miRDB file is gzip-compressed and is treated as a 4-column tab-separated table without a header:

1. `miRNA`
2. `refseq_id`
3. `prediction`
4. `ncbi_gene_id`

The pipeline reads and names all four columns, but gene-level mapping is attempted in this order:

1. `ncbi_gene_id`
2. `refseq_id` fallback for rows still unmapped after step 1

## What The Pipeline Does

The pipeline:

1. Downloads the miRDB v6.0 all-score raw predictions file.
2. Downloads miRBase `mature.fa` version 22.1.
3. Downloads a BioMart TSV with headers and unique rows containing:
   - `Gene stable ID`
   - `Gene name`
   - `RefSeq mRNA ID`
   - `NCBI gene (formerly Entrezgene) ID`
4. Loads the raw miRDB predictions file.
5. Drops rows with missing or invalid values in:
   - `miRNA`
   - `refseq_id`
   - `prediction`
   - `ncbi_gene_id`
6. Deduplicates exact duplicate raw rows on those same four columns.
7. Raises if the same `(miRNA, refseq_id)` pair has conflicting scores in the raw input.
8. Raises if the same `(miRNA, ncbi_gene_id)` pair has conflicting scores in the raw input.
9. Builds a miRNA mapping from human miRBase names (`hsa-*`) to `MIMAT` IDs.
10. Builds an NCBI Gene ID to `(Ensembl_ID, Gene_Name)` mapping from BioMart.
11. Builds a RefSeq mRNA to `(Ensembl_ID, Gene_Name)` mapping from the same BioMart file.
12. Maps:
    - `miRNA` to `miRNA_ID`
    - `ncbi_gene_id` to `Ensembl_ID` and `Gene_Name`
13. For rows still missing an Ensembl mapping, falls back to `refseq_id`.
14. Drops rows that still fail gene mapping after the fallback.
15. Converts `prediction` to numeric `Score`.
16. Drops and reports final `(Ensembl_ID, miRNA_ID)` pairs with conflicting scores after gene mapping.
17. Drops only exact duplicate final rows on:
    - `Ensembl_ID`
    - `miRNA_ID`
    - `Score`
18. Writes the standardized output table.

## Output Schema

The output TSV contains:

- `Ensembl_ID`
- `Gene_Name`
- `miRNA_ID`
- `miRNA_Name`
- `Score`

`miRNA_Name` is copied from the raw `miRNA` column. `Score` is the numeric form of `prediction`.

## Output Location

By default, the standardized file is written to:

```text
data/predictions/mirdb_mirtarget/mirdb_mirtarget_standardized.tsv
```

relative to the repository root.

## Run

From the repository root:

```bash
uv run pipelines/standardized_predictors/mirdb_mirtarget/pipeline.py
```

## CLI Arguments

Relative CLI paths are resolved from the repository root, so the current working directory does not change where inputs, resources, logs, or outputs are read or written.

```bash
uv run pipelines/standardized_predictors/mirdb_mirtarget/pipeline.py \
  --common-resources-dir data/common_resources \
  --data-dir pipelines/standardized_predictors/mirdb_mirtarget/data \
  --resources-dir pipelines/standardized_predictors/mirdb_mirtarget/resources \
  --standardized-output-file data/predictions/mirdb_mirtarget/mirdb_mirtarget_standardized.tsv \
  --log-file data/predictions/mirdb_mirtarget/mirdb_mirtarget_pipeline.log \
  --log-level INFO
```

## Logging

Logging is written both to stdout and to the log file passed via `--log-file`. By default, the log is written to `data/predictions/mirdb_mirtarget/mirdb_mirtarget_pipeline.log`.
