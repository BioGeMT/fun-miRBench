# microT-CNN

This directory contains the standardization pipeline for microT-CNN gene-level predictions.

## Files

- `pipeline.py`: CLI entrypoint for the pipeline.
- `utils.py`: helpers for logging, downloads, mapping, score conversion, and output construction.
- `data/predictions/microt_cnn/microt_cnn_pipeline.log`: example log from a completed run, relative to the repository root.

## What The Pipeline Does

The pipeline:

1. Downloads miRBase `mature.fa` version 22.1.
2. Downloads the microT-CNN all-score predictions file (via Zenodo record by default).
3. Loads the raw prediction columns used by the standardized output.
4. Maps miRNA human names (`hsa-*`) to `MIMAT` IDs using `mature.fa`.
5. Maps raw Ensembl transcript IDs (`ENST*`) to Ensembl gene IDs (`ENSG*`) using an Ensembl v115 transcript-to-gene cache. If the cache is missing, the pipeline downloads the Ensembl v115 GTF and builds it.
6. Drops rows that fail miRNA name -> MIMAT or transcript -> gene mapping.
7. Converts the prediction column to numeric `Score`.
8. Drops rows with non-numeric `Score` values after numeric conversion.
9. Collapses transcript rows to one row per miRNA-gene pair. In the current source file this is a no-op because each mapped miRNA-gene pair appears once; if duplicates appear, the strongest score is kept and the duplicate count is logged.
10. Builds the final output table and writes the standardized TSV to the default location.

## Output Schema

The output TSV contains:

- `Ensembl_ID`
- `Gene_Name`
- `miRNA_ID`
- `miRNA_Name`
- `Score`

`Ensembl_ID` is the mapped Ensembl gene ID, not the raw transcript ID. `miRNA_Name` is copied from the raw miRNA name column. `Score` is the numeric form of the raw prediction value. Rows with non-numeric scores are dropped.

## Output Location

By default, the standardized file is written to:

```
data/predictions/microt_cnn/microt_cnn_standardized.tsv
```

relative to the repository root.

## Resource Cache

The pipeline downloads external resources only when the expected cache files are missing. By default, the cache files are:

```
data/common_resources/mirbase/mature.fa
data/common_resources/ensembl/Homo_sapiens.GRCh38.115.gtf.gz
pipelines/standardized_predictors/microt_cnn/data/microT_CNN_prediction_result_human_all_scores_gene_level.tsv.gz
pipelines/standardized_predictors/microt_cnn/resources/ensembl115_tx2gene.tsv.gz
```

Local prediction-file overrides are not supported. The Zenodo source and cache path above are the canonical source for reproducible standardization.

The log notes whether each resource was reused from cache or downloaded.

## Run

From the repository root:

```bash
uv run pipelines/standardized_predictors/microt_cnn/pipeline.py
```

## CLI Arguments

```bash
uv run pipelines/standardized_predictors/microt_cnn/pipeline.py \
  --common-resources-dir data/common_resources \
  --data-dir pipelines/standardized_predictors/microt_cnn/data \
  --resources-dir pipelines/standardized_predictors/microt_cnn/resources \
  --standardized-output-file data/predictions/microt_cnn/microt_cnn_standardized.tsv \
  --log-file data/predictions/microt_cnn/microt_cnn_pipeline.log \
  --log-level INFO
```

Relative CLI paths are resolved from the repository root.

## Logging

Logging is written to stdout and to the file passed via `--log-file`. By default, the log is written to `data/predictions/microt_cnn/microt_cnn_pipeline.log`. Main processing stages are logged as numbered steps and row-count changes are logged in a `before -> after` format.
