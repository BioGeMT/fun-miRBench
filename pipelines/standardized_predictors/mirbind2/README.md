# miRBind2

This directory contains the standardization pipeline for miRBind2 human 3UTR predictions.

## Files

- `pipeline.py`: CLI entrypoint for the pipeline.
- `utils.py`: helpers for logging, cleaning, mapping, score conversion, and output construction.
- `data/predictions/mirbind2/mirbind2_pipeline.log`: example log from a completed run, relative to the repository root.

## Inputs

The raw input file is not tracked in Git because it is large. By default, the pipeline downloads it from Zenodo and caches it in this directory with this exact filename:

```text
pipelines/standardized_predictors/mirbind2/data/3utrs_mirbind2_predictions.tsv.gz
```

Source URL:

```text
https://zenodo.org/records/20609975/files/3utrs_mirbind2_predictions.tsv.gz?download=1
```

Local prediction-file overrides are not supported. The Zenodo file and cache path above are the canonical source for reproducible standardization.

Its header is:

```text
Transcript_ID	Gene_ID	Gene_Symbol	miRNA_Name	miRBind2_3UTR_prediction
```

The pipeline reads the raw miRBind2 score from `miRBind2_3UTR_prediction`, but this raw column name is not kept in the standardized output. The pipeline validates these raw columns:

- `Transcript_ID`
- `Gene_ID`
- `Gene_Symbol`
- `miRNA_Name`
- `miRBind2_3UTR_prediction`

The pipeline uses the existing shared miRBase resource:

```text
data/common_resources/mirbase/mature.fa
```

This is miRBase release 22.1. The pipeline downloads it when missing and reuses the shared cache
on later runs. The raw miRBind2 table already contains Ensembl gene IDs and gene symbols, so no
BioMart remapping is performed.

## What The Pipeline Does

The pipeline:

1. Loads the raw miRBind2 TSV.
2. Validates the raw miRBind2 columns listed above.
3. Drops rows with invalid values in `Gene_ID`, `miRNA_Name`, or `miRBind2_3UTR_prediction`.
4. Deduplicates exact raw rows on those same columns.
5. Builds a human miRNA name to `MIMAT` ID mapping from miRBase 22.1 `mature.fa`.
6. Annotates `miRNA_Name` to `miRNA_ID`.
7. Drops rows whose miRNA names cannot be mapped to `MIMAT` IDs.
8. Converts the raw prediction column to numeric `Score`.
9. Builds the final output table with the shared standardized schema.
10. Writes the standardized TSV to the default location.

## Output Schema

The output TSV contains:

- `Ensembl_ID`
- `Gene_Name`
- `miRNA_ID`
- `miRNA_Name`
- `Score`

`Gene_Name` is copied from the raw `Gene_Symbol` column. `miRNA_Name` is copied from the raw miRNA name column. `Score` is the numeric form of the raw miRBind2 prediction value. The raw `Transcript_ID` and `miRBind2_3UTR_prediction` columns are not retained in the standardized output.

The score is a signed raw miRBind2-3UTR repression prediction, not a probability.
Lower, more negative scores are treated as stronger predicted repression.

## Output Location

By default, the standardized file is written to:

```text
data/predictions/mirbind2/mirbind2_standardized.tsv
```

relative to the repository root.

## Run

From the repository root:

```bash
uv run pipelines/standardized_predictors/mirbind2/pipeline.py
```

## CLI Arguments

Relative CLI paths are resolved from the repository root.

```bash
uv run pipelines/standardized_predictors/mirbind2/pipeline.py \
  --common-resources-dir data/common_resources \
  --data-dir pipelines/standardized_predictors/mirbind2/data \
  --standardized-output-file data/predictions/mirbind2/mirbind2_standardized.tsv \
  --log-file data/predictions/mirbind2/mirbind2_pipeline.log \
  --log-level INFO
```

## Logging

Logging is written both to stdout and to the log file passed via `--log-file`. By default, the log is written to `data/predictions/mirbind2/mirbind2_pipeline.log`. Main processing stages are logged as numbered steps and row-count changes are logged in a `before -> after` format.
