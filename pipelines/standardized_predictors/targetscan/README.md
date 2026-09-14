# TargetScan

This directory contains the standardization pipeline for TargetScan v8 gene-level predictions.

## Files

- `pipeline.py`: CLI entrypoint for the pipeline.
- `utils.py`: shared helpers for logging, downloads, parsing, QC, mapping, and output construction.
- `data/predictions/targetscan/targetscan_pipeline.log`: example log from a completed run, relative to the repository root.
- `data/`: downloaded TargetScan prediction data.
- `resources/`: downloaded TargetScan metadata and derived Ensembl lookup tables.

Shared miRBase and Ensembl annotations are downloaded to:

```text
data/common_resources/mirbase/mature.fa
data/common_resources/ensembl/Homo_sapiens.GRCh38.115.gtf.gz
```

By default, `Summary_Counts.all_predictions.txt` and its ZIP archive are stored under `data/`.
`Gene_info.txt`, `miR_Family_Info.txt`, their ZIP archives, and derived Ensembl lookup tables are
stored under `resources/`.

## What The Pipeline Does

The pipeline:

1. Downloads or reuses the TargetScan v8 input files:
   - `Summary_Counts.all_predictions.txt`
   - `miR_Family_Info.txt`
   - `Gene_info.txt`
2. Builds a representative-transcript index from `Gene_info.txt` for human rows only (`Species ID = 9606`).
3. Keeps every transcript with `Representative transcript? == 1`.
4. Logs whether any gene has more than one representative transcript row.
5. Downloads or reuses the Ensembl v115 GTF for GRCh38.
6. Builds cached Ensembl lookup tables for:
   - transcript stable ID to gene stable ID
   - gene stable ID to gene name
7. Runs QC comparing TargetScan representative transcripts against Ensembl transcript and gene mappings.
8. Downloads or reuses miRBase `mature.fa` version 22.1.
9. Builds a human miRNA mapping from TargetScan `MiRBase ID` to:
   - `miRNA_ID` as miRBase accession
   - `miRNA_Name` as the miRBase mature name
10. Reads `Summary_Counts.all_predictions.txt` and keeps only human rows.
11. Drops rows where:
   - the transcript is not in the representative transcript set
   - the transcript does not map to an Ensembl gene
   - the TargetScan gene ID for that transcript disagrees with the Ensembl gene ID recovered from the transcript
   - the representative miRNA cannot be mapped through the miRNA annotation table
   - the score is null
12. Preserves the raw TargetScan `Cumulative weighted context++ score` as `Score`.
13. Checks whether any duplicate `(Ensembl_ID, miRNA_ID)` pairs remain after all mapping and mismatch filters.
14. Raises an error instead of silently choosing between duplicate scores if any such pairs survive.
15. Writes the standardized output table only when the filtered rows are already unique at the final gene-miRNA level.

## Output Schema

The output TSV contains:

- `Ensembl_ID`
- `Gene_Name`
- `miRNA_ID`
- `miRNA_Name`
- `Score`

`Score` is the raw numeric form of TargetScan `Cumulative weighted context++ score`.
No score-based duplicate collapse is applied in the final output step. If duplicate gene-miRNA pairs survive filtering, the pipeline fails and reports them.

## Output Location

By default, the standardized file is written to:

```text
data/predictions/targetscan/targetscan_standardized.tsv
```

relative to the repository root.

## Run

From the repository root:

```bash
uv run pipelines/standardized_predictors/targetscan/pipeline.py
```

## CLI Arguments

```bash
uv run pipelines/standardized_predictors/targetscan/pipeline.py \
  --common-resources-dir data/common_resources \
  --data-dir pipelines/standardized_predictors/targetscan/data \
  --resources-dir pipelines/standardized_predictors/targetscan/resources \
  --standardized-output-file data/predictions/targetscan/targetscan_standardized.tsv \
  --log-file data/predictions/targetscan/targetscan_pipeline.log \
  --log-level INFO
```

## Logging

Logging is written both to stdout and to the file passed via `--log-file`. By default, the log is written to `data/predictions/targetscan/targetscan_pipeline.log`.
The duplicate summary distinguishes provisional duplicate gene-miRNA pairs seen before the TargetScan/Ensembl mismatch filter from any duplicate pairs that survive after filtering.
