# Standardized Predictors

This directory contains pipelines that generate predictor outputs in a common standardized schema for downstream benchmarking.

Current predictor pipelines:

- `targetscan/`
- `microt_cnn/`
- `mirbind2/`
- `mirdb_mirtarget/`
- `miraw/`

## Standardized Schema

The predictor outputs are written in a shared TSV format with the columns:

- `Ensembl_ID`
- `Gene_Name`
- `miRNA_ID`
- `miRNA_Name`
- `Score`

The shared annotation schema uses Ensembl v115 (GRCh38) and miRBase release 22.1. Each standardized predictor table is expected to populate all five columns.

Shared downloaded annotations are cached under `data/common_resources/` at the repository root.
All predictor pipelines reuse `mirbase/mature.fa`; microT-CNN, miRAW, and TargetScan also reuse
`ensembl/Homo_sapiens.GRCh38.115.gtf.gz`.

Every pipeline accepts `--common-resources-dir`, `--data-dir`,
`--standardized-output-file`, `--log-file`, and `--log-level`. Pipelines with supporting or
derived predictor-specific resources also accept `--resources-dir`. Relative paths are resolved
from the repository root.

## Pipelines

- `targetscan/`
- `microt_cnn/`
- `mirbind2/`
- `mirdb_mirtarget/`
- `miraw/`

See the README in each pipeline directory for pipeline-specific inputs, processing steps, and outputs.
