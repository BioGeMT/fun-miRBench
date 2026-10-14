# Standardized Predictor Pipelines

This directory contains the **reproduction workflows** used to generate the standardized predictor
tables evaluated by fun-miRBench.

These pipelines are not required to reproduce the published benchmark. The published standardized
predictor tables can be downloaded directly from the project Zenodo release with:

```bash
uv run funmirbench-download-data
```

Run the pipelines in this directory when you want to regenerate a published standardized table
from its upstream source, audit the standardization procedure, or add a new predictor.

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

## Reproducibility contract

Each pipeline writes the same five-column output contract while preserving the predictor-specific
score itself. Predictor directionality is recorded separately in
`metadata/predictions_info.tsv`.

Generated standardized tables are written under `data/predictions/<tool>/`, matching the paths
used by the published metadata registry.

## Pipelines

- `targetscan/`
- `microt_cnn/`
- `mirbind2/`
- `mirdb_mirtarget/`
- `miraw/`

See the README in each predictor directory for its source data, processing decisions, command-line
options, and output definition.
