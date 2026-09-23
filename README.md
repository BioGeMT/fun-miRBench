# fun-miRBenvh

fun-miRBenvh is a benchmarking system for evaluating **functional miRNA target predictors**
against differential-expression (DE) experiments following miRNA perturbation.

The repository contains the benchmark software, curated experiments, predictor registries and
reproducibility pipelines.

## Published data

The curated benchmark inputs are archived in **Zenodo record 21671476**:

https://zenodo.org/records/21671476

The published data include:

- curated miRNA perturbation DE tables;
- standardized outputs for the registered target predictors.

Download the complete published input collection once with:

```bash
uv run funmirbench-download-data
```

The files are cached locally under:

```text
data/experiments/processed/21671476/
data/predictions/<tool>/
```

## Reproduce a benchmark run

### Requirements

- Python 3.10+
- `uv`

Clone the repository and install the locked environment:

```bash
git clone git@github.com:BioGeMT/fun-miRBench.git
cd fun-miRBench
uv sync
```

Download the published benchmark inputs:

```bash
uv run funmirbench-download-data
```

Run the included benchmark configuration:

```bash
uv run funmirbench --config benchmark.yaml
```

`benchmark.yaml` selects experiments and predictors from the versioned metadata registries. The
published data only need to be downloaded once; changing the config changes the subset evaluated.

The included configuration currently selects two curated experiment datasets and all five
registered predictors. Edit the filters in `benchmark.yaml` to evaluate a different subset of the
local published collection.

## Benchmark inputs

### Experiment registry

Experiment metadata are versioned in:

```text
metadata/mirna_experiment_info.tsv
```

Each selected row points to a benchmark-ready DE table. The canonical DE schema is:

| Column | Meaning |
|---|---|
| `gene_id` | Ensembl gene identifier |
| `logFC` | signed log fold change |
| `FDR` | adjusted p-value / q-value used for thresholded evaluation |
| `control_mean_normalized_count` | mean normalized control expression |
| `PValue` | optional raw p-value |

Common native DE column names are normalized when tables are read.

### Predictor registry

Predictor metadata are versioned in:

```text
metadata/predictions_info.tsv
```

Standardized predictor tables use the shared schema:

| Column | Meaning |
|---|---|
| `Ensembl_ID` | Ensembl gene identifier |
| `Gene_Name` | gene symbol |
| `miRNA_ID` | mature miRNA accession |
| `miRNA_Name` | mature miRNA name |
| `Score` | predictor-specific score |

The registry also records `score_direction`, so the benchmark knows whether larger or smaller
values represent stronger predictions.

The published standardized tables use Ensembl release 115 (GRCh38) and miRBase 22.1 annotations.

## Benchmark configuration

The benchmark is driven by one YAML file:

```yaml
experiments_tsv: metadata/mirna_experiment_info.tsv
predictions_tsv: metadata/predictions_info.tsv

experiments:
  id: [GSE169128_KO_miR_3662, GSE124411_OE_miR_150_5p_1]

predictors:
  tool_id: [targetscan, mirdb_mirtarget, microt_cnn, mirbind2, miraw]

evaluation:
  predictor_load_workers: 2
  fdr_threshold: 0.05
  effect_threshold: 1.0
  predictor_top_fraction: 0.10
  write_top_prediction_cdfs: true
  report_min_common_coverage: 0.10
  protein_coding_only: true
  protein_coding_gtf: data/common_resources/ensembl/Homo_sapiens.GRCh38.115.gtf.gz
  protein_coding_gene_cache: data/common_resources/ensembl/protein_coding_gene_ids.txt
  figure_dpi: 600

out_dir: results/
```

Experiment and predictor filters use **AND across columns** and **OR within a column's value list**.

Paths may be absolute or relative to the YAML file.

## Evaluation design

By default, evaluation is restricted to Ensembl release 115 protein-coding genes. The annotation
resource is downloaded or reused on the first protein-coding run and cached under
`data/common_resources/`.

For each experiment:

- predictor scores are joined to the DE table by miRNA-gene pair;
- missing predictor pairs remain missing and are not converted to zero-valued predictions;
- coverage is reported explicitly for each predictor;
- ground-truth positives are derived from the configured FDR and effect-size thresholds;
- predictor-specific score direction is respected;
- per-dataset plots use dataset-local tie-aware ranks;
- cross-dataset rank analyses use ranks derived from the full standardized predictor table;
- exact top-k selection with deterministic tie breaking is used for predictor-agreement analyses;
- combined PR, ROC, and GSEA comparisons use the common scored set for the predictors being
  compared.

This keeps predictive performance separate from prediction coverage and makes the evaluated gene
universe explicit.

## Outputs

Each run creates a date-based directory:

```text
results/YYYYMMDD_HHMMSS/
```

A run contains:

```text
results/YYYYMMDD_HHMMSS/
├── benchmark_config.yaml
├── README.md
├── REPORT.pdf
├── summary.json
├── datasets/
│   └── <dataset_id>/
│       ├── joined.tsv
│       ├── plots/
│       └── reports/
├── tables/
│   ├── per_experiment/
│   └── combined/
└── plots/
    └── combined/
```

The config snapshot, joined tables, metric tables, coverage summaries, reports, and figures make
each benchmark run self-describing and auditable.

## Using your own data

The Zenodo dataset is optional.

To benchmark custom experiments or predictors, create metadata tables with the same contracts and
point a YAML configuration to them:

```yaml
experiments_tsv: my_experiments.tsv
predictions_tsv: my_predictions.tsv
```

Then run:

```bash
uv run funmirbench --config my_benchmark.yaml
```

In this workflow, `funmirbench-download-data` is not required. The benchmark reads the files
referenced by your metadata and does not contact Zenodo for experiment or predictor inputs.

## Reproducing the published inputs

The published files can be used directly for benchmarking, but the repository also retains the
processing pipelines used to construct standardized inputs.

### Experiment processing

`pipelines/experiments/` converts a count matrix or FASTQ reads into a canonical DE table. Reads
mode runs FastQC, fastp, STAR, featureCounts, and DESeq2; count-matrix mode runs DESeq2 directly.

For GEO/SRA acquisition and config generation, see:

- [GEO download pipeline](pipelines/geo/README.md)
- [Experiment ingestion pipeline](pipelines/experiments/README.md)

These workflows require the additional Conda environment documented in the experiment pipeline.

### Predictor standardization

`pipelines/standardized_predictors/` contains the reproducibility pipelines for:

- TargetScan v8
- miRDB
- microT-CNN
- miRBind2-3UTR
- miRAW

These pipelines regenerate the common five-column predictor tables from the respective upstream
sources. They are not required when using the published standardized files from Zenodo.

See [standardized predictor pipelines](pipelines/standardized_predictors/README.md).

## Manuscript assets

Scripts used to generate manuscript-supporting figures and compact tables are under `scripts/`.
Final manuscript-facing assets are stored under `manuscript_assets/`.

See [manuscript figure scripts](scripts/README.md) for the exact regeneration commands.

## Repository structure

```text
funmirbench/                     benchmark package
metadata/                        versioned experiment and predictor registries
data/                            local downloaded/generated data cache
pipelines/geo/                   GEO/SRA acquisition and config generation
pipelines/experiments/           DE processing workflow
pipelines/standardized_predictors/
                                 predictor standardization workflows
scripts/                         manuscript-supporting analyses
manuscript_assets/               retained manuscript-facing outputs
results/                         local benchmark runs
benchmark.yaml                   included benchmark configuration
```

Large downloaded data, intermediate pipeline files, and benchmark run directories are intentionally
not versioned in Git.

## Command reference

```bash
# Download the complete published benchmark input collection once
uv run funmirbench-download-data

# Run a benchmark configuration
uv run funmirbench --config benchmark.yaml

# Validate registered experiment tables
uv run funmirbench-validate-experiments \
  --experiments-tsv metadata/mirna_experiment_info.tsv

# Download example inputs for the experiment-processing workflow
uv run funmirbench-experiments-download-examples

# Run the experiment-processing workflow
uv run funmirbench-experiments --config <experiment-config.yaml>

# Synchronize generated experiment metadata into the registry
uv run funmirbench-sync-metadata
```
