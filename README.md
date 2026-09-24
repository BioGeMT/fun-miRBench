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

Each selected row has a stable `id` for filenames and machine-readable outputs, plus explicit
human-facing identifiers such as `geo_accession`, `mirna_name`, `tested_cell_line`, and
`experiment_type`. Reports and plots combine those fields into readable experiment labels while
retaining the stable `id` for provenance.

The registry uses the following identity fields:

| Column | Meaning |
|---|---|
| `id` | stable experiment identifier used in filenames, directories, configs, and machine-readable outputs |
| `geo_accession` | GEO series accession, for example `GSE115646` |
| `mirna_name` | mature miRNA name |
| `experiment_type` | perturbation type, such as Overexpression, Knockout, or Knockdown |
| `tested_cell_line` | experimental cell line or cell context when available |
| `tissue` | tissue context when available |
| `organism` | source organism |
| `method` | expression profiling method |
| `pubmed_id` | PubMed identifier when available |
| `de_table_path` | local benchmark-ready DE table path |

The default reporting label is `GEO accession · miRNA · cell line · experiment type`. Missing
optional fields are omitted from the label rather than shown as `NA`.

Each row points to a benchmark-ready DE table. The canonical DE schema is:

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

The benchmark is controlled by a single YAML file. The included `benchmark.yaml` is a runnable
configuration for the published data:

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

### Configuration reference

| Key | Meaning |
|---|---|
| `experiments_tsv` | Experiment registry to use. Each selected row must point to a benchmark-ready DE table. |
| `predictions_tsv` | Predictor registry to use. Each selected row defines a standardized predictor file and its score direction. |
| `experiments` | Optional filters applied to the experiment registry. Omit this block to use all registered experiments. |
| `predictors` | Optional filters applied to the predictor registry. Omit this block to use all registered predictors. |
| `evaluation` | Evaluation and reporting settings described below. |
| `out_dir` | Root directory for benchmark outputs. Each run creates a timestamped subdirectory such as `results/20260923_120807/`. |

Filters can target any column in the corresponding metadata table. Different filter keys are
combined with **AND**; multiple values for one key are combined with **OR**. Paths may be absolute
or relative to the YAML file.

### Evaluation settings

| Setting | Example | Meaning |
|---|---:|---|
| `predictor_load_workers` | `2` | Number of workers used to load standardized predictor tables in parallel. This affects loading speed, not benchmark results. |
| `fdr_threshold` | `0.05` | FDR cutoff used when defining ground-truth positive genes. A positive must satisfy `FDR < threshold`. Set to `null` to define positives from effect size only. |
| `effect_threshold` | `1.0` | Minimum perturbation-aware effect used for ground-truth positives. The effect is `-logFC` for overexpression experiments and `+logFC` for knockout/knockdown experiments. |
| `predictor_top_fraction` | `0.10` | Fraction of each predictor's scored genes treated as its top predictions in agreement analyses. Selection uses an exact top-k with deterministic tie breaking. |
| `write_top_prediction_cdfs` | `true` | Whether to write the optional top-prediction effect CDF diagnostic plots. |
| `report_min_common_coverage` | `0.10` | Minimum fraction of experiment rows a predictor must score to be included in report-level common-set comparison plots. |
| `protein_coding_only` | `true` | Restrict the evaluation universe to Ensembl protein-coding genes. |
| `protein_coding_gtf` | Ensembl v115 GTF | GTF used to derive the protein-coding gene set. If the configured file is absent, the default Ensembl resource is downloaded. |
| `protein_coding_gene_cache` | `.../protein_coding_gene_ids.txt` | Local cache of protein-coding Ensembl gene IDs derived from the GTF. |
| `figure_dpi` | `600` | Resolution used for raster benchmark figures. |

## Evaluation behavior

The configuration above defines the evaluation choices that a user may change. The remaining
benchmark mechanics are fixed by the software:

- predictor scores are joined to each DE table by miRNA-gene pair, and missing predictions remain
  missing rather than being replaced with zero;
- coverage is reported separately from predictive performance, so a predictor is not penalized by
  silently treating unscored pairs as negative predictions;
- predictor score direction is normalized internally so that stronger predictions are interpreted
  consistently across tools;
- multi-predictor PR, ROC, and GSEA comparisons are calculated on the common set of genes scored by
  all predictors being compared.

The generated run `README.md`, `REPORT.pdf`, and `summary.json` record the settings and
evaluation details for that specific run.

## Outputs

Each run creates a timestamped directory:

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
