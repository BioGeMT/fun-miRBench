# FuNmiRBench

Benchmark functional miRNA target predictors against differential-expression tables.

## Install

Requirements:

- Python 3.10+
- `uv`
- `conda` for the experiments pipeline environment

Install `uv` on your machine first:

```bash
python -m pip install uv
```

Then clone the repo and install the Python package environment:

```bash
git clone git@github.com:BioGeMT/FuNmiRBench.git
cd FuNmiRBench
uv sync
```

If you want to use the experiments pipeline, create and activate the extra local environment after entering the repo:

```bash
conda env create -f pipelines/experiments/environment.yml
conda activate funmirbench-experiments
```

That environment also includes `uv`, so `uv run ...` keeps working after activation.

## Repo Layout

Main directories:

- `data/experiments/processed/`: root directory for processed experiment DE tables
- `data/experiments/processed/21671476/`: local cache for curated benchmark DE tables from Zenodo record `21671476`
- `data/experiments/raw/`: local raw GEO inputs such as count matrices and FASTQs
- `data/common_resources/`: shared, downloaded miRBase and Ensembl annotation cache
- `data/predictions/`: local cache for published standardized predictor TSVs; predictor pipelines can also regenerate these files
- `metadata/mirna_experiment_info.tsv`: experiment registry
- `metadata/predictions_info.tsv`: predictor registry
- `pipelines/experiments/`: experiment-ingestion backend files and example configs
- `pipelines/standardized_predictors/`: reproducibility pipelines used to regenerate standardized predictor tables from upstream sources
- `manuscript_assets/`: final manuscript figure exports and compact supporting tables
- `results/`: benchmark outputs

The benchmark reads file paths from the two metadata TSVs. `data/` holds the real files. `results/`
is only for benchmark output.

## Quick Start

If you just want to run the benchmark, you do not need the experiments pipeline first. The repo
already includes:

- experiment metadata in `metadata/mirna_experiment_info.tsv`
- predictor metadata in `metadata/predictions_info.tsv`

To reproduce the published FuNmiRBench benchmark, download the published benchmark inputs once:

```bash
uv run funmirbench-download-data
```

This downloads all curated experiment DE tables and all published standardized predictor tables
from Zenodo record `21671476` into the local `data/` cache.

Then run the benchmark:

```bash
uv run funmirbench --config benchmark.yaml
```

The benchmark command itself does not download experiment or predictor data. The config only
selects and filters the local files registered in the metadata tables.

The default config already points at:

- `metadata/mirna_experiment_info.tsv`
- `metadata/predictions_info.tsv`

and selects a small set of curated experiment datasets plus the registered external predictors.

## Workflow

### 1. Add Experiment Data

Experiment DE tables are produced by a dedicated ingestion pipeline that runs DESeq2 on a count
matrix or on FASTQs (via FastQC + fastp + STAR + featureCounts), then syncs the result into
`metadata/mirna_experiment_info.tsv`. Raw data can also be fetched directly from GEO and turned
into a ready-to-review YAML config automatically. See the dedicated pipelines for the full
workflow:

- [GEO download pipeline](pipelines/geo/README.md) — fetch FASTQs/count matrices from GEO and
  auto-generate a YAML config
- [Experiment ingestion pipeline](pipelines/experiments/README.md) — run DESeq2 (count matrix or
  reads mode) and sync results into the registry

For the curated benchmark datasets tracked in `metadata/mirna_experiment_info.tsv`, the metadata
rows stay versioned in the repo while the corresponding DE tables are published in Zenodo record
`21671476`. Running `uv run funmirbench-download-data` downloads the full curated experiment
collection into `data/experiments/processed/21671476/`. The benchmark then only filters and reads
those local tables according to the selected config.

### 2. Add Predictor Data

Predictor score files live under `data/predictions/` and are discovered through
`metadata/predictions_info.tsv`. The published standardized predictor tables are part of Zenodo
record `21671476`. Running `uv run funmirbench-download-data` downloads the published predictor
archive, verifies its Zenodo checksum, and prepares all registered standardized predictor tables
under `data/predictions/<tool>/`.

The pipelines under `pipelines/standardized_predictors/` are the reproducibility path for
regenerating those published standardized tables from their upstream sources. They are not a
prerequisite for using the published Zenodo artifacts.

### 3. Run The Benchmark

The default benchmark config is `benchmark.yaml`. The primary benchmark universe is
Ensembl protein-coding genes only; FuNmiRBench filters each DE table to that universe
before joining predictor scores.

Benchmark config summary:

- `experiments_tsv`: experiment metadata table
- `predictions_tsv`: predictor metadata table
- `experiments`: which experiment rows to include
- `predictors`: which predictor rows to include
- `evaluation`: thresholds, ranking settings, and the protein-coding benchmark universe
- `out_dir`: results root directory; each benchmark run creates its own subfolder under this root

Run it with:

```bash
uv run funmirbench --config benchmark.yaml
```

Before the first run with the published FuNmiRBench inputs, prepare the complete local dataset once:

```bash
uv run funmirbench-download-data
```

After that, `benchmark.yaml` only controls which experiments and predictors are selected from the
local metadata registries. On the first protein-coding benchmark run, FuNmiRBench may still
download or reuse the Ensembl v115 annotation resource and cache the protein-coding gene set at
`data/common_resources/ensembl/protein_coding_gene_ids.txt`.

During evaluation, DE rows are first restricted to Ensembl protein-coding genes from Ensembl
release 115. This keeps the evaluation universe aligned with mRNA target predictors and avoids
penalizing tools for non-coding genes they are not designed to score. Each predictor is then
scored only on miRNA-gene pairs that exist in that predictor's standardized file. Missing pairs
are not filled with zero for metrics. Each run writes coverage information to
`tables/per_experiment/coverage_per_experiment.tsv`, and the per-predictor Markdown/PDF reports
also record total rows, scored rows, missing rows, and coverage. For per-dataset heatmaps and agreement plots,
FuNmiRBench uses a dataset-local tie-aware rank over the scored rows. For cross-dataset
rank-distribution plots, it keeps a separate global tie-aware rank derived from each predictor's
full standardized file. Predictor-agreement top fractions use an exact top-k selection per
predictor with a deterministic tie-break instead of a quantile threshold. Combined PR, ROC, and
GSEA comparison plots are computed on the common set of genes scored by all compared predictors.

YAML paths can be:

- absolute paths
- relative to the YAML file
- repo-root-relative paths such as `data/...` and `metadata/...`

The default config shape is:

```yaml
experiments_tsv: metadata/mirna_experiment_info.tsv
predictions_tsv: metadata/predictions_info.tsv

experiments:
  id: [GSE169128_KO_miR_3662, GSE124411_OE_miR_150_5p_1]

predictors:
  tool_id: [targetscan, mirdb_mirtarget, microt_cnn, mirbind2, miraw]

evaluation:
  fdr_threshold: 0.05
  effect_threshold: 1.0
  predictor_top_fraction: 0.10
  write_top_prediction_cdfs: true
  report_min_common_coverage: 0.10
  protein_coding_only: true
  protein_coding_gtf: data/common_resources/ensembl/Homo_sapiens.GRCh38.115.gtf.gz
  protein_coding_gene_cache: data/common_resources/ensembl/protein_coding_gene_ids.txt

out_dir: results/
```

Filter behavior:

- different keys are combined with AND
- values inside one list are combined with OR

### Using Your Own Data

Zenodo is optional. To benchmark your own experiments or predictor outputs, provide your own
metadata TSV files and point a config at them:

```yaml
experiments_tsv: my_experiments.tsv
predictions_tsv: my_predictions.tsv
```

Then run:

```bash
uv run funmirbench --config my_benchmark.yaml
```

In this workflow, `funmirbench-download-data` is not needed and the benchmark does not contact
Zenodo for experiment or predictor inputs.

`benchmark.yaml` already includes other experiment and predictor columns as commented rows, so the
normal workflow is just to edit or uncomment filters.

## Manuscript Assets

Manuscript figure scripts live under `scripts/`; see `scripts/README.md` for
the exact commands. The repo tracks final manuscript-facing assets rather than
all generated intermediates. For Figure 2, this means tracking the combined
figure exports and compact TSVs that support the manuscript values. Individual
panel image exports, large raw support tables such as gene-level conservation,
downloaded resources, diagnostic outputs, and benchmark run folders are local
generated artifacts and remain gitignored.

## Outputs

After a benchmark run, `results/` contains one new run folder named with the run date. For example:

- `results/20260510/`

If another run is created on the same date, FuNmiRBench keeps the date-based naming and adds a
collision suffix such as `results/20260510__r2/`.

Detailed dataset, miRNA, predictor, perturbation, cell-line, evaluation-threshold, and
protein-coding-universe metadata is recorded inside the run folder.

Inside each run folder you get:

- `README.md`: human-readable run guide and quick-start map for the output folder
- `REPORT.pdf`: main run-level PDF report with explanations and selected combined plots
- `datasets/<dataset_id>/joined.tsv`: joined DE + predictor score table for that dataset
- `datasets/<dataset_id>/plots/predictors/<tool_id>/`: per-tool plots for that dataset
- `datasets/<dataset_id>/plots/comparisons/`: multi-predictor comparison plots for that dataset
- `datasets/<dataset_id>/plots/heatmaps/`: dataset-level heatmaps
- `datasets/<dataset_id>/reports/`: per-dataset Markdown/PDF reports and correlation TSVs
- `tables/per_experiment/`: per-experiment metric tables
- `tables/combined/`: cross-dataset predictor summary table
- `plots/combined/metrics/`, `plots/combined/ranks/`, `plots/combined/combinations/`: cross-dataset comparison plots grouped by theme
- `summary.json`: run summary

When 2 or more predictors are selected, each dataset gets:

- one score-vs-expected-effect scatter per predictor
- one combined PR curve on common scored pairs
- one combined ROC curve on common scored pairs
- one algorithms-vs-genes heatmap
- common-prediction overlap summaries for comparable predictor sets

## Commands

```bash
uv run funmirbench-download-data
uv run funmirbench --config benchmark.yaml
uv run funmirbench-validate-experiments --experiments-tsv metadata/mirna_experiment_info.tsv
uv run funmirbench-experiments-download-examples
uv run funmirbench-experiments --config config.yaml
uv run funmirbench-sync-metadata
```

## Tests

```bash
uv run python -m unittest discover -s tests
```
