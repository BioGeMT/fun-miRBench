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
- `data/experiments/processed/18745741/`: local cache for curated benchmark DE tables from Zenodo record `18745741`; the repo currently ships the 3 default benchmark TSVs here
- `data/experiments/raw/`: local raw GEO inputs such as count matrices and FASTQs
- `data/predictions/`: local generated predictor TSVs
- `metadata/mirna_experiment_info.tsv`: experiment registry
- `metadata/predictions_info.tsv`: predictor registry
- `pipelines/experiments/`: experiment-ingestion backend files and example configs
- `pipelines/standardized_predictors/`: predictor pipelines
- `manuscript_assets/`: final manuscript figure exports and compact supporting tables
- `results/`: benchmark outputs

The benchmark reads file paths from the two metadata TSVs. `data/` holds the real files. `results/`
is only for benchmark output.

## Quick Start

If you just want to run the benchmark, you do not need the experiments pipeline first. The repo
already includes:

- experiment metadata in `metadata/mirna_experiment_info.tsv`
- predictor metadata in `metadata/predictions_info.tsv`

Then run the default benchmark:

```bash
uv run funmirbench --config benchmark.yaml
```

Before benchmarking, `funmirbench` syncs the selected curated experiment DE tables from Zenodo
into `data/experiments/processed/18745741/` as needed. The repo currently ships the 3 TSVs used by
the default benchmark config, while other curated benchmark DE tables are treated as fetched local cache.

The default config already points at:

- `metadata/mirna_experiment_info.tsv`
- `metadata/predictions_info.tsv`

and selects a small set of curated experiment datasets plus the registered external predictors.

## Workflow

### 1. Add Experiment Data

The experiment-ingestion pipeline creates DE tables under `data/experiments/processed/` for local
workflow use.

For the curated benchmark datasets tracked in `metadata/mirna_experiment_info.tsv`, the expected
workflow is different: those metadata rows stay versioned in the repo, and the corresponding DE
tables live under the local `data/experiments/processed/18745741/` cache. The repo currently ships
the 3 default benchmark TSVs there, and other curated tables are fetched from Zenodo when needed.

Experiment config summary:

- top-level:
  `dataset_id`, `mirna_name`, `experiment_type`, optional `gse`
- `source`:
  `mode: count_matrix` or `mode: reads`
- `comparison`:
  control vs treated columns or explicit control vs treated samples
- `metadata`:
  fields that will later be synced into `metadata/mirna_experiment_info.tsv`

Supported inputs:

- count matrix: counts matrix + control columns + treated columns -> DESeq2
- reads: local FASTQs + local reference files + explicit sample groups -> `FastQC + fastp + STAR + featureCounts + DESeq2`

RNA-seq pipeline summary:

- `count_matrix` mode validates the configured count matrix and runs DESeq2 directly on the selected
  control and treated columns.
- `reads` mode runs read QC/trimming, aligns reads to the configured reference with STAR, counts
  genes with featureCounts, then runs the same DESeq2 step.
- DE outputs preserve the original method-level `PValue` and `FDR` values. Benchmarking derives
  `benchmark_FDR` and `plot_FDR` internally without modifying those RNA-seq outputs.

This path expects the `funmirbench-experiments` environment from `pipelines/experiments/environment.yml` to be
active so `fastqc`, `fastp`, `STAR`, `featureCounts`, and `Rscript` are available on `PATH`.

Download the shipped real example inputs:

```bash
uv run funmirbench-experiments-download-examples
```

That downloader fetches:

- the real `GSE253003` count matrix
- the real `GSE93717` FASTQ files
- the shared Homo sapiens Ensembl v115 genome FASTA and GTF used by the reads example

Run the real count-matrix example:

```bash
uv run funmirbench-experiments --config pipelines/experiments/configs/gse253003.count_matrix.example.yaml
```

Run the reads example the same way:

```bash
uv run funmirbench-experiments --config pipelines/experiments/configs/gse93717.reads.example.yaml
```

Tracked example configs:

- `pipelines/experiments/configs/gse253003.count_matrix.example.yaml`
- `pipelines/experiments/configs/gse93717.reads.example.yaml`

Reads configs can either:

- use local `reads_1` and optional `reads_2`; each value may be one path or a
  list of technical-run paths for the same biological sample
- use local `genome_fasta_path` and `gtf_path`

Technical runs listed for one sample are combined before QC and trimming.
Single-end and paired-end biological samples may coexist in one comparison;
featureCounts processes the two layouts separately and the pipeline combines
their gene-count columns before DESeq2.

So the practical reads flow is:

1. activate `funmirbench-experiments`
2. run `uv run funmirbench-experiments-download-examples`
3. run `uv run funmirbench-experiments --config pipelines/experiments/configs/gse93717.reads.example.yaml`

The shipped reads example now points at the downloaded Ensembl v115 reference source files under
`data/experiments/raw/refs/ensembl_v115/`, so it builds the derived STAR index automatically.

Each run writes:

- `data/experiments/processed/<dataset_id>.tsv`
- `pipelines/experiments/runs/<timestamp>_<dataset_id>/candidate_metadata.tsv`
- `pipelines/experiments/runs/<timestamp>_<dataset_id>/run_manifest.json`

The processed DE TSV keeps `PValue` and `FDR` unchanged from the DE method, including missing
adjusted p-values and exact zero values. During benchmarking, rows with `PValue` but no adjusted
p-value get `benchmark_FDR = 1`, so they remain non-significant background genes. Rows with both
`PValue` and `FDR` missing remain excluded from FDR-thresholded evaluation. Rows with `FDR = 0`
use the smallest positive floating-point value for `plot_FDR`, computed with `np.nextafter(0, 1)`,
for finite `-log10` plotting while keeping `benchmark_FDR = 0`.

Benchmark DE tables use a canonical schema independent of the DE method:

- `gene_id`: Ensembl gene identifier
- `logFC`: signed log fold change
- `FDR`: adjusted p-value or q-value used for FDR-thresholded evaluation
- `PValue`: optional raw p-value used to classify rows with missing adjusted p-values

Common native headers are normalized when benchmark tables are read. For example, `log2FoldChange`
maps to `logFC`, `padj`, `adj.P.Val`, and `qvalue` map to `FDR`, and `pvalue` or `P.Value` map to
`PValue`.

The reads example uses a reproduced dataset id, `GSE93717_OE_miR_941_deseq2`, so syncing it creates
a separate variant instead of overwriting the curated `GSE93717_OE_miR_941` registry row.

### 2. Sync Experiment Metadata

The ingestion pipeline does not edit `metadata/mirna_experiment_info.tsv` by itself. It writes a
`candidate_metadata.tsv` under `pipelines/experiments/runs/<timestamp>_<dataset_id>/` first.
Then sync it into the registry with:

```bash
uv run funmirbench-sync-metadata
```

That command auto-discovers all `candidate_metadata.tsv` files under `pipelines/experiments/runs/`
and upserts them into the registry. Re-running is safe — existing rows with matching `id` values are
replaced, not duplicated.

To sync a specific file instead:

```bash
uv run funmirbench-sync-metadata --input pipelines/experiments/runs/<run_dir>/candidate_metadata.tsv
```

### 3. Add Predictor Data

Predictor score files live under `data/predictions/` and are discovered through
`metadata/predictions_info.tsv`.

The repo includes standardization pipelines for the curated external predictors under
`pipelines/standardized_predictors/`. Run the predictor-specific pipeline for each score file you
want to generate, then register the resulting standardized TSV in `metadata/predictions_info.tsv`.

### 4. Run The Benchmark

The default benchmark config is `benchmark.yaml`. The primary benchmark universe is
Ensembl protein-coding genes only; FuNmiRBench filters each DE table to that universe
before joining predictor scores.

Benchmark config summary:

- `experiments_tsv`: experiment metadata table
- `predictions_tsv`: predictor metadata table
- `experiments`: which experiment rows to include
- `predictors`: which predictor rows to include
- `evaluation`: thresholds, ranking settings, and the protein-coding benchmark universe
- `tags`: optional labels included in the per-run output folder name
- `out_dir`: results root directory; each benchmark run creates its own subfolder under this root

Run it with:

```bash
uv run funmirbench --config benchmark.yaml
```

That command automatically syncs only the experiment DE tables selected by your benchmark config
from Zenodo into the local `data/experiments/processed/18745741/` cache before joining predictions.
On the first protein-coding run, it also downloads or reuses the Ensembl v115 GTF and caches the
protein-coding gene set at `data/resources/ensembl/protein_coding_gene_ids.txt`.

If you want to prefetch the full curated experiment cache yourself, you can also run:

```bash
uv run funmirbench-experiments-store
```

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
  protein_coding_gtf: data/resources/ensembl/Homo_sapiens.GRCh38.115.gtf.gz
  protein_coding_gene_cache: data/resources/ensembl/protein_coding_gene_ids.txt

out_dir: results/
```

Filter behavior:

- different keys are combined with AND
- values inside one list are combined with OR

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
uv run funmirbench --config benchmark.yaml
uv run funmirbench-experiments-store
uv run funmirbench-validate-experiments --experiments-tsv metadata/mirna_experiment_info.tsv
uv run funmirbench-experiments-download-examples
uv run funmirbench-experiments --config config.yaml
uv run funmirbench-sync-metadata
```

## Tests

```bash
uv run python -m unittest discover -s tests
```
