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
- `data/experiments/processed/21671476/`: local cache for curated benchmark DE tables from Zenodo record `21671476`; the repo currently ships the 3 default benchmark TSVs here
- `data/experiments/raw/`: local raw GEO inputs such as count matrices and FASTQs
- `data/common_resources/`: shared, downloaded miRBase and Ensembl annotation cache
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
into `data/experiments/processed/21671476/` as needed. The repo currently ships the 3 TSVs used by
the default benchmark config, while other curated benchmark DE tables are treated as fetched local cache.

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

For the curated benchmark datasets tracked in `metadata/mirna_experiment_info.tsv`, the expected
workflow is different: those metadata rows stay versioned in the repo, and the corresponding DE
tables live under the local `data/experiments/processed/21671476/` cache. The repo currently ships
the 3 default benchmark TSVs there, and other curated tables are fetched from Zenodo when needed.

### 2. Add Predictor Data

Predictor score files live under `data/predictions/` and are discovered through
`metadata/predictions_info.tsv`.

The repo includes standardization pipelines for the curated external predictors under
`pipelines/standardized_predictors/`. Run the predictor-specific pipeline for each score file you
want to generate, then register the resulting standardized TSV in `metadata/predictions_info.tsv`.

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
- `tags`: optional provenance labels recorded in the run `summary.json`
- `out_dir`: results root directory; each benchmark run creates its own subfolder under this root

Run it with:

```bash
uv run funmirbench --config benchmark.yaml
```

That command automatically syncs only the experiment DE tables selected by your benchmark config
from Zenodo into the local `data/experiments/processed/21671476/` cache before joining predictions.
On the first protein-coding run, it also downloads or reuses the Ensembl v115 GTF and caches the
protein-coding gene set at `data/common_resources/ensembl/protein_coding_gene_ids.txt`.

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

If validation finds that an experiment has no ground-truth positives at the configured thresholds,
the benchmark records it as skipped and continues (unless `zero_positive_policy: continuous_only`
is set; see Evaluation Options below). The same happens when an experiment has no rows
scored by the selected predictors, or no predictor remains evaluable after DE and ground-truth
filtering. Experiments that reach predictor joining retain their `joined.tsv` for diagnosis. Every
skip reason is reported in the run `README.md`, `REPORT.pdf`, and `summary.json`. A run in which every
selected experiment is skipped still produces these three summary files; metric tables and
evaluation plots are produced only for experiments that can be evaluated.

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

### Evaluation Options

The defaults of these optional `evaluation` keys retain the legacy GT label rule and dataset
selection. Other changes are active regardless of these keys:

- Rank repair: genes without a score keep a missing dataset-local rank. Previously, if a predictor had a
  single distinct score in a dataset, every gene, including unscored ones, received rank 1.0. This can
  change combination, heatmap, and rank-plot outputs for such predictors.
- Added outputs: DE covariate columns in `joined.tsv`, continuous metrics for every evaluated
  predictor and reference-baseline tables.

The evaluator still uses the legacy binary rule (`FDR < fdr_threshold` and oriented effect
`> effect_threshold`; everything else usable is a negative). Confidence-aware labels, for example
the `label3` columns in the v2 candidate DE tables, are not used yet; integrating them is future
work.

Options:

- `fdr_na_policy` (`negative` | `exclude`, default `negative`): genes with a missing FDR but a
  reported PValue (for example DESeq2 independent filtering). `negative` treats them as FDR 1.0,
  so they can become GT negatives; `exclude` leaves them out of both GT classes. Excluded genes are
  counted per predictor as `rows_excluded_fdr_na` / `rows_excluded_fdr_na_scored`.
- `zero_positive_policy` (`exclude` | `continuous_only`, default `exclude`): experiments without GT
  positives. `exclude` skips them during validation. `continuous_only` evaluates them with
  Spearman/Pearson only; AUROC/APS/PR-AUC are NA and `class_metrics_na_reason` says why. The same
  applies to a predictor whose scored genes contain no positives (or no negatives): by default it
  is listed as skipped, with `continuous_only` it gets a metric row with NA class metrics.
  Spearman of score vs expected effect is always written, for every predictor with scored genes in
  every experiment that reaches evaluation (including predictors without positives), to
  `tables/per_experiment/continuous_metrics_per_experiment.tsv`.
- `reference_baselines` (default `true`) and `reference_seed` (default `0`): reference scorers
  computed inside the evaluator from each DE table: `reference_random` (seeded, reproducible,
  independent of row order) and `reference_expression_high` / `reference_expression_low`
  (`control_mean_normalized_count`, descending and ascending). They use the predictors' GT labels
  and are scored on every usable gene. They never enter predictor metric tables, the per-predictor
  cross-dataset summary, or plots.

DE covariates present in the DE table (`baseMean`, `lfcSE`, `stat`,
`control_mean_normalized_count`, `treated_mean_normalized_count`, edgeR `logCPM`/`F`, the
normalized/raw group means, and any `*status` column) are carried into `joined.tsv` for
stratified diagnostics. They are never treated as predictor scores; only `score_<tool_id>`
columns from predictor files are scored.

Registry `experiment_type` values are case-insensitive and accept the abbreviations `OE`, `KO`,
and `KD` for `Overexpression`, `Knockout`, and `Knockdown`.

Filter behavior:

- different keys are combined with AND
- values inside one list are combined with OR

`benchmark.yaml` already includes other experiment and predictor columns as commented rows, so the
normal workflow is just to edit or uncomment filters.

## v2 Evaluation (three-state labels)

`funmirbench-v2-evaluate` scores tools on the overexpression (OE) comparisons of the v2 final tables made by the `v2-processing` branch (`pipelines/final_data/`), on all genes. KO comparisons in the tables are not scored.

```bash
uv run funmirbench-download-data                       # published standardized predictors
uv run python pipelines/standardized_predictors/targetscan/pipeline.py --expand-family-members
uv run funmirbench-v2-evaluate --final-dir <final-v2> --utr-fasta <MANE Select 3'UTR FASTA> --out results/v2/
```

- **Labels:** `state_d05` / `state_d10`: responsive, no response above margin, or unresolved, from the shrunken log2FC (FSOS s-value < 0.05 against d; whole shrunken 95% interval below d). Unresolved genes are never negatives. The raw states `state_raw_d05` / `state_raw_d10` are scored too, as a check (metrics named `*_raw_*`).
- **Genes:** every ranking is scored on the same 19,226-gene list. A gene a tool does not score ranks last; `coverage` is the share of genes it scores.
- **Scores per comparison:**
  - Spearman between the score and the oriented log2FC (raw and shrunken), overall and within tenths of control expression;
  - AUROC of responsive against no-response genes, plain and within tenths of control expression (Mann-Whitney over pairs in the same tenth), only where each state has at least 30 genes.
- **Reference rows:**
  - `random`;
  - `expression_high` / `expression_low` (control expression);
  - `seed_site`: 3 × 8mer + 2 × 7mer-m8 + 7mer-A1 + 0.5 × 6mer in the MANE 3′UTR;
  - `seed_plus_expression`: seed score with expression breaking ties.
- **Decoy check:** for each tool and the two seed baselines, the same Spearman for its scores of the other benchmark miRNAs with a different nt 2–8 seed (`metadata/v2/benchmark_mirnas.tsv`). `decoy_z` (raw log2FC) and `decoy_z_shrunk` are the cognate miRNA against its decoys.
- **Summaries:**
  - on all comparisons and on unflagged ones;
  - mean over comparisons with a 95% interval from 2,000 bootstraps over dependency groups;
  - `paired_vs_reference.tsv` compares each ranking with `seed_plus_expression` on the same comparisons.
- **Predictors:** `metadata/v2/predictions.tsv` lists the five published tools plus `targetscan_family_expanded`, which gives every family member the family score.

- **Plots**, in `plots/` as PNG and SVG:
  - `coverage`: genes each tool scores, the experiments it is retained in, and the responsive and no-response pairs it covers;
  - `oe_all_genes`: AUROC at d = 0.5, Spearman with the shrunken log2FC and the experiments each tool wins, for every tool on all genes;
  - `recovery`: recall by prediction budget and each experiment's best precision at 50% recall relative to random, on the comparisons with at least 30 genes in each state.

  `targetscan_family_expanded` takes the place of `targetscan` in the plots.
- **Report:** `REPORT.pdf`: run settings, key results, the mean coverage, AUROC and Spearman per tool, the experiments evaluated and skipped, and the plots.

Outputs: `REPORT.pdf`, `per_comparison.tsv`, `decoy.tsv`, `summary.tsv`, `paired_vs_reference.tsv`, `coverage.tsv`, `recovery.tsv`, `precision_at_recall.tsv`, `run_info.json`, `plots/`.

`python -m funmirbench.v2_margin_report --final-dir <final dir> --utr-fasta <MANE Select 3'UTR FASTA> --out <dir>` counts the states for every d of the final tables, the genes that moved the other way by more than d (the responsive rule mirrored), their ratio to the responsive genes, and the share of genes with a seed site on each side (`margin_grid.tsv`, and `margin_summary.tsv` at d = 0.5).

## Manuscript Assets

Manuscript figure scripts live under `scripts/`; see `scripts/README.md` for
the exact commands. The repo tracks final manuscript-facing assets rather than
all generated intermediates. For Figure 2, this means tracking the combined
figure exports and compact TSVs that support the manuscript values. Individual
panel image exports, large raw support tables such as gene-level conservation,
downloaded resources, diagnostic outputs, and benchmark run folders are local
generated artifacts and remain gitignored.

## Outputs

After a benchmark run, `results/` contains one new timestamped run folder. For example:

- `results/20260510_143022/`

Detailed dataset, miRNA, predictor, perturbation, cell-line, evaluation-threshold, and
protein-coding-universe metadata is recorded inside the run folder.

Inside each run folder you get:

- `README.md`: human-readable run guide and quick-start map for the output folder
- `REPORT.pdf`: main run-level PDF report with explanations and selected combined plots
- `datasets/<dataset_id>/joined.tsv`: joined DE and predictor scores, retained for evaluated datasets
  and skipped datasets that reached predictor joining
- `datasets/<dataset_id>/plots/predictors/<tool_id>/`: per-tool plots for evaluated datasets
- `datasets/<dataset_id>/plots/comparisons/`: multi-predictor comparison plots for evaluated datasets
- `datasets/<dataset_id>/plots/heatmaps/`: heatmaps for evaluated datasets
- `datasets/<dataset_id>/reports/`: per-dataset Markdown/PDF reports and correlation TSVs for
  evaluated datasets
- `tables/per_experiment/`: per-experiment metric tables when at least one dataset is evaluated,
  plus `continuous_metrics_per_experiment.tsv` and `reference_baselines_per_experiment.tsv`
- `tables/combined/`: cross-dataset predictor summaries when evaluation metrics are available.
  `cross_dataset_predictor_summary.tsv` averages each predictor over its own scored genes and
  evaluated datasets. `cross_dataset_reference_summary.tsv` averages each reference over its own
  datasets.
- `plots/combined/metrics/`, `plots/combined/ranks/`, `plots/combined/combinations/`:
  cross-dataset comparison plots grouped by theme when evaluation metrics are available
- `summary.json`: machine-readable run summary, including configured tags and skipped-dataset reasons

When 2 or more predictors are selected, each evaluated dataset gets:

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
