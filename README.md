# FuNmiRBench

Benchmark functional miRNA target predictors on miRNA overexpression (OE) experiments.

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

- `benchmark.yaml`: the benchmark config
- `funmirbench/`: the code
- `metadata/`: lists of experiments, predictors and miRNA sequences
- `pipelines/`: pipelines that make experiment DE tables and predictor scores
- `data/benchmark_data/`: the per-gene benchmark tables, the benchmark's input
- `data/predictions/`: predictor scores, the benchmark's input
- `data/common_resources/`: downloaded miRBase and Ensembl annotation
- `data/experiments/`: raw experiment inputs and DE tables of `pipelines/experiments/`
- `results/benchmark_results/`: one folder per benchmark run
- `benchmark_results/`: `REPORT.pdf`, plots and tables of our benchmark run on the 21 OE test experiments

## Run The Benchmark

`funmirbench` scores predictors on the overexpression (OE) comparisons of the per-gene benchmark tables, on all genes. Other comparisons in the tables, such as KOs, are not scored; a warning names them.

```bash
uv run funmirbench-download-data
uv run funmirbench --config benchmark.yaml
```

`funmirbench-download-data` puts the benchmark tables in `data/benchmark_data/` (Zenodo [10.5281/zenodo.23824917](https://doi.org/10.5281/zenodo.23824917)) and the predictor scores in `data/predictions/`, and keeps files that are already there. To build them yourself, use `pipelines/benchmark_data/` for the tables and `pipelines/standardized_predictors/targetscan/pipeline.py --expand-family-members` for the TargetScan scores.

`benchmark.yaml` sets the benchmark tables (`input_dir`), which experiments and predictors to run, and `out_dir`. Each run creates its own folder under `out_dir`, with a copy of the config.

- **Labels:** `state_d05`: responsive, no response above margin, or unresolved at d = 0.5, from the shrunken log2FC (FSOS s-value < 0.05 against d; whole shrunken 95% interval below d). Unresolved genes are never negatives.
- **Genes:** every ranking is scored on the same 19,226-gene list. A gene a tool does not score ranks last; `coverage` is the share of genes it scores.
- **Scores per comparison:**
  - Spearman between the score and the oriented shrunken log2FC, overall and within tenths of control expression;
  - AUROC of responsive against no-response genes, plain and within tenths of control expression (Mann-Whitney over pairs in the same tenth), only where each state has at least 30 genes.
- **Reference rows:**
  - `random`;
  - `expression_high` / `expression_low` (control expression).
- **Decoy check:** for each tool, the same Spearman for its scores of the other benchmark miRNAs with a different nt 2–8 seed (`metadata/benchmark_mirnas.tsv`). `decoy_z_shrunk` is the cognate miRNA against its decoys.
- **Summaries:**
  - on all comparisons and on unflagged ones;
  - mean over comparisons with a 95% interval from 2,000 bootstraps over dependency groups;
- **Predictors:** the default config runs the five published tools. TargetScan is `targetscan_family_expanded`, which gives every family member the family score; plain `targetscan` stays in the registry and can be selected instead.
- **Plots**, in `plots/` as PNG and SVG:
  - `coverage`: genes each tool scores, the experiments it is retained in, and the responsive and no-response pairs it covers;
  - `oe_all_genes`: AUROC at d = 0.5, Spearman with the shrunken log2FC and the experiments each tool wins, for every tool on all genes;
  - `recovery`: recall by prediction budget and each experiment's best precision at 50% recall relative to random, on the comparisons with at least 30 genes in each state.
- **Report:** `REPORT.pdf`: run settings, key results, the mean coverage, AUROC and Spearman per tool, the experiments evaluated and skipped, and the plots.

Outputs: `REPORT.pdf`, `per_comparison.tsv`, `summary.tsv`, `coverage.tsv`, `recovery.tsv`, `precision_at_recall.tsv`, `benchmark_config.yaml`, `plots/`.

## Add Your Own Data

- **Predictors:** put a standardized score TSV under `data/predictions/`, register it in `metadata/predictions_info.tsv`, and add it to `predictors` in `benchmark.yaml`. `pipelines/standardized_predictors/` holds the pipelines of the five published tools.
- **Experiments:** [`pipelines/experiments/`](pipelines/experiments/README.md) runs DESeq2 from a YAML config and syncs the result into `metadata/mirna_experiment_info.tsv`.

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
