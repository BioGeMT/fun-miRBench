# FuNmiRBench

Benchmark functional miRNA target predictors on miRNA overexpression (OE) experiments.

## Install

You need Python 3.10+ and `uv`.

```bash
python -m pip install uv
git clone git@github.com:BioGeMT/FuNmiRBench.git
cd FuNmiRBench
uv sync
```

## Repo Layout

- `benchmark.yaml`: the benchmark config
- `funmirbench/`: the code
- `metadata/`: the comparisons, gene list, sample sheets, predictor list and miRNA sequences
- `pipelines/processing/`: from reads to the per-gene benchmark tables
- `pipelines/standardized_predictors/`: one pipeline per predictor
- `tests/`: the tests
- `benchmark_results/`: `REPORT.pdf`, plots and tables of our benchmark run on the 21 OE test experiments

Not tracked by git:

- `data/benchmark_data/`: the per-gene benchmark tables, the benchmark's input
- `data/predictions/`: predictor scores, the benchmark's input
- `data/common_resources/`: downloaded miRBase and Ensembl annotation
- `data/experiments/`: reads, references and DE tables of `pipelines/processing/`
- `results/<date_time>/`: one folder per benchmark run

## Run The Benchmark

`funmirbench` scores predictors on the overexpression (OE) comparisons of the per-gene benchmark tables, on all genes. Other comparisons in the tables, such as KOs, are not scored; a warning names them.

```bash
uv run funmirbench-download-data
uv run funmirbench --config benchmark.yaml
```

`funmirbench-download-data` puts the benchmark tables in `data/benchmark_data/` (Zenodo [10.5281/zenodo.23824917](https://doi.org/10.5281/zenodo.23824917)) and the predictor scores in `data/predictions/`, and keeps files that are already there. To build them yourself, use `pipelines/processing/` for the tables and `pipelines/standardized_predictors/targetscan/pipeline.py --expand-family-members` for the TargetScan scores.

`benchmark.yaml` sets the benchmark tables (`input_dir`), which experiments and predictors to run, and `out_dir`. Each run creates its own folder under `out_dir`, with a copy of the config.

- **Labels:** `state_d05`: responsive, no response above margin, or unresolved at d = 0.5, from the shrunken log2FC (FSOS s-value < 0.05 against d; whole shrunken 95% interval below d). Unresolved genes are never negatives.
- **Genes:** every ranking is scored on the same 19,226-gene list. A gene a tool does not score ranks last; `coverage` is the share of genes it scores.
- **Scores per comparison:**
  - Spearman between the score and the oriented shrunken log2FC, overall and within tenths of control expression;
  - AUROC of responsive against no-response genes, plain and within tenths of control expression (Mann-Whitney over pairs in the same tenth), only where each state has at least 30 genes.
- **Decoy check:** for each tool, the same Spearman for its scores of the other benchmark miRNAs with a different nt 2–8 seed (`metadata/benchmark_mirnas.tsv`). `decoy_z_shrunk` is the cognate miRNA against its decoys.
- **Summaries:**
  - on all comparisons and on unflagged ones;
  - mean over comparisons with a 95% interval from 2,000 bootstraps over dependency groups;
- **Predictors:** the default config runs the five published tools. TargetScan is `targetscan_family_expanded`, which gives every family member the family score; plain `targetscan` stays in the registry and can be selected instead.
- **Plots**, in `plots/` as PNG:
  - `coverage`: genes each tool scores, the experiments it is retained in, and the responsive and no-response pairs it covers;
  - `oe_all_genes`: AUROC at d = 0.5, Spearman with the shrunken log2FC, decoy z-score and the experiments each tool wins, for every tool on all genes, with y axes that follow the values of the run;
  - `recovery`: recall by prediction budget and each experiment's best precision at 50% recall relative to random, on the comparisons with at least 30 genes in each state.
- **Report:** `REPORT.pdf`: run settings, key results, the mean coverage, AUROC and Spearman per tool, the experiments evaluated and skipped, and the plots.

Outputs: `REPORT.pdf`, `per_comparison.tsv`, `summary.tsv`, `coverage.tsv`, `recovery.tsv`, `precision_at_recall.tsv`, `benchmark_config.yaml`, `plots/`.

## Predictors

A standardized score file is a TSV with the columns `Ensembl_ID`, `Gene_Name`, `miRNA_ID`, `miRNA_Name` and `Score`, on Ensembl 115 (GRCh38) and miRBase 22.1. To add your own predictor, put its file under `data/predictions/`, register it in `metadata/predictions_info.tsv`, and add it to `predictors` in `benchmark.yaml`.

`funmirbench-download-data` fetches the scores of the five published tools. To rebuild one, run `uv run python pipelines/standardized_predictors/<tool>/pipeline.py`; downloads are cached under `data/common_resources/` and the tool's folder.

- **`targetscan`:** TargetScan v8 `Summary_Counts.all_predictions.txt`, human rows on representative transcripts that agree with Ensembl 115. `Score` is the cumulative weighted context++ score; lower is stronger. `--expand-family-members` gives each family score to every member miRNA without recomputing it, which writes the `targetscan_family_expanded` file used in the benchmark.
- **`mirdb_mirtarget`:** miRDB v6.0 (MirTarget v4.0) all-score file, mapped to genes through the NCBI gene ID, then the RefSeq ID, with BioMart. Higher is stronger.
- **`microt_cnn`:** microT-CNN all-score file from Zenodo, with transcripts mapped to genes through Ensembl 115. Higher is stronger.
- **`mirbind2`:** miRBind2-3UTR predictions from Zenodo record 20609975. The score is a signed repression prediction; lower is stronger.
- **`miraw`:** miRAW site-level predictions from Figshare (10.6084/m9.figshare.32982218), keeping the highest score per miRNA and gene. Higher is stronger.

## Build The Benchmark Tables

`pipelines/processing/` goes from reads to the tables, for anyone who wants to rebuild them or add a comparison. It needs `conda` for the aligners and R; see its [README](pipelines/processing/README.md).

## Tests

```bash
uv run python -m unittest discover -s tests
```
