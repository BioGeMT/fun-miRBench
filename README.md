# fun-miRBench

fun-miRBench evaluates functional miRNA target predictors against gene-expression changes in miRNA overexpression (OE) experiments. It provides a reproducible benchmark run, per-experiment reports, and the pipelines used to build the benchmark tables.

## Install

You need Python 3.10 or later and [`uv`](https://docs.astral.sh/uv/).

```bash
python -m pip install uv
git clone https://github.com/BioGeMT/fun-miRBench.git
cd fun-miRBench
uv sync
```

## Repository contents

- `benchmark.yaml` and `benchmark_publication.yaml`: example and publication run configurations.
- `funmirbench/`: benchmark, report, and data-preparation code.
- `metadata/`: benchmark genes, experiment metadata, predictor registry, and miRNA sequences.
- `pipelines/processing/`: RNA-seq processing and benchmark table construction.
- `pipelines/standardized_predictors/`: pipelines for rebuilding predictor files.
- `benchmark_results/`: committed reference report, plots, and result tables.

Downloaded inputs and intermediate files are kept outside version control: `data/benchmark_data/`, `data/predictions/`, `data/common_resources/`, `data/experiments/`, and `pipelines/processing/work/`. New benchmark runs are written under `results/` by default.

## Quick start

Download the published benchmark tables and predictor scores, then run the example configuration:

```bash
uv run funmirbench-download-data
uv run funmirbench --config benchmark.yaml
```

The command creates a timestamped directory under `results/`. It contains `REPORT.pdf`, the result tables, plots, and a copy of the configuration used. The example config evaluates the five published predictors and creates detailed experiment reports. Edit `benchmark.yaml` to select experiments, predictors, or another output directory.

`funmirbench-download-data` downloads the current 72-experiment benchmark tables and family-expanded TargetScan file from Zenodo [record 23988869](https://zenodo.org/records/23988869). It gets the standardized predictor archive from [record 21671476](https://zenodo.org/records/21671476), which also contains an older collection of differential-expression tables. Existing local files are kept; pass `--force` to replace the files managed by the downloader. To rebuild the benchmark tables or predictor files, see [Build the benchmark tables](#build-the-benchmark-tables) and [Predictors](#predictors).

## Reproduce the publication analysis

`benchmark.yaml` is the general example configuration. Its optional 3′UTR feature comparison is off. `benchmark_publication.yaml` selects the five published predictors, enables that comparison, and disables per-experiment detailed reports. It writes the main report, plots, and tables to `results/publication/<timestamp>/`.

After installing the base environment, install the optional feature dependency and run the publication config:

```bash
uv sync --extra gene-features
uv run funmirbench-download-data
uv run funmirbench --config benchmark_publication.yaml
```

On the first run with gene features enabled, fun-miRBench downloads the Ensembl 115 GTF and hg38 phastCons100way BigWig to the paths in the config. The BigWig is about 5.5 GB; both files are cached there for subsequent runs.

The Intersection Gene Set (IGS) contains benchmark genes with at least one score from every selected predictor across the evaluated miRNAs. The non-IGS group contains the remaining genes in the same benchmark universe. For each gene, the feature table uses the longest 3′UTR among its Ensembl 115 protein-coding transcripts. Conservation is the mean hg38 phastCons100way score over bases with available scores in that UTR. Lengths above 10,000 nt are capped at 10,000 in the plot; tables and statistical tests use full lengths. Genes without an annotated UTR or conservation values are retained, with missing feature values. The report includes group sizes, matched feature counts, medians, median differences, and two-sided asymptotic Mann–Whitney U tests. The summary TSV includes a formatted p-value and log10 p-value so very small results remain readable when the ordinary floating-point p-value is below representable precision. These summaries describe differences between coverage groups.

## Benchmark design

- **Gene universe:** every predictor is evaluated on the same 19,226 genes. A missing predictor score ranks below all scored genes; coverage is reported separately.
- **Labels:** `state_d05` classifies genes as responsive, no response above margin, or unresolved using the shrunken log2 fold change and a margin of 0.5. Unresolved genes are excluded from AUROC and recovery analyses.
- **Per-comparison metrics:** Spearman correlation with the oriented shrunken log2 fold change, overall and within control-expression tenths; and AUROC of responsive versus no-response genes, overall and within expression tenths. AUROC is calculated only when both classes contain at least 30 genes.
- **Decoy comparison:** each predictor is also evaluated using scores for other benchmark miRNAs with a different nt 2–8 seed. `spearman_above_decoys_shrunk` is the cognate miRNA's Spearman correlation minus the mean correlation of its decoys. It is reported in correlation units and is unavailable when no decoys can be scored.
- **Summaries:** metrics are reported over all comparisons and unflagged comparisons. Means and 95% intervals are calculated by bootstrapping dependency groups 2,000 times.
- **Default predictors:** TargetScan family-expanded, miRDB/MirTarget, microT-CNN, miRBind2, and miRAW. The family-expanded TargetScan score is assigned to each miRNA in that family; the unexpanded TargetScan entry remains available in the registry.

## Outputs

Each run writes the following files under its timestamped output directory:

- `REPORT.pdf`: run settings, key results, plots, result tables, and the evaluated and skipped experiments.
- `per_comparison.tsv`, `summary.tsv`, `coverage.tsv`, `recovery.tsv`, and `precision_at_recall.tsv`.
- `benchmark_config.yaml`: the configuration used for the run.
- `plots/`: PNG figures for coverage, performance, and recovery. The recovery plot is omitted when no selected comparison has enough responsive and no-response genes.
- `gene_set_features.tsv` and `gene_set_feature_summary.tsv` when `reports.gene_set_features.enabled` is true.

The run-level plots show:

- **Coverage:** genes scored by each predictor, experiments retained, and responsive and no-response miRNA–gene pairs covered.
- **Performance:** AUROC, Spearman correlation, Spearman correlation above decoys, and the number of experiments each predictor ranks first.
- **Recovery:** recall at prediction budgets and each experiment's best precision at 50% recall relative to random. Genes tied at a budget contribute in proportion to the share of the tie included.
- **Gene features (optional):** 3′UTR length and conservation distributions for IGS and non-IGS genes.

### Detailed reports

Detailed reports are enabled by default and saved in a separate `reports/` directory:

- `reports/experiments/<experiment>/REPORT.pdf`: metrics for the full fixed gene-list pair set (FPS), pairs scored by every predictor (IPS), and each predictor's own scored pairs; precision–recall and ROC curves; enrichment; false-positive burden; top-target effect distributions; recovery; gene heatmaps; and predictor agreement.
- `reports/experiments/<experiment>/metrics.tsv`: metrics for that experiment. `plots/overview/` contains cross-predictor figures; `plots/predictors/<tool>/` contains predictor-specific figures; `plots/data/` contains the plotted values as compressed TSV files. Shared gene scores, ranks, and response annotations are written once per experiment.
- `reports/experiments/<experiment>/predictors/<tool>/REPORT.pdf`: predictor-specific classification, score-versus-response, and enrichment results.

The main `REPORT.pdf` summarizes run-level results. Detailed reports live under `reports/` and are not linked from the main PDF. No cross-experiment predictor reports or supplementary report index are generated. Set `reports.detailed: false` to skip detailed reports; this does not affect the optional gene-feature page in the main report. Detailed reports are generated in parallel; `reports.workers` sets the worker count (default 2), and `reports.top_n` sets the top-effect distribution budget (default 100).

## Predictors

Predictor files are tab-separated tables with `Ensembl_ID`, `Gene_Name`, `miRNA_ID`, `miRNA_Name`, and `Score` columns. They use Ensembl 115 (GRCh38) gene IDs and miRBase 22.1 miRNA names. To add a predictor, place its file under `data/predictions/`, register it in `metadata/predictions_info.tsv`, and select it in `benchmark.yaml`.

The published predictor files can be downloaded with `funmirbench-download-data`. To rebuild one, run `uv run python pipelines/standardized_predictors/<tool>/pipeline.py`. Downloads and shared resources are cached under `data/common_resources/` and the predictor's directory.

- **TargetScan:** TargetScan v8 cumulative weighted context++ score; lower scores indicate stronger predictions. `--expand-family-members` assigns each family score to every family member and produces the `targetscan_family_expanded` input used by the default configs.
- **miRDB/MirTarget:** miRDB 6.0 / MirTarget 4.0 scores mapped to genes using NCBI gene IDs, RefSeq IDs, and BioMart. Higher scores are stronger.
- **microT-CNN:** scores distributed through Zenodo and mapped from transcripts to genes using Ensembl 115. Higher scores are stronger.
- **miRBind2:** miRBind2-3UTR predictions from [Zenodo record 20609975](https://zenodo.org/records/20609975). Lower signed scores are stronger.
- **miRAW:** site-level predictions from [Figshare article 32982218](https://doi.org/10.6084/m9.figshare.32982218). The highest site score is kept for each miRNA–gene pair. Higher scores are stronger.

## Build the benchmark tables

The RNA-seq processing pipeline rebuilds the per-gene tables from reads. It requires conda for the aligners and R; see the [processing README](pipelines/processing/README.md). To add one experiment without rebuilding the full dataset, see [Add your own experiment](pipelines/processing/ADD_EXPERIMENT.md).

## Tests

```bash
uv run python -m unittest discover -s tests
```

## License

The code is licensed under the [MIT License](LICENSE). The published benchmark tables on Zenodo are licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
