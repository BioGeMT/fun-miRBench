# Single-cell pipeline

Turns a single-cell RNA-seq experiment in which one miRNA was raised (mimic, overexpression) or lowered (inhibitor)
into the same per-gene table as a bulk experiment. The cells of each sample are filtered and summed into one
pseudobulk sample, and the summed counts go through the bulk DESeq2 and labelling steps
(`../processing/`) unchanged.

## Steps

1. **Config.** One YAML per comparison in `configs/`: the miRNA, the experiment type and one entry per sample with
   its group (`control` or `treated`), an optional batch, and its 10x files (`matrix`, `features`, `barcodes`; local
   paths, or URLs that are downloaded into `data/single_cell/raw/<dataset_id>/`). A sample is one culture, donor,
   animal, day or hashed arm. For arms pooled in one lane, `cell_labels` selects the sample's cells from a label
   table (`path`, `barcode_column`, `label_column`, `label`, and optionally `label_regex`, whose first group is the
   label). Matrices with gene symbols instead of Ensembl IDs need `feature_ids: symbol` and `gene_symbols_gtf` (the
   Ensembl 115 GTF); symbols that name more than one gene are left out.
2. **Pseudobulk.** `uv run funmirbench-single-cell --config pipelines/single_cell/configs/<id>.yaml` keeps the cells
   that pass the `qc` thresholds (`min_genes`, `min_counts`, `max_counts`, `max_mito_fraction`; mitochondrial genes are
   the symbols starting with `qc.mito_prefix`, default `MT-`), sums their counts per sample, and fails if a sample keeps
   fewer than `min_cells_per_sample` cells (default 20). Features that are not human Ensembl genes (for example a
   reporter transgene) are left out and counted in `sample_summary.tsv`, and only genes listed in every sample's
   features file are kept (samples quantified against different annotations). Fractional counts (multi-mapping
   reads shared between genes) are summed and then rounded. Doublets are not removed: both cells of a doublet come
   from the same sample, so they change its depth, not its composition.
   With `pseudo_replicates: k` the kept cells of each sample are split at random into k columns. These are technical
   replicates (cell sampling, not culture-to-culture variation), for designs with one sample per arm, which DESeq2
   cannot fit otherwise.
3. **DESeq2.** Add `--run-de` (with `Rscript`, DESeq2 and apeglm on `PATH`, as for `../processing/`) to run
   `funmirbench-experiments` on the written `experiments_config.yaml`. When the samples carry different batches the
   design is `~ batch + condition` (for example one control and one treated sample per day).
4. **Labels.** `../processing/src/build_benchmark_data.py` on the fit (`complete_analysis.rds` in the DE run directory), with
   the comparison listed in a comparisons table, gives the final table with the shrunken log2FC and the states.

## Outputs (`runs/<stamp>_<dataset_id>/`, not tracked)

| File | Content |
|---|---|
| `pseudobulk_counts.tsv.gz` | summed counts, unversioned Ensembl gene IDs x samples (or pseudo-replicates) |
| `design.tsv` | column, sample, group and batch of each count column |
| `cell_qc.tsv.gz` | every cell: sample, barcode, total counts, detected genes, mitochondrial fraction, passed |
| `sample_summary.tsv` | per sample: group, batch, cells in the file and kept, summed counts, median counts per cell, features left out |
| `experiments_config.yaml` | the count_matrix-mode config for `funmirbench-experiments` |
| `run_manifest.json` | inputs, thresholds and, with `--run-de`, the DE run directory and table |

The DE table is written to `data/experiments/processed/<dataset_id>.tsv` as for bulk experiments.

## Selected experiments

`selected_comparisons.tsv` lists the 70 comparisons used for evaluation (`main` 5, `machinery` 62, `supplement` 3), `builds.tsv` the 18 datasets they come from and `sources.tsv` the 24 source files, with their addresses, sizes and MD5 where known. These commands install, download and run them, in an environment with `conda` and `uv`:

```bash
pipelines/single_cell/setup.sh
pipelines/single_cell/download.sh --groups main
pipelines/single_cell/run.sh --groups main
```

Without `--groups`, both scripts take all three groups. `main` is about 2 GB of downloads, `supplement` about 3 GB, `machinery` about 228 GB. Complete files and finished datasets are kept, so both scripts can be rerun after an interruption. `SINGLE_CELL_DATA` and `SINGLE_CELL_FINAL` set the download and output folders (default `data/single_cell/raw` and `data/single_cell/final`), `RUN_WORKERS` the fits at a time per dataset (default 8). The result is `<final>/comparisons.tsv` with a `group` column, one folder per dataset, and `genes/<comparison_id>.tsv.gz` inside it.

## Evaluation

[EXPERIMENTS.md](EXPERIMENTS.md) lists every selected experiment and how it was processed. To score the predictors on the tables that `run.sh` writes, set the paths in `single_cell_eval.yaml` at the repository root and run:

```bash
uv run funmirbench-single-cell-eval --config single_cell_eval.yaml
```

It needs `<final>/comparisons.tsv`, the standardized predictions and `per_comparison.tsv` of a `funmirbench` run on the bulk experiments. `evaluation/` holds its fixed tables: `mirna_profiles.tsv` (miRNA family shares per cell line), `mirna_families.tsv` (miRBase 22.1), `experiment_labels.tsv`, and `seed_check.tsv`, `seed_check_bulk.tsv` and `control_check.tsv`, which were computed separately because they need the 3'UTR sequences and the cell stores. The output is `REPORT.pdf`, with two figures of grouped bars and tables of the same numbers, the figures in `plots/`, `per_comparison.tsv`, `machinery_per_comparison.tsv`, `machinery_family_coverage.tsv`, `data_checks.tsv` (with `in_figure`), and the bulk evaluator's own run in `evaluator/`.

## Final data

`funmirbench-single-cell-final` builds the final tables of every dataset from one cell store and one comparison list. The per-experiment steps above stay for single configs. Adding a dataset takes two files and no change anywhere else.

### Add a dataset

1. **Adapter**, `adapters/<dataset>.py`: reads the deposit and calls `write_cell_store` from `funmirbench/single_cell_store.py`, which maps features to Ensembl 115 gene IDs through the GTF, validates and writes `<dataset>.h5ad`. The same module has `read_10x_mtx`, `read_10x_h5`, `raw_counts_from_lognorm` (HKU files), `pseudo_replicate_units` and `assign_lane_sets`.
2. **Comparison list**, `comparisons/<dataset>.tsv`, one row per comparison.

Then run, in an environment with DESeq2 and apeglm on `PATH` (`environment.yml` pins it):

`uv run funmirbench-single-cell-final --store <dataset>.h5ad --comparisons pipelines/single_cell/comparisons/<dataset>.tsv --genes metadata/benchmark_genes.tsv.gz --out <final>/<dataset> --workers 8`

The output folder must not exist yet. The command exits non-zero only on infrastructure errors: a bad store or list, no R, or a comparison whose R process was killed.

### Cell store

| Part | Content |
|---|---|
| `X` | raw integer counts, cells x features, sparse; a deposit with only fractional counts (reads shared between genes) keeps them, with `uns["fractional_counts"] = True` |
| `var` | index = feature as deposited; `gene_id` = unversioned Ensembl 115 gene ID where the feature maps uniquely, otherwise empty; `symbol` |
| `obs` | `dataset`, `sample`, `batch`, `replicate_unit`, `replicate_type` (`biological`, `lane_set`, `technical`), `arm`, `is_control`, `control_type` (empty for cells that are not controls), `guide`, `n_guides`, `passed_qc` |
| `uns["source"]` | where the counts come from; copied to every comparison |

Fractional counts are rounded after summing, once for each replicate unit of a comparison, as DESeq2 takes integers.

For screens, store only cells of the listed targets plus all control cells of the batches where those cells occur.

Replicate units, in order of preference: real replicates (donor, independent culture); lane sets for screens (`assign_lane_sets`, with `replicate_type = lane_set`); with one sample per arm, three pseudo-replicates (`pseudo_replicate_units`, with `replicate_type = technical`), whose labels mean "changed in this culture".

At high guide load (median `n_guides` of the passing cells above 1.5), `arm` lists the targets of the cell's guides, separated by `;`, for example `DICER1;DROSHA`.

### Comparison list

Columns: `comparison_id`, `dataset`, `category` (`raise`, `lower`, `pooled_mirna`, `pathway`, `host_gene`, `reporter`, `viral_mirna`, `ko`, `other`), `mirna_or_gene`, `orientation` (`OE` when the miRNA activity goes up, `KO` when it goes down), `treated_arm`, `control_arm`, `min_cells` (per replicate unit), `caveats`. Optional: `other_pathway_arms`, separated by `;`, used only at high guide load; `target_gene_id` (Ensembl IDs separated by `;`, the first in the store counts), used for the target's own change.

Eligibility comes from design and data only: no seed or target signal is used to keep or drop a comparison.

### What a comparison does

1. Passing cells are summed per `batch`, `replicate_unit` and arm. These sums are computed once per store and shared, so control cells are not read again for each comparison.
2. Units below `min_cells` are dropped, then units in batches that do not hold both arms.
3. With at least two units per side, the sums go to `run_deseq2_counts.R` (`~ batch + condition` when batches differ, otherwise `~ condition`), then `label_states.R` with the gene list.

At high guide load, treated is the cells carrying the target's guides. Control is the cells of the same batch carrying no guide of the target and no guide of `other_pathway_arms`, and `control_arm` is not used. The caveat is added to the comparison.

### Outputs

| File | Content |
|---|---|
| `comparisons.tsv` | one row per comparison: the list's columns, `source`, cells and units per side, units below `min_cells`, `replicate_type`, `control_type`, `control_selection`, `median_guides_per_cell`, `design`, `status`, `reason`, `table`, gene and d = 0.5 state counts; with `target_gene_id`, the target's own change as information only: `target_in_store`, `target_used_gene_id`, `target_treated_sum`, `target_control_sum`, `target_log2_ratio_of_sums`, `target_log2fc`, `target_log2fc_shrunk` |
| `genes/<comparison_id>.tsv.gz` | one row per listed gene. `fitted`: the table described in `../processing/README.md`. Otherwise: `treated_sum`, `control_sum`, `log2_ratio_of_sums` (library-normalised, 0.5 added to each sum), `analysis_status = not_fitted` |
| `comparisons/<comparison_id>/` | `pseudobulk_counts.tsv.gz`, `design.tsv`, and for attempted fits the DESeq2 outputs and label logs |
| `units.tsv` | `comparison_id`, `ok`, `skipped` or `failed`, reason, as `long-runs` reads it (`--units` writes it to the run directory) |
| `manifest.json` | inputs with hashes, status counts |

`status` is `fitted`, `not_fitted` (too few cells or units, no shared batch; the reason says which), `fit_failed` (DESeq2 or `label_states.R` stopped; the reason holds its last lines) or `error` (infrastructure). All but `fitted` keep the ratio table if the comparison got that far. `not_fitted` and `fit_failed` are results and leave the exit code at 0.

`funmirbench-single-cell-combine --root <final>` writes `<final>/comparisons.tsv`, the concatenated `comparisons.tsv` of the dataset folders with `final_dir` added. A folder with a `SUPERSEDED` file, naming its replacement, stays out.

### Relabel

`uv run python -m funmirbench.single_cell_relabel --root <final> --out <new root> --genes metadata/benchmark_genes.tsv.gz --workers 24` runs `label_states.R` again on the saved fits, without reading a store or fitting anything. Only the `svalue_d*` and `state_d*` columns may change; any other difference makes the comparison an `error`. Finished tables under `--out` are kept, so a second run adds what is new. Stores and fits stay in `--root`; `manifest.json` of each dataset names them.

### Screens

`adapters/screen_common.py` holds the shared steps and `adapters/target_genes.tsv` the targets: 19 pathway genes, non-coding hosts, miRNA loci and protein-coding hosts, with the aliases safe to match. Protein-coding hosts go to `comparisons/<dataset>_pc_hosts.tsv`, run from the same store into `<final>/<dataset>_pc_hosts/`. `adapters/run_screen.sh` runs one screen: adapter, one smoke fit, all fits, outcome row.

### Co-occurrence datasets

Small RNA and mRNA of the same cells use their own step, no comparison list:

`uv run funmirbench-single-cell-cooccurrence --mirna <cells x miRNAs> --mrna <cells x genes> --cells <cells.tsv> --group-column cell_line --out <final>/<dataset>`

`--mirna` and `--mrna` are raw integer counts as h5ad (`var.gene_id` is used where present) or TSV with the cell in the first column; `--allow-fractional` accepts non-integer values (multi-mapper shares, normalised deposits). `cells.tsv` has `cell` and the grouping column. `--library-totals` is a TSV with `cell` and `library_total`, the sum of every deposited feature of the library that holds the miRNA features. Genes are divided by the cell's mRNA table total. miRNAs are divided by the library total, and cells without a miRNA count stay as zeros; without `--library-totals`, miRNAs are divided by the cell's miRNA table total and those cells are left out. `groups.tsv` records the rule per group. Then within each group, for each miRNA and gene detected in at least 10% of the group's cells (`--min-detected-fraction`), the Spearman correlation across cells and its Bonett-Wright 95% interval are computed. Groups below 10 cells (`--min-cells`) are listed as `not_correlated`.

| File | Content |
|---|---|
| `groups.tsv` | per group: cells, cells left out, cells with no miRNA count, expressed miRNAs and genes, status, reason, miRNA normalisation |
| `expressed_mirnas.tsv` | per group: expressed miRNAs, cells detected, fraction detected, total counts |
| `correlations/<group>.tsv.gz` | `mirna`, `gene`, `spearman_r`, `ci95_low`, `ci95_high`, `n_cells` |
| `per_cell/<group>_mirna.tsv.gz`, `<group>_mrna.tsv.gz` | raw counts of the expressed features per cell, with `library_size`, the denominator used |
