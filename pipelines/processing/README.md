# Build the benchmark tables

This pipeline processes the 72 published miRNA overexpression (OE) comparisons from FASTQ files through DESeq2 to the per-gene benchmark tables. The full read set contains 738 FASTQ files (about 1.3 TB before intermediate files are written).

## Run

```bash
pipelines/processing/setup.sh
pipelines/processing/download.sh
pipelines/processing/run.sh
```

- `setup.sh` creates the conda environment from `environment.yml` in the ignored `pipelines/processing/work/env` directory (or `ENV_PREFIX`), syncs the root Python environment, and checks the required tools. It also downloads the miRBase names used to validate configs. The environment pins the rerun versions: FastQC 0.12.1, fastp 1.3.6, STAR 2.7.11b, featureCounts 2.1.1, R 4.5.3, DESeq2 1.50.2, and apeglm 1.32.0.
- `download.sh` downloads the genome, GTF, and reads. `--references-only` downloads just the genome and GTF. Completed downloads are kept, so the script can resume after an interruption; other arguments are passed to `src/download_reads.py` (for example, `--study GSE124411`).
- `run.sh` runs the configured comparisons, up to `RUN_WORKERS` at a time (default 2), and skips comparisons with a completed fit. It writes the fit list to `pipelines/processing/work/fits.tsv` and builds tables in `data/benchmark_data` (or `BENCHMARK_DATA_DIR`) once all 72 comparisons have completed. Arguments can select a subset of runs, but table generation requires all 72.

## Workflow

1. **Reads.** `pipelines/processing/src/download_reads.py` downloads the 738 FASTQ files (1,291.6 GB) listed in `pipelines/processing/configs/fastq_manifest.tsv` from ENA, and checks each file's md5 and size.
2. **References.** `download.sh` retrieves the soft-masked genome (`Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz`) and `Homo_sapiens.GRCh38.115.gtf.gz` from Ensembl release 115, verifies their checksums, and stores the decompressed files in `data/experiments/raw/refs/ensembl_v115/`. The genome is saved as `Homo_sapiens.GRCh38.dna.primary_assembly.fa`, the path used in the configs and by `funmirbench-experiments-download-examples`.
3. **DESeq2 fits.** Each comparison runs FastQC, fastp, STAR, featureCounts, and DESeq2. Use `run.sh` to process the full set, or run an individual comparison with `uv run funmirbench-experiments --config pipelines/processing/configs/<id>.yaml`. Configs list control and treated samples, read paths, and any sample-level batch terms (GSE115646 uses TR blocks).

   Some paired-end runs also have an unpaired FASTQ (reads without a mate), which the configs include: two samples of GSE210778, two of GSE158317 and all twelve of the two GSE78078 comparisons. DESeq2 writes `complete_analysis.rds` with the fitted object.
4. **Gene list.** `src/build_gene_list.py` lists the Ensembl 115 protein-coding genes with a MANE Select transcript: 19,226 genes, in `metadata/benchmark_genes.tsv.gz`.
5. **Benchmark tables.** `src/build_benchmark_data.py --fits fits.tsv --genes metadata/benchmark_genes.tsv.gz --out <dir>` runs `src/label_states.R` on every fit and copies `metadata/comparisons.tsv` into the output. Here, `fits.tsv` lists `experiment_id` and `fit_rds`.

Each DESeq2 run is saved under `pipelines/processing/work/runs/<timestamp>_<id>/` with `complete_analysis.rds` and `run_manifest.json`. Its DE table is written to `data/experiments/processed/<id>.tsv`.

To add one experiment without rebuilding the 72, see [ADD_EXPERIMENT.md](ADD_EXPERIMENT.md).

## Comparison configs

`pipelines/processing/configs/<id>.yaml` has:

- `dataset_id`, `gse`, `mirna_name`, and `experiment_type` identify the experiment;
- `source` records the input mode, reference paths, and thread counts. For read-based inputs, use `mode: reads`;
- `comparison` lists control and treated samples, with `sample_id`, `reads_1`, optional `reads_2`, and batch assignments when needed;
- `metadata` records the organism, cell line, treatment, and method.

For a count matrix, set `source.mode: count_matrix`, `count_matrix_path`, and `gene_id_column`; list sample columns under `comparison.control_columns` and `comparison.treated_columns`. This mode skips read alignment.

## Benchmark table format

Each comparison produces one table with 19,226 Ensembl 115 genes. The main columns are:

| Columns | Meaning |
|---|---|
| `logFC`, `lfcSE`, `stat`, `PValue`, `FDR` | DESeq2 Wald test against 0; FDR over the listed genes |
| `control_mean_normalized_count`, `treated_mean_normalized_count`, `*_nonzero_samples` | expression |
| `control_normalized_count_<sample>`, `treated_normalized_count_<sample>` | normalized count of each sample, one column per sample (so the number of columns differs between comparisons); the means above average them |
| `logFC_shrunk`, `lfcSE_shrunk`, `svalue`, `apeglm_conv` | apeglm shrinkage (prior from the full fit); `lfcSE_shrunk` is the posterior SD |
| `oriented_logFC`, `oriented_ci95_low/high`, `oriented_logFC_shrunk`, `oriented_ci95_low/high_shrunk` | sign flipped, so positive means "moved like a target" |
| `svalue_<d>` | FSOS ("false sign or small") s-value against d in the direction of repression, over the listed genes, for d in `d00`, `d01`, `d025`, `d05`, `d075`, `d10`, `d15`, `d20`, `d30` (0 to 3) |
| `state_<d>` | `responsive`, `no_response_above_margin` or `unresolved`, from the shrunken log2FC, for the same d values |
| `thr_pvalue_d05`, `thr_padj_d05`, `thr_pvalue_d10`, `thr_padj_d10` | one-sided DESeq2 test against margin d = 0.5 or 1.0 on the raw log2FC; FDR over the listed genes |
| `state_raw_d05`, `state_raw_d10` | the earlier raw states, kept as a check |
| `eligible`, `ineligible_reason`, `analysis_status` | why a gene has no usable estimate (not fitted, Cook's outlier, …) |

**States** (`state_<d>`, from the shrunken log2FC and its posterior SD).
- **responsive:** the oriented shrunken change is positive and the FSOS s-value against d is below 0.05. The per-gene probability is the normal tail of the oriented shrunken estimate below d, so that the 5% limit holds for the genes called responsive; the s-value is `apeglm::svalue` over the listed genes. A gene whose whole shrunken 95% interval is below d is no_response_above_margin even if its s-value passes, since the s-value averages over the genes ranked before it.
- **no_response_above_margin:** the whole shrunken 95% interval of the oriented change is below d. This also covers genes that moved the opposite way.
- **unresolved:** everything else.

**Raw states** (`state_raw_d05`, `state_raw_d10`): responsive if the raw threshold test passes at FDR < 0.05; no response if the whole raw 95% interval is below d.

Genes with zero counts in every sample keep a row, with no statistics.

## Normalization

Each comparison is normalized against its own control samples. DESeq2's size factors are computed against the geometric mean of the controls instead of all samples, and rescaled to a geometric mean of 1 over the controls. A control's normalized counts therefore do not depend on the treated samples, and a control shared by comparisons with the same controls has the same normalized counts in each.

## Metadata

These files are in `metadata/`.

- **`benchmark_genes.tsv.gz`:** the fixed gene list of step 4, built from the GTF alone.
- **`comparisons.tsv`:** one row per comparison.
  - design;
  - dependency group (comparisons that share controls or samples);
  - `action` (keep or flag) and `flags`.

  Human RNA-seq samples of the study that no comparison uses are listed with role `not_used`.
- **`predictions_info.tsv`, `benchmark_mirnas.tsv`:** the predictor list and miRNA sequences that the evaluation (`benchmark.yaml`) reads.

## Flags

The same rules apply to every comparison. Metrics are computed on the benchmark tables over the listed genes. Flags drop nothing; scores should be reported with and without flagged comparisons.

| Flag | Rule |
|---|---|
| `low_depth` (set aside) | median counted reads per sample below 100,000 |
| `global_shift` | median log2FC beyond ±0.15 |
| `most_genes_changed` | more than 75% of genes at FDR < 0.05 |
| `gc_shift` | Spearman ρ between log2FC and CDS GC beyond ±0.3 |
| `qc_samples_k_of_n` | k samples with STAR unique mapping below 60% or featureCounts assignment below 20% |

Result: 49 comparisons are kept and 23 flagged.
- QC samples: 7 comparisons, 5 of them GSE115646. GSE115646 libraries have heavy multimapping but still 7–10 million assigned reads per sample.
- GC shift: 16 comparisons.
- Most genes changed: GSE168484.

## Publication report feature comparison

The publication config can add 3′UTR length and conservation plots for predictor coverage groups to the main benchmark report. This step uses the Ensembl 115 annotation and hg38 phastCons100way track; see [Reproduce the publication analysis](../../README.md#reproduce-the-publication-analysis) for setup and run instructions. It is separate from rebuilding the RNA-seq benchmark tables.

## Provenance

The tables were built on 2026-11-05, using control-only size factors and counts from reruns with this pipeline's reads, references, and tool versions.
