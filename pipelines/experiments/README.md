# Experiment Ingestion Pipeline

This pipeline turns raw experiment data (a pre-existing count matrix or FASTQs) into a
benchmark-ready differential-expression (DE) table under `data/experiments/processed/`, using
DESeq2 (`funmirbench/experiments_pipeline.py`). It is typically the second stage of the full
experiment workflow, after the [GEO download pipeline](../geo/README.md) has produced a YAML
config — though a config can also be written by hand.

---

## Overview of the full pipeline

```
0. (Optional) Run the GEO download pipeline   →   auto-generates a YAML config
          ↓
1. Review / write pipelines/experiments/configs/{dataset_id}.yaml
          ↓
2. Run funmirbench-experiments --config ...   →   data/experiments/processed/{dataset_id}.tsv
          ↓
3. Run funmirbench-sync-metadata   →   upserts row into metadata/mirna_experiment_info.tsv
```

---

## Step 0 (Optional) - Acquire data from GEO

Raw FASTQs or a count matrix, plus a ready-to-review YAML config, can be generated automatically
from a GEO accession. See the [GEO download pipeline](../geo/README.md) for the full workflow.

## Step 1 - Write or review the YAML config

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

Tracked example configs:

- `pipelines/experiments/configs/gse253003.count_matrix.example.yaml`
- `pipelines/experiments/configs/gse93717.reads.example.yaml`

Reads configs can either:

- use local `reads_1` and optional `reads_2`
- use local `genome_fasta_path` and `gtf_path`

---

## Step 2 - Set up the environment

```bash
conda env create -f pipelines/experiments/environment.yml
conda activate funmirbench-experiments
```

That environment also includes `uv`, so `uv run ...` keeps working after activation, and provides
`fastqc`, `fastp`, `STAR`, `featureCounts`, and `Rscript` on `PATH`.

---

## Step 3 - Run the pipeline

RNA-seq pipeline summary:

- `count_matrix` mode validates the configured count matrix and runs DESeq2 directly on the selected
  control and treated columns.
- `reads` mode runs read QC/trimming, aligns reads to the configured reference with STAR, counts
  genes with featureCounts, then runs the same DESeq2 step.
- DE outputs preserve the original method-level `PValue` and `FDR` values. Benchmarking derives
  `benchmark_FDR` and `plot_FDR` internally without modifying those RNA-seq outputs.

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

The reads example uses a reproduced dataset id, `GSE93717_OE_miR_941_deseq2`, so syncing it creates
a separate variant instead of overwriting the curated `GSE93717_OE_miR_941` registry row.

### Output DE table schema

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
- `control_mean_normalized_count`: mean of the DESeq2-normalized counts across the control samples
- `PValue`: optional raw p-value used to classify rows with missing adjusted p-values

Common native headers are normalized when benchmark tables are read. For example, `log2FoldChange`
maps to `logFC`, `padj`, `adj.P.Val`, and `qvalue` map to `FDR`, and `pvalue` or `P.Value` map to
`PValue`.

---

## Step 4 - Sync into the experiment registry

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

---

## Outputs

| Path | Description |
|---|---|
| `data/experiments/processed/<dataset_id>.tsv` | Final DE table |
| `pipelines/experiments/runs/<timestamp>_<dataset_id>/candidate_metadata.tsv` | Registry row candidate, synced via `funmirbench-sync-metadata` |
| `pipelines/experiments/runs/<timestamp>_<dataset_id>/run_manifest.json` | Run manifest (commands, inputs, timestamps) |
