# Add an experiment

Add a miRNA overexpression experiment from raw reads or from a count matrix. Both go through the same DESeq2 fit and labels as the published experiments.

```bash
pipelines/processing/setup.sh
pipelines/processing/download.sh --references-only # required for raw-read input
uv run funmirbench-download-data                    # to include the published experiments
uv run funmirbench-add --config my_experiment.yaml
uv run funmirbench --config benchmark.yaml
```

Run these commands from the repository root. `setup.sh` prepares the conda environment for R and the aligners and syncs the Python environment. The references-only download is needed for raw reads; count-matrix input does not use genome alignment. Download the published data before adding your experiment if you want to score it alongside the published experiments. `funmirbench-add` adds the new table and metadata to `data/benchmark_data/`; the benchmark config then includes it with the selected experiments. Paths in the experiment config are relative to the repository root unless they are absolute.

With raw reads:

```yaml
dataset_id: my_experiment
mirna_name: hsa-miR-375-3p
experiment_type: Overexpression
source:
  mode: reads
  genome_fasta_path: data/experiments/raw/refs/ensembl_v115/Homo_sapiens.GRCh38.dna.primary_assembly.fa
  gtf_path: data/experiments/raw/refs/ensembl_v115/Homo_sapiens.GRCh38.115.gtf
comparison:
  control_samples:
  - {sample_id: ctrl_1, reads_1: reads/ctrl_1_R1.fastq.gz, reads_2: reads/ctrl_1_R2.fastq.gz}
  - {sample_id: ctrl_2, reads_1: reads/ctrl_2_R1.fastq.gz, reads_2: reads/ctrl_2_R2.fastq.gz}
  treated_samples:
  - {sample_id: oe_1, reads_1: reads/oe_1_R1.fastq.gz, reads_2: reads/oe_1_R2.fastq.gz}
  - {sample_id: oe_2, reads_1: reads/oe_2_R1.fastq.gz, reads_2: reads/oe_2_R2.fastq.gz}
metadata:
  tested_cell_line: HeLa
```

For single-end reads, leave out `reads_2`. For a sample sequenced in several runs, give `reads_1` and `reads_2` as lists.

With a count matrix, one row per gene and one column per sample:

```yaml
dataset_id: my_experiment
mirna_name: hsa-miR-375-3p
experiment_type: Overexpression
source:
  mode: count_matrix
  count_matrix_path: counts.tsv
  gene_id_column: gene_id
comparison:
  control_columns: [ctrl_1, ctrl_2]
  treated_columns: [oe_1, oe_2]
metadata:
  tested_cell_line: HeLa
```

- **Counts:** provide raw integer counts, one row per Ensembl gene ID without a version suffix. Do not use normalized counts or other fractional values. Counts made with another aligner or annotation are not processed exactly like the published experiments; `source_mode` in `comparisons.tsv` records the input route for each added experiment.
- **Samples:** at least two per side. Add `batch_by_sample` under `comparison`, with the batch of each sample, for a `~ batch + condition` design.
- **Shared samples:** experiments that share samples should get the same `dependency_group` in the config, so that the bootstrap treats them as one group.
- **miRNA:** a mature human miRNA name from miRBase 22.1. A miRNA new to the benchmark gets its sequence added to `metadata/benchmark_mirnas.tsv`, which also adds it to the decoys of the other experiments. A predictor without scores for it ranks every gene last.
- **Preserving added experiments:** do not run `funmirbench-download-data --force` after adding experiments. That command replaces the downloaded benchmark tables and removes locally added rows.
