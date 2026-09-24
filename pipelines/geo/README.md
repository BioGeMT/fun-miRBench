# GEO Download Pipeline

This is an **optional experiment reproduction workflow**. It is not required to run the published
benchmark; curated benchmark inputs are available through `uv run funmirbench-download-data`.

The default metadata source is:

```text
metadata/mirna_experiment_info.tsv
```

## Methodology

Control and perturbed-condition samples are **explicitly curated by the user**. The GEO pipeline
does not infer or replace these assignments.

For each experiment to reproduce, review the metadata row and fill:

```text
control_samples
condition_samples
```

with comma-separated GSM/SRR accessions. Rows where either field is empty are skipped with an INFO
message, so the registry can be refined incrementally.

The workflow is:

```text
metadata/mirna_experiment_info.tsv
        ↓
review control_samples / condition_samples
        ↓
geo_download.py
        ↓
FASTQs + pipelines/experiments/configs/<dataset_id>.yaml
        ↓
review generated YAML
        ↓
funmirbench/experiments_pipeline.py
```

## Environment

The GEO workflow has external dependencies that are intentionally kept in its Conda environment
(Biopython, SRA tools, pigz):

```bash
conda env create -f pipelines/geo/environment.yml
conda activate funmirbench-geo
```

Run this workflow with `python` inside that environment, not with the root `uv run` environment.

## Canonical metadata fields

The GEO workflow requires these columns:

```text
id
geo_accession
mirna_name
experiment_type
control_samples
condition_samples
```

The canonical registry also carries cell/tissue context, organism, method, PubMed ID, and the
processed DE-table path.

Example sample assignments:

```text
control_samples:   GSM3692987,GSM3692988,GSM3692989
condition_samples: GSM3692990,GSM3692991,GSM3692992
```

These assignments are part of the experiment-reproduction curation and should be checked against
the GEO record/publication before running the pipeline.

## Run one experiment

From the repository root:

```bash
python pipelines/geo/geo_download.py \
    --dataset-id GSE129076_OE_miR_450a_5p_1 \
    --entrez-email your@email.com
```

`--tsv` is optional and defaults to `metadata/mirna_experiment_info.tsv`.

You can also set:

```bash
export NCBI_ENTREZ_EMAIL="your@email.com"
```

For GEO/SRA mode, an Entrez email is required when a processable experiment is downloaded.

## Multiple experiments

Repeat `--dataset-id`:

```bash
python pipelines/geo/geo_download.py \
    --dataset-id GSE129076_OE_miR_450a_5p_1 \
    --dataset-id GSE129076_OE_miR_450a_5p_2 \
    --entrez-email your@email.com
```

If no dataset ID is provided, all rows with complete control/condition assignments are processed.

## Custom TSVs

A custom TSV can be supplied with:

```bash
python pipelines/geo/geo_download.py \
    --tsv path/to/custom_experiments.tsv \
    --entrez-email your@email.com
```

`pipelines/geo/input_experiments.tsv` is kept as a filled example of the same explicit-assignment
workflow.

Custom TSVs may also use:

- `raw_data_dir` for local FASTQ files;
- `count_matrix_path` and `gene_id_column` for a pre-existing count matrix.

Both modes still use explicit `control_samples` and `condition_samples`.

## Optional metadata helper

`fetch_geo_metadata.py` remains an optional helper for inspecting/prefilling GEO metadata.
Any suggested biological metadata or sample grouping must be reviewed before being used in the
canonical registry.

## Outputs

| Path | Description |
|---|---|
| `data/experiments/raw/<GSE>/*.fastq.gz` | downloaded FASTQ files |
| `data/experiments/raw/<GSE>/manifest.json` | sample-to-run/file manifest |
| `pipelines/experiments/configs/<dataset_id>.yaml` | generated experiment-processing config |

Always review the generated YAML before running `funmirbench/experiments_pipeline.py`.
