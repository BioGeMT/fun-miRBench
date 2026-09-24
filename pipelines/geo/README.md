# GEO Download Pipeline

This is an **optional experiment acquisition workflow** for extending or reproducing experiment
processing in fun-miRBench. It is not required to run the published benchmark; the curated
publication inputs are available through `uv run funmirbench-download-data`.

The default source of experiment identity is the canonical registry:

```text
metadata/mirna_experiment_info.tsv
```

That registry contains the stable experiment ID, GEO accession, miRNA, perturbation, cell/tissue
context, publication metadata, and processed DE-table path. It intentionally does **not** store
raw control/condition sample assignments.

## Workflow

```text
metadata/mirna_experiment_info.tsv
        ↓
geo_download.py discovers GEO samples
        ↓
pipelines/geo/sample_assignments/<dataset_id>.tsv
        ↓
review group = control / condition / exclude
        ↓
rerun geo_download.py
        ↓
FASTQs + pipelines/experiments/configs/<dataset_id>.yaml
        ↓
review YAML
        ↓
run funmirbench/experiments_pipeline.py
```

The sample-assignment review step is deliberate. GEO sample titles and characteristics are used to
suggest groups, but those suggestions are not treated as ground truth.

## 1. Discover samples from the canonical registry

Process one experiment:

```bash
python pipelines/geo/geo_download.py \
    --dataset-id GSE115646_OE_miR-18a-5p_HepG2
```

The `--tsv` option is optional. By default the script reads:

```text
metadata/mirna_experiment_info.tsv
```

On the first run, if no reviewed sample assignment exists, the script fetches GEO SOFT metadata,
classifies samples heuristically, and writes:

```text
pipelines/geo/sample_assignments/<dataset_id>.tsv
```

No FASTQ files are downloaded at this stage.

The assignment file contains:

| Column | Meaning |
|---|---|
| `sample_id` | GEO sample accession, usually GSM |
| `title` | GEO sample title |
| `suggested_group` | heuristic suggestion: control, condition, or uncertain |
| `group` | review field that must be set manually to `control`, `condition`, or `exclude` |
| `control_score` | heuristic control score |
| `condition_score` | heuristic condition score |
| `organism` | sample organism from GEO metadata |
| `cell_line` | inferred cell-line context when available |
| `tissue` | inferred tissue context when available |

## 2. Review the sample assignment

Open the generated TSV and fill the `group` column for every row:

```text
control
condition
exclude
```

The workflow will not continue while any sample has an empty or unsupported group.

At least one control and one condition sample are required.

## 3. Download FASTQs and generate the experiment config

After reviewing the sample assignment, rerun the same command and provide an Entrez email:

```bash
python pipelines/geo/geo_download.py \
    --dataset-id GSE115646_OE_miR-18a-5p_HepG2 \
    --entrez-email your@email.com
```

You can also set:

```bash
export NCBI_ENTREZ_EMAIL="your@email.com"
```

The workflow then:

1. resolves reviewed GSM samples to SRA runs;
2. downloads the FASTQs;
3. writes `data/experiments/raw/<GSE>/manifest.json`;
4. generates `pipelines/experiments/configs/<dataset_id>.yaml`.

The generated YAML must still be reviewed before running the RNA-seq processing pipeline.

## Selecting experiments

Repeat `--dataset-id` to process multiple experiments:

```bash
python pipelines/geo/geo_download.py \
    --dataset-id GSE115646_OE_miR-18a-5p_HepG2 \
    --dataset-id GSE169128_KO_miR_3662
```

If no `--dataset-id` is supplied, all rows in the selected metadata TSV are considered. For normal
manual reproduction work, selecting experiments explicitly is recommended.

## Custom TSVs

A custom TSV can still be supplied:

```bash
python pipelines/geo/geo_download.py \
    --tsv path/to/custom_experiments.tsv \
    --entrez-email your@email.com
```

The required identity columns are:

```text
id
geo_accession
mirna_name
experiment_type
```

Custom TSVs may additionally contain `control_samples` and `condition_samples`. If both are
present for an experiment, the review-file discovery step is skipped and those explicit assignments
are used directly.

The older `pipelines/geo/input_experiments.tsv` format remains usable as a custom-input workflow,
but it is no longer the canonical source of experiment identity.

## Local reads and count-matrix modes

Custom TSVs may also provide:

- `raw_data_dir` plus explicit `control_samples` / `condition_samples` for local FASTQ files;
- `count_matrix_path`, `gene_id_column`, and explicit sample assignments for count-matrix mode.

These modes require explicit sample assignments because GEO discovery is not involved.

## Outputs

| Path | Description |
|---|---|
| `pipelines/geo/sample_assignments/<dataset_id>.tsv` | review-required control/condition assignment |
| `data/experiments/raw/<GSE>/*.fastq.gz` | downloaded FASTQ files |
| `data/experiments/raw/<GSE>/manifest.json` | sample-to-run/file manifest |
| `pipelines/experiments/configs/<dataset_id>.yaml` | generated experiment-processing config |

## Notes

- The canonical registry remains focused on experiment identity and benchmark provenance.
- Sample grouping stays in a separate review artifact because it is a reproduction decision, not a
  benchmark identity field.
- GEO-derived group suggestions are heuristic and must be reviewed.
- `--entrez-email` is needed only after sample assignments are reviewed and SRA resolution begins.
- If any run download fails, config generation is aborted for that experiment so incomplete inputs
  are not propagated.

