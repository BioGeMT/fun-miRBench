# GEO Download Pipeline

This pipeline is the first stage of the FuNmiRBench experiment processing workflow. It downloads raw FASTQ files from GEO/SRA (or locates local FASTQ files) or verifies a pre-existing count matrix and automatically generates YAML configuration files for the RNA-seq pipeline (`funmirbench/experiments_pipeline.py`).

---

## Overview of the full pipeline

```
0. (Optional) Run fetch_geo_metadata.py   →   auto-fills a row in input_experiments.tsv
          ↓
1. Review / fill in pipelines/geo/input_experiments.tsv
          ↓
2. Run geo_download.py   →   pipelines/experiments/configs/{dataset_id}.yaml  (auto-generated)
          ↓
3. Review and fine-tune the generated YAML config  ← important, see note below
          ↓
4. Run funmirbench/experiments_pipeline.py --config pipelines/experiments/configs/{dataset_id}.yaml
```

---

## Step 0 (Optional) - Auto-fill the TSV from GEO

`fetch_geo_metadata.py` fetches metadata for a GEO series from the NCBI SOFT API
and appends a pre-filled row to `input_experiments.tsv`. It can optionally use the
**Gemini Flash LLM** to fill in fields that require biological interpretation
(miRNA name, experiment type, treatment description).

### Option A - Rule-based only

Uses keyword matching and majority voting across sample metadata. No API key needed.

```bash
conda activate funmirbench-geo

python pipelines/geo/fetch_geo_metadata.py --gse-url GSE93717
```

### Option B - With Gemini Flash LLM (recommended)

Adds LLM-assisted extraction of `mirna_name`, `experiment_type`, `treatment`,
`tested_cell_line`, `tissue`, and `organism`. Also validates the suggested miRNA
name against a local miRBase cache.

#### Get a free Gemini API key

1. Go to **[aistudio.google.com](https://aistudio.google.com)** and sign in with a Google account.
2. Click **"Get API key"** → **"Create API key"**.
3. Copy the key (starts with `AIza...`).

No billing required — the free tier is sufficient for this use case.

#### Run with the LLM

Set the Gemini API key as an environment variable first. This is preferred over
passing the key directly on the command line, because command-line arguments can
be saved in shell history or exposed in process listings.

```bash
export GEMINI_API_KEY="AIzaXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"

conda activate funmirbench-geo

python pipelines/geo/fetch_geo_metadata.py \
    --gse-url GSE93717 \
    --llm gemini
```

You can still pass `--gemini-key` explicitly if needed, but avoid doing so on
shared systems.

## Step 1 - Fill in the input TSV

Edit `pipelines/geo/input_experiments.tsv`. For each experiment you want to process,
fill in the relevant columns depending on the mode.

> **All columns are mandatory.** If a value is not available for a given field
> (e.g. `article_pubmed_id`), fill it with `NA` rather than leaving it empty.

| Column | Description |
|---|---|
| `control_samples` | Comma-separated list of control sample identifiers (SRR/GSM accessions, file base-names, or count matrix column names) |
| `condition_samples` | Comma-separated list of condition sample identifiers (same format as above) |
| `raw_data_dir` | *Reads mode only.* Path to a local directory with pre-existing FASTQ files. Leave empty to download from GEO. |
| `count_matrix_path` | *Count matrix mode only.* Path to the count matrix file. If set, count matrix mode is used. See [example_count_matrix.csv](example_count_matrix.csv) for the expected format. |
| `gene_id_column` | *Count matrix mode only.* Column name for gene IDs in the matrix (e.g. `ENSEMBL`). Required when `count_matrix_path` is set. |

Rows with empty `control_samples` or `condition_samples` are skipped automatically,
so you can fill in the TSV incrementally.

> The file `pipelines/geo/input_experiments.tsv` already contains a real filled-in row
> as a reference. The examples below are illustrative — **replace them with your own data**.

### GEO mode - download FASTQ files from GEO/SRA

Leave `raw_data_dir` and `count_matrix_path` empty. The pipeline auto-detects the
accession type from the prefix and routes accordingly.

**Option A — SRR accessions (recommended)**

Provide SRR run accessions directly. The pipeline fetches the library layout (single/paired)
from NCBI and downloads immediately — no GSM resolution step needed.

```
id:                  GSE129076_OE_miR_21_5p
mirna_name:          hsa-miR-21-5p
article_pubmed_id:   https://pubmed.ncbi.nlm.nih.gov/34685526
organism:            Homo sapiens
tested_cell_line:    A549
treatment:           Overexpression of miR-21-5p
tissue:              Lung
method:              RNA-seq
experiment_type:     Overexpression
gse_url:             https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129076
raw_data_dir:        (empty)
control_samples:     SRR8816234,SRR8816235,SRR8816236
condition_samples:   SRR8816237,SRR8816238,SRR8816239
count_matrix_path:   (empty)
gene_id_column:      (empty)
```

**Option B — GSM accessions**

Provide GSM sample accessions. The pipeline resolves each GSM to its SRR run(s) via
the NCBI API before downloading. Useful when you want to copy accessions directly from
the GEO page, or when one sample has multiple runs that you prefer not to list individually.
Multiple runs belonging to one GSM remain one biological sample in the generated
config and are combined before QC and trimming.

```
id:                  GSE129076_OE_miR_21_5p
mirna_name:          hsa-miR-21-5p
article_pubmed_id:   https://pubmed.ncbi.nlm.nih.gov/34685526
organism:            Homo sapiens
tested_cell_line:    A549
treatment:           Overexpression of miR-21-5p
tissue:              Lung
method:              RNA-seq
experiment_type:     Overexpression
gse_url:             https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129076
raw_data_dir:        (empty)
control_samples:     GSM6437108,GSM6437109,GSM6437110
condition_samples:   GSM6437113,GSM6437114,GSM6437115
count_matrix_path:   (empty)
gene_id_column:      (empty)
```

### Local reads mode — use pre-existing FASTQ files

Set `raw_data_dir` to the directory containing your FASTQ files and provide
**sample base-names** (without extension) as identifiers.

```
id:                  GSE999999_KO_miR_155_5p
mirna_name:          hsa-miR-155-5p
article_pubmed_id:   NA
organism:            Homo sapiens
tested_cell_line:    HaCaT
treatment:           Knockdown of miR-155-5p
tissue:              Skin
method:              RNA-seq
experiment_type:     Knockdown
gse_url:             https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE999999
raw_data_dir:        /path/to/your/fastq/files
control_samples:     ctrl_rep1,ctrl_rep2,ctrl_rep3
condition_samples:   treated_rep1,treated_rep2,treated_rep3
count_matrix_path:   (empty)
gene_id_column:      (empty)
```

The script looks for:
- Single-end: `{raw_data_dir}/{name}.fastq.gz`
- Paired-end (auto-detected): `{raw_data_dir}/{name}_1.fastq.gz` + `{raw_data_dir}/{name}_2.fastq.gz`

### Count matrix mode - use a pre-existing count matrix

> The matrix must be a CSV with gene annotation columns (`ENSEMBL`, `name`, `type`, `chr`, `start`, `end`, `str`) followed by per-sample count columns. See [example_count_matrix.csv](example_count_matrix.csv) for the expected format.

Set `count_matrix_path` to the matrix file and provide the **column names** from
the matrix as sample identifiers.

```
id:                  GSE253003_OE_miR_323_3p
mirna_name:          hsa-miR-323-3p
article_pubmed_id:   https://pubmed.ncbi.nlm.nih.gov/12345678
organism:            Homo sapiens
tested_cell_line:    HSVSMC
treatment:           Overexpression of miR-323-3p
tissue:              Vascular
method:              RNA-seq
experiment_type:     Overexpression
gse_url:             https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253003
raw_data_dir:        (empty)
control_samples:     HSVSMC_IP_miRCTRL_1,HSVSMC_IP_miRCTRL_2,HSVSMC_IP_miRCTRL_3
condition_samples:   HSVSMC_IP_MIR323_1,HSVSMC_IP_MIR323_2,HSVSMC_IP_MIR323_3
count_matrix_path:   data/experiments/raw/GSE253003/GSE253003_Count.csv.gz
gene_id_column:      ENSEMBL
```

The script only verifies that the matrix file exists and then generates the YAML config.
No downloading takes place.

---

## Step 2 - Set up the environment

```bash
conda env create -f pipelines/geo/environment.yml
conda activate funmirbench-geo
```

---

## Step 3 - Run the pipeline

Run from the **repo root**:

```bash
python pipelines/geo/geo_download.py \
    --tsv pipelines/geo/input_experiments.tsv \
    --entrez-email your@email.com
```

FASTQs are always saved to `data/experiments/raw/` and YAML configs to
`pipelines/experiments/configs/`.

Key options:

| Flag | Default | Description |
|---|---|---|
| `--tsv` | *(required)* | Path to the input TSV |
| `--entrez-email` | *(required)* | Email address for the NCBI Entrez API |
| `--threads` | `4` | Threads for `fasterq-dump` |

---

## Step 4 - Review the generated YAML config

> **Before running the RNA-seq pipeline, always open and review the generated YAML.**

The YAML is auto-generated from the TSV and uses sensible defaults, but you should
verify and adjust the following before proceeding:

- **Genome reference paths** — the defaults point to Ensembl v115 (GRCh38). If your
  experiment uses a different genome or Ensembl version, update `genome_fasta_path`
  and `gtf_path` accordingly.
- **Thread counts** — `fastqc_threads`, `fastp_threads`, `star_threads`, and
  `featurecounts_threads` are set to fixed defaults. Adjust them to match the
  resources available on your machine.
- **Sample identifiers** — check that `sample_id` values and file paths look correct.
  Multiple SRR paths for one GSM should appear under one biological sample.
- **Metadata fields** — organism, cell line, tissue, and treatment are copied from
  the TSV; verify they are accurate before archiving results.

The generated YAML is at:
```
pipelines/experiments/configs/{dataset_id}.yaml
```

---

## Step 5 - Run the RNA-seq pipeline

Once the YAML config is reviewed and ready, refer to the [main project README](../../README.md)
for instructions on running the RNA-seq pipeline (`funmirbench/experiments_pipeline.py`).

---

## Outputs

| Path | Description |
|---|---|
| `data/experiments/raw/{GSE}/{SRR}.fastq.gz` | Downloaded FASTQ files (GEO mode) |
| `data/experiments/raw/{GSE}/manifest.json` | Maps sample → SRR → file paths and group (GEO mode) |
| `pipelines/experiments/configs/{dataset_id}.yaml` | Auto-generated RNA-seq pipeline config |

---

## Notes

- `--entrez-email` is required for GEO mode (NCBI policy); it is not used in local or count matrix mode but is still a required argument.
- If any SRR download fails, the entire experiment is skipped — no YAML config is written for it to avoid referencing missing files.
