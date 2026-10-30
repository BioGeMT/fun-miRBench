# Final benchmark data

From reads to the final per-gene tables for the 21 FuNmiRBench OE test comparisons.

## Run

```bash
pipelines/final_data/setup.sh
pipelines/final_data/download.sh
pipelines/final_data/run.sh
```

- `setup.sh` creates the conda environment of `environment.yml` in the ignored `pipelines/final_data/work/env` (or `ENV_PREFIX`), syncs the Python environment and checks the tools. `environment.yml` pins the tool versions of the reruns (FastQC 0.12.1, fastp 1.3.6, STAR 2.7.11b, featureCounts 2.1.1, R 4.5.3, DESeq2 1.50.2) and apeglm 1.32.0.
- `download.sh` runs steps 1 and 2 and gets the miRBase names the pipeline checks configs against. It keeps complete files, so it can be rerun after an interruption; its arguments go to `src/download_reads.py` (for example `--study GSE124411`).
- `run.sh` runs step 3 for every comparison, `RUN_WORKERS` at a time (default 2), and skips comparisons with a finished run. It then writes the fit list to `pipelines/final_data/work/fits.tsv` and runs step 5 into `results/final_data` (or `FINAL_DIR`). Its arguments limit the runs to those comparisons; the final tables need all 21.

## Steps

1. **Reads.** `pipelines/final_data/src/download_reads.py` downloads the 178 FASTQ files (317.7 GB) listed in `pipelines/final_data/configs/fastq_manifest.tsv` from ENA, and checks each file's md5 and size.
2. **References.** Download the soft-masked genome `Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz` (https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/) and `Homo_sapiens.GRCh38.115.gtf.gz` (https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/) into `data/experiments/raw/refs/ensembl_v115/`, decompress both, and save the genome as `Homo_sapiens.GRCh38.dna.primary_assembly.fa`, the name the configs point to. This is the genome the reruns used, and the one `funmirbench-experiments-download-examples` fetches.
3. **One DESeq2 fit per comparison.** Run `uv run funmirbench-experiments --config pipelines/final_data/configs/<id>.yaml` for each comparison: FastQC, fastp, STAR, featureCounts, DESeq2. Configs record the control and treated samples, the read files, and a batch per sample where the design needs one (GSE115646: TR blocks).

   GSE210778 includes the unpaired reads of two samples. DESeq2 writes `complete_analysis.rds` with the fitted object.
4. **Gene list.** `src/build_gene_list.py` lists the Ensembl 115 protein-coding genes with a MANE Select transcript: 19,226 genes, in `metadata/benchmark_genes.tsv.gz`.
5. **Final tables.** `src/build_final_data.py --fits fits.tsv --genes metadata/benchmark_genes.tsv.gz --out <dir>` runs `src/label_states.R` on every fit and copies `metadata/comparisons.tsv` into the output. Here, `fits.tsv` lists `experiment_id` and `fit_rds`.

## Final table (one per comparison, 19,226 rows)

| Columns | Meaning |
|---|---|
| `logFC`, `lfcSE`, `stat`, `PValue`, `FDR` | DESeq2 Wald test against 0; FDR over the listed genes |
| `control_mean_normalized_count`, `treated_mean_normalized_count`, `*_nonzero_samples` | expression |
| `logFC_shrunk`, `lfcSE_shrunk`, `svalue`, `apeglm_conv` | apeglm shrinkage (prior from the full fit); `lfcSE_shrunk` is the posterior SD |
| `oriented_logFC`, `oriented_ci95_low/high`, `oriented_logFC_shrunk`, `oriented_ci95_low/high_shrunk` | sign flipped, so positive means "moved like a target" |
| `svalue_<d>` | FSOS ("false sign or small") s-value against d, over the listed genes, for d in `d00`, `d01`, `d025`, `d05`, `d075`, `d10`, `d15`, `d20`, `d30` (0 to 3) |
| `state_<d>` | `responsive`, `no_response_above_margin` or `unresolved`, from the shrunken log2FC, for the same d values |
| `thr_pvalue_d05`, `thr_padj_d05`, `thr_pvalue_d10`, `thr_padj_d10` | one-sided DESeq2 test against margin d = 0.5 or 1.0 on the raw log2FC; FDR over the listed genes |
| `state_raw_d05`, `state_raw_d10` | the earlier raw states, kept as a check |
| `eligible`, `ineligible_reason`, `analysis_status` | why a gene has no usable estimate (not fitted, Cook's outlier, …) |

**States** (`state_<d>`, from the shrunken log2FC and its posterior SD).
- **responsive:** the oriented shrunken change is positive and the FSOS s-value against d is below 0.05. The per-gene probability is the normal tail of the shrunken estimate below d, which is what `lfcShrink(lfcThreshold = d)` computes; the s-value is then taken over the listed genes.
- **no_response_above_margin:** the whole shrunken 95% interval of the oriented change is below d. This also covers genes that moved the opposite way.
- **unresolved:** everything else.

**Raw states** (`state_raw_d05`, `state_raw_d10`): responsive if the raw threshold test passes at FDR < 0.05; no response if the whole raw 95% interval is below d.

Genes with zero counts in every sample keep a row, with no statistics.

`metadata/comparisons.tsv` and `metadata/sample_sheets/` describe the comparisons and samples; see `metadata/README.md`.
