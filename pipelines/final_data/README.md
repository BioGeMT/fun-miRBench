# v2 final benchmark data

From reads to the final per-gene tables for the 30 FuNmiRBench test comparisons.

## Steps

1. **Reads.** `pipelines/experiments/download_v2_reads.py` downloads the 260 FASTQ files (514.5 GB) listed in `pipelines/experiments/configs/v2/fastq_manifest.tsv` from ENA, and checks each file's md5 and size. For GSE199332, ENA holds both mates in one file as two positional blocks; the script splits them with `pipelines/experiments/split_mate_blocks.py`, which accepts the split only if true pairs overlap far more often than shuffled pairs.
2. **References.** Download `Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz` (https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/) and `Homo_sapiens.GRCh38.115.gtf.gz` (https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/) into `data/experiments/raw/refs/ensembl_v115/`, then decompress both. The configs point to these uncompressed files, which the reruns used; note that `funmirbench-experiments-download-examples` fetches the soft-masked genome instead.
3. **One DESeq2 fit per comparison.** Run `uv run funmirbench-experiments --config pipelines/experiments/configs/v2/<id>.yaml` for each comparison: FastQC, fastp, STAR, featureCounts, DESeq2. Configs record the control and treated samples, the read files, and a batch per sample where the design needs one:
   - GSE115646: TR blocks;
   - GSE51217: the two experimental sets.

   GSE210778 includes the unpaired reads of two samples. DESeq2 writes `complete_analysis.rds` with the fitted object.
4. **Gene list.** `build_gene_list.py` lists the Ensembl 115 protein-coding genes with a MANE Select transcript: 19,226 genes, in `metadata/v2/benchmark_genes.tsv.gz`.
5. **Final tables.** `build_final_data.py --fits fits.tsv --genes metadata/v2/benchmark_genes.tsv.gz --out <dir>` runs `label_states.R` on every fit and copies `metadata/v2/comparisons.tsv` into the output. Here, `fits.tsv` lists `experiment_id`, `perturbation` (OE or KO) and `fit_rds`.

## Final table (one per comparison, 19,226 rows)

| Columns | Meaning |
|---|---|
| `logFC`, `lfcSE`, `stat`, `PValue`, `FDR` | DESeq2 Wald test against 0; FDR over the listed genes |
| `control_mean_normalized_count`, `treated_mean_normalized_count`, `*_nonzero_samples` | expression |
| `logFC_shrunk`, `lfcSE_shrunk`, `svalue`, `apeglm_conv` | apeglm shrinkage (prior from the full fit) |
| `oriented_logFC`, `oriented_ci95_low/high`, `oriented_logFC_shrunk` | sign flipped for OE, so positive means "moved like a target" |
| `thr_pvalue_d05`, `thr_padj_d05`, `thr_pvalue_d10`, `thr_padj_d10` | one-sided DESeq2 test against margin d = 0.5 or 1.0; FDR over the listed genes |
| `state_d05`, `state_d10` | `responsive`, `no_response_above_margin` or `unresolved` |
| `eligible`, `ineligible_reason`, `analysis_status` | why a gene has no usable estimate (not fitted, Cook's outlier, …) |

**States.**
- **responsive:** the threshold test passes at FDR < 0.05.
- **no_response_above_margin:** the whole 95% interval of the oriented change is below d. This also covers genes that moved the opposite way.
- **unresolved:** everything else.

Genes with zero counts in every sample keep a row, with no statistics.

`metadata/v2/comparisons.tsv` and `metadata/v2/sample_sheets/` describe the comparisons and samples; see `metadata/v2/README.md`.
