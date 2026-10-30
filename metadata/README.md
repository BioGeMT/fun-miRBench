# Metadata of the 21 OE test comparisons

- **`benchmark_genes.tsv.gz`:** the fixed gene list. It has 19,226 Ensembl 115 protein-coding genes with a MANE Select transcript, built by `pipelines/benchmark_data/src/build_gene_list.py` from the GTF alone.
- **`comparisons.tsv`:** one row per comparison.
  - design, sample numbers;
  - dependency group (comparisons that share controls or samples) and exact-seed group;
  - QC metrics, `action` (keep or flag) and `flags`.
- **`sample_sheets/<GSE>.tsv`:** one row per sample per comparison.
  - role, design, batch, donor, clone, dose and source block;
  - read files (repository paths) and GEO title and characteristics;
  - the number of comparisons sharing the sample.

  Human RNA-seq samples of the study that no comparison uses are listed with role `not_used`.
- **`mirna_experiment_info.tsv`, `predictions_info.tsv`:** the experiment and predictor lists that the evaluation (`benchmark.yaml`) reads.

## Flags

The same rules apply to every comparison. Metrics are computed on the benchmark tables over listed genes with a control mean of at least 10 normalised reads. Flags drop nothing; scores should be reported with and without flagged comparisons.

| Flag | Rule |
|---|---|
| `low_depth` (set aside) | median counted reads per sample below 100,000 |
| `global_shift` | median log2FC beyond ±0.15 |
| `most_genes_changed` | more than 75% of genes at FDR < 0.05 |
| `gc_shift` | Spearman ρ between log2FC and CDS GC beyond ±0.3 |
| `qc_samples_k_of_n` | k samples with STAR unique mapping below 60% or featureCounts assignment below 20% |

Result: 9 comparisons are kept and 12 flagged.
- QC samples: all 6 GSE115646 and both GSE129076. GSE115646 libraries have heavy multimapping but still 7–10 million assigned reads per sample.
- GC shift: 5 comparisons.
- Global shift and most genes changed: GSE168484.

## Provenance

The tables were built on 2026-10-17 from the fresh reruns of 2026-10-11 and the input QC review. `seed_signal_audit` is the October 2026 audit's seed-response category; it is descriptive and was computed on the old labels.
