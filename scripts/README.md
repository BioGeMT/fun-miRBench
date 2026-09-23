# Manuscript Figure Scripts

This folder contains scripts used to generate manuscript figure assets and the
tables that support them.

## Output Layout

Stable manuscript outputs are grouped by artifact type:

```text
manuscript_assets/figure2/
manuscript_assets/tables/
```

Figure-specific image assets belong under `manuscript_assets/figure<N>/`.
Reusable or supplementary TSV tables belong under `manuscript_assets/tables/`.

## Prepare Conservation Scores

Panel F uses gene-level mean phastCons100way scores over each gene's longest
protein-coding 3'UTR. If the conservation BigWig is not already available
locally, download it first:

```bash
uv run python scripts/download_phastcons100way.py
```

This downloads UCSC `hg38.phastCons100way.bw` to:

```text
data/resources/conservation/hg38.phastCons100way.bw
```

The BigWig is ignored by git because it is large.

Then compute the raw gene-level conservation table:

```bash
uv run python scripts/figure2_utr_conservation.py
```

By default this writes:

```text
manuscript_assets/tables/figure2F_utr_conservation_raw.tsv
```

The table contains `gene_id` and `utr3_mean_conservation`, plus supporting
transcript, scored-base, and 3'UTR-length columns. The main Figure 2 script
consumes this precomputed table so normal figure regeneration does not need to
re-read the multi-GB BigWig. This raw support table is generated locally and is
not tracked in git.

## Generate Figure 2

Run from the repository root:

```bash
uv run python scripts/figure2_coverage.py \
  --run-dir results/<YYYYMMDD_HHMMSS> \
  --panel all
```

This writes generated panel images and the combined six-panel figure to:

```text
manuscript_assets/figure2/
```

It writes the panel TSVs and supplementary Figure 2 tables to:

```text
manuscript_assets/tables/
```

The Figure 2 script reads the ground-truth FDR and effect thresholds from the
completed benchmark run's config snapshot:

```text
<run-dir>/benchmark_config.yaml
```

Regenerate the benchmark run with the current pipeline if that file is missing.

The combined figure is:

```text
manuscript_assets/figure2/figure2_combined.png
manuscript_assets/figure2/figure2_combined.svg
```

The main table outputs are:

```text
manuscript_assets/tables/figure2A_experiment_coverage.tsv
manuscript_assets/tables/figure2B_gene_set_coverage.tsv
manuscript_assets/tables/figure2C_positive_coverage.tsv
manuscript_assets/tables/figure2D_background_coverage.tsv
manuscript_assets/tables/figure2E_utr_length.tsv
manuscript_assets/tables/figure2F_utr_conservation.tsv
manuscript_assets/tables/figure2_supplementary_coverage_table.tsv
manuscript_assets/tables/figure2_supplement_gene_set_overlap_mirna_coverage.tsv
```

Only the combined Figure 2 image exports and compact manuscript-facing tables
are tracked. Individual panel images, the raw conservation support table, and
the verbose gene-set overlap support table are generated locally and ignored by
git.

The manuscript-focused supplement table
`figure2_supplement_gene_set_overlap_mirna_coverage.tsv` reports predictor
gene-set intersections and unions, plus the mean number of tested miRNAs scored
per covered gene for each predictor in the combination.

Panels:

- A: experiment coverage by predictor
- B: FGS/IGS gene-set coverage
- C: positive miRNA-gene pair coverage
- D: background miRNA-gene pair coverage
- E: longest protein-coding transcript 3'UTR length for IGS vs non-IGS genes
- F: mean 3'UTR conservation for IGS vs non-IGS genes

Panel E derives longest 3'UTR length from:

```text
data/resources/ensembl/Homo_sapiens.GRCh38.115.gtf.gz
```

For each gene, the script keeps one value: the longest 3'UTR among
protein-coding transcripts, after merging transcript-level 3'UTR intervals.

Panel F reads the precomputed gene-level conservation table:

```text
manuscript_assets/tables/figure2F_utr_conservation_raw.tsv
```

The table must contain:

```text
gene_id
utr3_mean_conservation
```
