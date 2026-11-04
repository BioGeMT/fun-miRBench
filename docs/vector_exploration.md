# Exploratory vector analysis from Zenodo

This workflow downloads a pinned published data snapshot, runs the vector analyses and builds
standalone English HTML and PDF reports. It does not run RNA-seq processing or change the
main evaluation command/report. Its data, matrices and results use separate folders.

## Full run

From the repository root:

```bash
git switch feature/response-vector-analysis
git pull --ff-only
uv sync

uv run funmirbench-vector-explore --config benchmark_vectors_zenodo.yaml --published-only --details 21
```

By default this selects the 21 OE comparisons, five predictors and 62 registry miRNAs on
the shared 19,226-gene universe. Filters in the exploration config can select a subset.
Flagged comparisons are retained and identified; KO and set-aside comparisons are skipped.

The data sources are:

| Input | Published source |
|---|---|
| Experimental tables and expanded-family TargetScan | Zenodo record [23839843](https://zenodo.org/records/23839843) |
| miRDB, microT-CNN, miRBind2 and miRAW predictions | Standardized predictions ZIP in record [21671476](https://zenodo.org/records/21671476) |
| Selected miRNA sequences, tool names and score directions | Repository metadata TSVs copied into each run |

The first download is approximately 2.22 GB for all five predictors: 85.7 MB of tables,
120.4 MB of TargetScan predictions and a 2.01 GB prediction archive. No raw reads or
alignments are downloaded. Archives are checksum-verified, retained and reused.

Without `--published-only`, existing files referenced by the predictor registry are reused
when present. Their provenance is explicitly labelled local and not verified against Zenodo.
Custom predictors need a local standardized TSV and registry entry.

For a smaller first run, use only TargetScan (about 206 MB of published downloads):

```bash
uv run funmirbench-vector-explore --predictors targetscan_family_expanded --published-only
```

## Background run

```bash
nohup uv run funmirbench-vector-explore \
  --config benchmark_vectors_zenodo.yaml --published-only --details 21 \
  > vector_explore.log 2>&1 < /dev/null &

tail -f vector_explore.log
```

This exploration is independent of any RNA-seq jobs. Closing `tail` does not stop it.

## Outputs

Each run prints its output directory:

```text
results/vector_explore/<timestamp>/response_vectors/
```

Open **REPORT.html** for the complete interactive document: overview, figures, expandable
per-experiment tables, missing-value reasons, methods and limitations. It embeds its images,
so it can be copied as one file. **REPORT.pdf** contains methods, descriptive overview tables
and the same figures for sharing.

| File | Contents |
|---|---|
| `overview.tsv` | Descriptive medians and valid experiment counts, on all and unflagged comparisons |
| `specificity.tsv` | Cognate alignment, coverage, decoy ranks/fractions/z-scores, panel counts and NA reasons |
| `candidate_alignment.tsv` | Every experiment–predictor–candidate miRNA alignment |
| `expression_bins.tsv` | Cognate correlations per expression bin, usable gene counts and NA reasons |
| `comparisons.tsv`, `skipped_comparisons.tsv` | Included and skipped comparisons with metadata |
| `provenance.json`, `analysis.json` | Data records/checksums, prediction sources and cache identity |
| `timings.json` | Measured runtime and analysis-cache hit status |
| `plots/` | Reusable PNG figures and their captions |

Per-gene vector exports are optional because they increase writing time and output size:

```bash
uv run funmirbench-vector-explore --published-only --export-vectors
```

## Report questions and plots

| Question | Figure |
|---|---|
| Do predictions agree with repression after accounting for baseline expression? | Paired overall and within-bin alignments; points are experiments, black symbols are medians |
| Is the agreement miRNA specific? | Cognate alignment versus the fraction of better decoys, overall and within-bin |
| Which experiments carry the signal? | Alignment and specificity heatmaps with flagged rows and grey NA cells |
| How much information supports the results? | Coverage among usable-response genes and comparable decoy counts |
| Does performance depend on baseline expression? | Median per-bin alignment profiles, backed by the full bin table |
| Does the cognate stand out from other miRNAs in an individual experiment? | Decoy-correlation dot plots with the cognate marked by a black diamond |

The default draws detailed decoy plots for the first three comparisons in metadata order,
not the three best results. All comparisons and candidates remain in the numerical tables.
Use `--details 21` for every current experiment or `--details 0` for overview figures only.
No violin or rank-distribution plots are generated.

## Rank calculations

For OE, the observed vector is `-logFC_shrunk`. Predictor scores are oriented so higher
means stronger targeting. Non-finite scores are excluded before duplicate predictions are
aggregated, keeping the strongest finite score per miRNA and gene.

The same usable-response gene mask is used for every candidate. Unscored genes form one
bottom tie. **Both** observed responses and predictions are converted to average ranks,
centered and normalized to unit length. Their dot product is Spearman correlation.
The fast engine is checked against the original per-candidate implementation.

Within-expression alignment is the gene-count-weighted mean of defined bin correlations,
with at least 20 usable genes per bin. Equal baseline values stay together. Unknown baseline
values have no bin, but usable responses still enter the overall metric. Within-bin decoys
must use the same contributing bins as the cognate.

Specificity uses different nucleotide 2–8 seeds. Ties count half. Decoy z-scores are unavailable
with fewer than two usable decoys or no decoy variation. Counts refer to miRNAs, not independent
seed families. These comparisons are descriptive, not calibrated p-values. Shared-control
comparisons are kept together in the metadata; no independent-experiment significance is claimed.

## Fast repeats and report-only runs

The first pass streams just three prediction columns and filters each chunk into a compact
miRNA-by-gene matrix. ZIP/gzip inputs are read directly without extracting enormous predictor TSVs.
Rank calculations are batched across miRNAs and reused for repeated gene masks, with a bounded
rank-cache memory budget. Predictors are processed sequentially to keep memory use predictable.

Repeat the full command to reuse downloads, matrices, analysis tables and figures. Input file
size/mtime, genes, miRNAs, directions, filters and algorithm version determine cache identity.
Changes to inputs invalidate the relevant cache. Do not edit cached files in place while preserving
their original size and timestamps. On POSIX systems concurrent explorers serialize cache writes.

For report regeneration from an existing run, without network or prediction parsing:

```bash
uv run funmirbench-vector-explore --replot results/vector_explore/<timestamp>/response_vectors
```

The saved numerical tables are sufficient for replotting. For faster draft plots use `--dpi 90`;
the default is 130 and `--dpi 250` is available for higher-resolution exports. Plot DPI/detail
choices do not change the numerical analyses. `--refresh-metadata` refreshes the pinned records'
file registries; it does not select a newer record implicitly.

## Checks that are still pending

Replicate stability, controls-only contrasts and a seed-site baseline need extra sample-level
or UTR/site inputs. They are described as pending in each report and never presented as
computed results. Published tables are a snapshot, not proof that the final RNA-seq processing
has been reproduced. Agreement with transcriptomic responses does not establish direct binding.
