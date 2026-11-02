# Detailed report methods

The detailed reports extend the v2 benchmark. They do not restore the retired DE/FDR labelling rules. The primary evaluation remains the full pair set on the fixed gene list and `state_d05` from the processing pipeline. There are no new model downloads or calls to the processing pipeline.

## Report structure and styling

The main PDF ends with an index of experiment, predictor and supplementary reports. Each detailed PDF uses the same A4 page size, margins, blue headings, body font, table rendering and plot placement helpers as the main PDF. Plot sizes, axis fonts and predictor colours are imported from the main plotting module, including single-predictor plots. Heatmaps use a value colour scale rather than predictor colours.

Experiment reports begin with the experiment ID, miRNA, cell line, QC action, dependency group and gene count in Run settings. They then explain the evaluation, show the plots, list their metric details, and link the predictor-experiment reports. Predictor reports inventory their per-experiment results. Missing curves remain visible as unavailable, with reasons and state counts in the metrics tables.

## Pair sets and labels

| Pair set | Genes ranked | Missing predictions |
|---|---|---|
| FPS (primary) | All genes in the fixed benchmark universe | One tied block below every scored gene |
| IPS (sensitivity) | Genes scored by every selected tool for the experiment's miRNA | None; an empty intersection is reported explicitly |
| Own scored (coverage diagnostic) | Genes scored by the individual tool for the experiment's miRNA | Omitted |

IPS is an intersection of **pairs in each experiment**, not an intersection of genes that happened to receive scores for different miRNAs in other experiments. Selecting a different tool set can change IPS. Own-scored curves can have different prevalences and are not a common-universe leaderboard.

For PR, ROC, enrichment and false-positive burden, positives are `responsive` and negatives are `no_response_above_margin`. Unresolved genes are excluded. The configured `evaluation.min_per_state` applies separately in every pair set; at least one gene per state is always required. No-response includes opposite-direction responses and is not proof of an unchanged transcript.

Continuous diagnostics use finite oriented shrunken log2FC values with `eligible = true`. They include eligible unresolved genes. A larger oriented effect means a stronger change in the repression direction. These measurements describe transcriptomic response, not direct binding.

## Restored plots

| Plot | Definition and interpretation |
|---|---|
| PR curves | Whole score blocks enter at one threshold. Dashed baselines are the responsive fraction of the corresponding labelled pair set. |
| ROC curves | Same labelled rows and thresholds as PR; the diagonal is random ordering. |
| Average precision and PR-AUC | AP is the sum of recall increments times precision at the attained threshold. PR-AUC uses trapezoidal integration, with the conventional (recall 0, precision 1) endpoint. They are exported separately; interpolation can inflate trapezoidal PR-AUC for coarse scores, so AP is primary. |
| Score versus response | Within-experiment average score percentile against oriented shrunken log2FC. Responsive points are coloured; other states are grey. The tables also retain raw oriented scores. |
| GSEA-style enrichment | Unweighted running hit fraction minus miss fraction, evaluated at whole score-block boundaries. This is a descriptive signed ES, without permutation p-values or normalized ES. A constant ranking has ES zero. |
| Top-target effect CDF | Expected inclusion among the top `reports.top_n` genes, averaging over orders at tied cutoffs. The black dashed reference is all eligible genes. Selection uses the full ranking before unavailable effects are removed; selected/usable weights and medians are in `top_effect_summary.tsv`. A right shift indicates stronger repression-consistent effects. The plotted x-range is restricted for readability; exported CDFs retain the tails. |
| Recovery | Fraction of responsive genes recovered at a prediction budget on the complete ranking. Unresolved and ineligible genes consume budget. A cutoff inside a tied block counts its expected share. |
| Precision and false-positive burden | All predictors are displayed. Precision counts labelled genes only; score ties enter together, so achieved recall may exceed 50%. False positives per true positive are `(selected - TP) / TP`. |
| Gene heatmaps | Up to 40 genes selected by maximum percentile across tools, and separately responsive genes in any tool's top decile. Gene IDs only resolve display selection ties. These capped images are not used to compute metrics. |
| Predictor correlations | Spearman agreement on all eligible genes, and separately on each pair's shared top decile. Cells with fewer than three genes or constant ranks are NA. Agreement does not demonstrate biological accuracy. |
| Cross-experiment distributions | AP, PR-AUC, AUROC, signed Spearman, Pearson and ES, separately for each universe. Points are experiments; all/unflagged group-bootstrap summaries are in `universe_summary.tsv`. |
| Combination frontier | Exploratory maximum of within-tool scored percentiles for each subset, with missing scores zero, evaluated against fixed FPS labels. Single tools retain their predictor colours; combinations are grey diamonds. Lines show nondominated mean coverage/AP candidates. Exhaustive combinations are limited to eight selected tools. These are descriptive, untrained ensembles without a claim of generalization or significance. |

Mean metric uncertainty is computed by resampling dependency groups, using the same bootstrap method as the main benchmark.

## Legacy and unmerged PR audit

Sources inspected: the pre-v2 `evaluate_plots.py`, `evaluate_reports.py`, `benchmark_reports.py`, and PRs [#113](https://github.com/BioGeMT/fun-miRBench/pull/113), [#116](https://github.com/BioGeMT/fun-miRBench/pull/116), [#117](https://github.com/BioGeMT/fun-miRBench/pull/117), [#118](https://github.com/BioGeMT/fun-miRBench/pull/118).

- **#113:** retains the distinction between own-scored, IPS and FPS, separate AP/PR-AUC summaries, metric distributions and coverage/performance diagnostics. The new IPS uses the actual per-experiment pair intersection. We retain signed Spearman and do not label squared rank correlation as explained variance.
- **#116:** retains a dashed all-gene CDF reference, compact effect-range presentation and exported top-effect summaries. Top-N selection now averages score ties rather than using gene ID order.
- **#117:** restores per-experiment recovery and an all-predictor precision/false-positive burden comparison, rather than showing only the winner. Separate experiment PDFs replace rasterized multi-experiment manuscript mosaics.
- **#118:** retains FPS recovery, precision at 50% recall and false positives per true positive. These agree with the existing main recovery outputs; the new detailed curves add the complete threshold path.
- **Legacy reports:** restores individual predictor-experiment classification/scatter/enrichment diagnostics, common/own-scored comparisons, top-target CDFs, algorithms-versus-genes and top-positive heatmaps, predictor rank correlations, and exploratory combination frontiers.

Local/global rank distribution CDFs and count histograms are intentionally omitted. Score percentiles are retained only for scatter plots, heatmaps and the explicitly defined ensemble diagnostic.

No old pipeline or figure script is merged wholesale. Old OE/KO harmonization, raw-DE/FDR labels, unresolved-as-background treatment, gene-order tie breaking, and squared-rank-correlation-as-explained-variance interpretations are not reintroduced.

## Reproducible outputs

`reports/index.tsv` lists all generated PDFs. `universe_metrics.tsv` contains each experiment/tool/pair-set record, including counts and unavailable-curve reasons. `universe_summary.tsv` contains all/unflagged means, dependency-group bootstrap intervals and contributing experiment/group counts. `top_effect_summary.tsv` and `combination_metrics.tsv` contain the continuous and ensemble diagnostics. Individual plot directories hold threshold paths, CDF coordinates, gene matrices and predictor-correlation matrices.

The main metrics and original output tables are unchanged. Configure `reports.detailed: false` to omit the detailed reports entirely.
