"""Experiment, predictor and supplementary reports on the v2 benchmark inputs."""
from __future__ import annotations

import hashlib
import itertools
import logging
import pathlib
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

from funmirbench import subreport_plots as plots
from funmirbench.plots import draw_box, drawn_tools, tool_colors, two_panels
from funmirbench.report_pdf import (
    BOTTOM, INVENTORY_WRAP_SCALE, LEFT, ROW_HEIGHT, TOP, WIDTH,
    _bullet_block, _draw_basic_table, _header, _new_page, _paginate_table_rows,
    _save_page, _table_row_weight, draw_plot,
)
from funmirbench.subreport_metrics import (
    UNIVERSES, UNIVERSE_LABELS, coverage_frontier, metrics, threshold_curves, top_weights, weighted_cdf,
)

logger = logging.getLogger(__name__)


def component(name: str) -> str:
    """Stable filenames without allowing an input ID to escape the report folder."""
    safe = re.sub(r"[^A-Za-z0-9_.-]", "_", str(name)).strip(".")
    if not safe or safe != str(name) or len(safe) > 120:
        safe = safe[:100] + "_" + hashlib.sha256(str(name).encode()).hexdigest()[:10]
    return safe


def export(frame: pd.DataFrame, path: pathlib.Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False)


def number(value) -> str:
    return "NA" if pd.isna(value) else f"{value:.3f}"


def pdf_report(path, title, settings, specs, frame, *, links=()):
    """Use the main PDF's pages, bullets, tables, margins and image placement."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(path) as pdf:
        fig, ax = _new_page()
        y = _header(ax, title)
        y = _bullet_block(ax, "Run settings", settings, y=y)
        y = _bullet_block(ax, "Evaluation design", [
            "Positive: responsive at d = 0.5. Negative: no response above margin; unresolved genes are excluded.",
            "FPS ranks unscored genes last; IPS uses pairs scored by every selected predictor in this experiment.",
            "Own-scored curves use different gene sets and baselines; they are coverage diagnostics.",
            "Continuous diagnostics use eligible shrunken estimates. Enrichment is descriptive, without an NES or p-value.",
        ], y=y)
        _bullet_block(ax, "Output tables", ["metrics.tsv and the plot-data TSVs beside each figure contain the plotted numbers.",
                     "Average precision (AP) and trapezoidal PR-AUC are distinct; AP is the primary PR summary."], y=y)
        _save_page(pdf, fig)
        for start in range(0, len(specs), 2):
            fig, ax = _new_page()
            y = TOP
            for image, label, lines in specs[start:start + 2]:
                y = draw_plot(fig, ax, image.parent, image.stem, label, lines, y)
                if y < BOTTOM:
                    raise ValueError(f"Report page overflow: {path}: {label}")
            _save_page(pdf, fig)
        columns = ["Experiment", "Predictor", "Pair set", "n (+ / -)", "AP", "AUROC", "Spearman"]
        rows = [[r.experiment_id, r.predictor, r.universe, f"{r.n_responsive} / {r.n_no_response}", number(r.average_precision),
                 number(r.auroc), number(r.spearman)] for r in frame.itertuples(index=False)]
        widths = [.29, .17, .1, .14, .1, .1, .1]
        for page in _paginate_table_rows(rows, col_widths=widths, max_weight=18):
            fig, ax = _new_page()
            y = _header(ax, "Metric details")
            weight = 1 + sum(_table_row_weight(row, col_widths=widths, wrap_scale=96) for row in page)
            h = min(y - BOTTOM, ROW_HEIGHT * weight)
            _draw_basic_table(ax, page, columns=columns, col_widths=widths, bbox=[LEFT, y - h, WIDTH, h])
            _save_page(pdf, fig)
        if links:
            report_index_pages(pdf, links, base=path.parent)


def report_index_pages(pdf, entries, *, base):
    """Clickable PDF index, also used in the main report."""
    import os
    rows = [[e["name"], os.path.relpath(e["path"], base)] for e in entries]
    targets = dict((row[1], row[1]) for row in rows)
    widths = [.55, .45]
    for page in _paginate_table_rows(rows, col_widths=widths, wrap_scale=INVENTORY_WRAP_SCALE, max_weight=18):
        fig, ax = _new_page()
        y = _header(ax, "Detailed report index")
        weight = 1 + sum(_table_row_weight(row, col_widths=widths, wrap_scale=INVENTORY_WRAP_SCALE) for row in page)
        h = min(y - BOTTOM, ROW_HEIGHT * weight)
        table = _draw_basic_table(ax, page, columns=["Report", "Relative path"], col_widths=widths,
                                 bbox=[LEFT, y - h, WIDTH, h], wrap_scale=INVENTORY_WRAP_SCALE)
        for i, row in enumerate(page, 1):
            table[i, 1].get_text().set_url(targets[row[1]])
        _save_page(pdf, fig)


def write_subreports(tables, comps, preds, tools, out, *, min_per_state=30, top_n=100, bootstrap=2000, seed=0):
    # Imported at call time, after the benchmark module has loaded.
    from funmirbench.benchmark import group_bootstrap, rank_last, recall_at_budgets
    labels, colors = drawn_tools(tools), tool_colors(list(preds))
    report_root = out / "reports"
    index, records, cdf_records = [], [], []
    pair_records, combination_records = [], []
    for comparison in comps.itertuples(index=False):
        exp = comparison.experiment_id
        logger.info("Detailed reports: %s", exp)
        table = tables[exp]
        folder = report_root / "experiments" / component(exp)
        raw = pd.DataFrame({tool: by_mirna.get(comparison.mirna_name, pd.Series(dtype=float)).reindex(table.index)
                            for tool, by_mirna in preds.items()})
        ranked = pd.DataFrame({tool: rank_last(raw[tool], table.index) for tool in labels})
        labelled = table.state_d05.isin(["responsive", "no_response_above_margin"])
        pos = table.loc[labelled, "state_d05"].eq("responsive").to_numpy()
        if min(pos.sum(), (~pos).sum()) >= max(1, min_per_state) and len(labels) <= 8:
            # An exploratory max-percentile ensemble, not a trained/validated model.
            own_percentiles = raw.rank(method="average", pct=True).fillna(0)
            for count in range(1, len(labels) + 1):
                for members in itertools.combinations(labels, count):
                    curve = threshold_curves(own_percentiles.loc[labelled, list(members)].max(axis=1).to_numpy(), pos)
                    combination_records.append({"experiment_id": exp, "members": "+".join(members), "n_members": count,
                                                "action": comparison.action, "coverage": raw[list(members)].notna().any(axis=1).mean(),
                                                "average_precision": float(np.sum(np.diff(np.r_[0, curve.recall]) * curve.precision))})
        specs, exp_records, pair_links = [], [], []
        for universe in UNIVERSES:
            curves, recs = {}, {}
            for tool in labels:
                rec, curve = metrics(table, raw, ranked[tool], tool, universe, min_per_state)
                rec.update(experiment_id=exp, action=comparison.action, dependency_group=comparison.dependency_group)
                recs[tool], curves[tool] = rec, curve
                records.append(rec); exp_records.append(rec)
                export(curve, folder / "plots" / f"{universe}_{component(tool)}_curves.tsv")
            image = plots.classification(curves, recs, labels, folder / "plots", f"{universe}_classification", universe)
            specs.append((image, UNIVERSE_LABELS[universe], ["PR baselines are dashed; tied genes enter at one threshold.",
                          "Missing curves have too few labelled genes; counts and reasons are in metrics.tsv."]))
            image = plots.enrichment(curves, labels, folder / "plots", f"{universe}_enrichment_burden")
            specs.append((image, "Enrichment and false-positive burden", [f"Pair set: {UNIVERSE_LABELS[universe]}.",
                          "Running ES uses score blocks, so gene order inside ties cannot create enrichment."]))
            if universe == "fps":
                for tool in labels:
                    pair = folder / "predictors" / component(tool)
                    pair_frame = pd.DataFrame([recs[tool]])
                    single = {tool: labels[tool]}
                    first = plots.classification({tool: curves[tool]}, {tool: recs[tool]}, single, pair / "plots",
                                                 "classification", universe, palette=colors)
                    second = plots.scatter_and_enrichment(table, raw[tool], ranked[tool], curves[tool], tool,
                                                           labels[tool], colors[tool], pair / "plots")
                    export(pair_frame, pair / "metrics.tsv")
                    export(curves[tool], pair / "plots" / "classification.tsv")
                    export(pd.DataFrame({"gene_id": table.index, "score": raw[tool], "ranked_score": ranked[tool],
                                         "effect": table.oriented_logFC_shrunk, "eligible": table.eligible,
                                         "state": table.state_d05}), pair / "plots" / "score_and_enrichment.tsv")
                    pdf_report(pair / "REPORT.pdf", "Predictor-experiment report",
                               [f"Experiment: {exp}", f"miRNA: {comparison.mirna_name}; predictor: {labels[tool]}",
                                f"QC action: {comparison.action}; pair set: FPS"],
                               [(first, "Classification", ["Unresolved genes are excluded from PR and ROC."]),
                                (second, "Score and enrichment", ["Coloured points: responsive; grey: no-response or unresolved."])], pair_frame)
                    pair_links.append({"kind": "predictor_experiment", "name": labels[tool], "path": pair / "REPORT.pdf"})
        cdfs = {"background": weighted_cdf(table.loc[table.eligible, "oriented_logFC_shrunk"].to_numpy(float))}
        recovery = []
        for tool in labels:
            weights = top_weights(ranked[tool], top_n)
            valid_weights = weights.where(table.eligible, 0).to_numpy(float)
            cdfs[tool] = weighted_cdf(table.oriented_logFC_shrunk.to_numpy(float), valid_weights)
            cdf_records.append({"experiment_id": exp, "predictor": tool, "budget": top_n,
                                "selected_weight": weights.sum(), "usable_weight": cdfs[tool].weight.sum(),
                                "median_effect": float(cdfs[tool].loc[cdfs[tool].cdf >= .5, "effect"].iloc[0]) if len(cdfs[tool]) else np.nan})
            export(cdfs[tool], folder / "plots" / f"top_effect_cdf_{component(tool)}.tsv")
            if min((table.state_d05 == "responsive").sum(), (table.state_d05 == "no_response_above_margin").sum()) >= max(1, min_per_state):
                recovery += [{"predictor": tool, "budget": b, "recall": r}
                             for b, r in recall_at_budgets(ranked[tool].to_numpy(), (table.state_d05 == "responsive").to_numpy()).items()]
        recovery = pd.DataFrame(recovery, columns=["predictor", "budget", "recall"])
        export(recovery, folder / "plots" / "recovery.tsv")
        export(cdfs["background"], folder / "plots" / "effect_cdf_background.tsv")
        image = plots.effect_and_recovery(cdfs, recovery, labels, folder / "plots")
        specs.append((image, "Effect distributions and practical recovery", [f"Top {top_n}: fractional inclusion at tied cutoffs; dashed CDF: all eligible genes.",
                      "CDF view uses the central effect range; TSVs retain the tails. Budgets include unresolved genes."]))
        local = ranked.rank(method="average", pct=True)
        usable = table.eligible.to_numpy(bool) & np.isfinite(table.oriented_logFC_shrunk)
        gene_heatmaps = []
        for positive_only in [False, True]:
            candidate = local.loc[table.state_d05 == "responsive"] if positive_only else local.loc[usable]
            if positive_only:
                candidate = candidate.loc[candidate.ge(.9).any(axis=1)]
            order = pd.DataFrame({"strength": candidate.max(axis=1), "gene": candidate.index.astype(str)}, index=candidate.index)
            ids = order.sort_values(["strength", "gene"], ascending=[False, True]).head(40).index
            matrix = local.loc[ids].T
            gene_heatmaps.append(matrix)
            export(matrix.rename_axis("predictor").reset_index(), folder / "plots" / f"{'top_positive' if positive_only else 'algorithms_vs_genes'}.tsv")
        image = plots.heatmaps(gene_heatmaps, ["Algorithms versus genes", "Top-decile responsive genes"], labels, folder / "plots", "gene_heatmaps")
        specs.append((image, "Gene-level ranks", ["At most 40 genes displayed, selected by the strongest percentile across tools.",
                      "Gene IDs resolve display-only ties; numerical diagnostics use score blocks or average ranks."]))
        correlations = [local.loc[usable].corr(method="spearman"), pd.DataFrame(np.nan, index=list(labels), columns=list(labels))]
        for a, b in itertools.product(labels, repeat=2):
            common = usable & local[a].ge(.9) & local[b].ge(.9)
            if common.sum() >= 3 and local.loc[common, a].nunique() > 1 and local.loc[common, b].nunique() > 1:
                correlations[1].loc[a, b] = local.loc[common, a].corr(local.loc[common, b], method="spearman")
        for name, matrix in zip(["all", "shared_top_decile"], correlations):
            export(matrix.rename_axis("predictor").reset_index(), folder / "plots" / f"predictor_correlations_{name}.tsv")
        image = plots.heatmaps(correlations, ["All eligible genes", "Shared top-decile genes"], labels, folder / "plots", "predictor_correlations")
        specs.append((image, "Predictor rank agreement", ["Spearman correlations describe agreement, not target validity.",
                      "Shared-top cells are NA with fewer than three genes or constant ranks."]))
        exp_frame = pd.DataFrame(exp_records)
        export(exp_frame, folder / "metrics.tsv")
        pdf_report(folder / "REPORT.pdf", "Experiment report", [f"Experiment: {exp}",
                   f"miRNA: {comparison.mirna_name}; cell line: {comparison.cell_line}",
                   f"QC action: {comparison.action}; dependency group: {comparison.dependency_group}",
                   f"Genes: {len(table):,}; top-effect budget: {top_n}"], specs, exp_frame, links=pair_links)
        index.append({"kind": "experiment", "name": exp, "path": folder / "REPORT.pdf"})
        pair_records.extend(pair_links)
    frame = pd.DataFrame(records)
    export(frame, report_root / "universe_metrics.tsv")
    export(pd.DataFrame(cdf_records), report_root / "top_effect_summary.tsv")
    summaries = []
    rng = np.random.default_rng(seed)
    for subset, sub in [("all", frame), ("unflagged", frame[frame.action == "keep"])]:
        for (tool, universe), group in sub.groupby(["predictor", "universe"]):
            for metric in ["average_precision", "pr_auc", "auroc", "spearman", "pearson", "enrichment_score"]:
                v = group[metric]
                lo, hi = group_bootstrap(v, group.dependency_group, n=bootstrap, rng=rng)
                summaries.append({"subset": subset, "predictor": tool, "universe": universe, "metric": metric,
                                  "mean": v.mean(), "ci95_low": lo, "ci95_high": hi,
                                  "n_comparisons": int(v.notna().sum()), "n_groups": group.loc[v.notna(), "dependency_group"].nunique()})
    export(pd.DataFrame(summaries), report_root / "universe_summary.tsv")
    export(frame, report_root / "diagnostics" / "metrics.tsv")
    overview_specs = []
    for universe in UNIVERSES:
        for pair in [("average_precision", "pr_auc"), ("auroc", "spearman"), ("pearson", "enrichment_score")]:
            fig, axes = two_panels()
            for ax, metric in zip(axes, pair):
                draw_box(ax, frame[frame.universe == universe], list(labels), colors, metric)
                ax.set_title(metric.replace("_", " ").title())
            image = plots.finish(fig, axes, labels, report_root / "diagnostics" / "plots", f"{universe}_{pair[0]}_{pair[1]}")
            overview_specs.append((image, UNIVERSE_LABELS[universe], ["Each dot is an experiment; use universe_summary.tsv for group-bootstrap intervals."]))
    precision = pd.read_csv(out / "precision_at_recall.tsv", sep="\t")
    fig, axes = two_panels()
    for ax, metric in zip(axes, ["precision", "false_positives_per_true_positive"]):
        draw_box(ax, precision, list(labels), colors, metric)
        ax.set_title("Precision at 50% recall" if metric == "precision" else "False positives / true positive")
    image = plots.finish(fig, axes, labels, report_root / "diagnostics" / "plots", "all_predictor_precision_burden")
    overview_specs.append((image, "All-predictor precision and burden", ["Every predictor is shown, rather than only the best predictor in each experiment.",
                          "Tied score blocks enter together; precision excludes unresolved genes."]))
    export(precision, report_root / "diagnostics" / "plots" / "all_predictor_precision_burden.tsv")
    combinations = pd.DataFrame(combination_records)
    if len(combinations):
        export(combinations, report_root / "combination_metrics.tsv")
        fig, axes = two_panels()
        for ax, subset in zip(axes, ["all", "unflagged"]):
            rows = combinations if subset == "all" else combinations[combinations.action == "keep"]
            mean = rows.groupby(["members", "n_members"])[["coverage", "average_precision"]].mean().reset_index()
            for r in mean.itertuples(index=False):
                ax.scatter(r.coverage, r.average_precision, s=75, marker="o" if r.n_members == 1 else "D",
                           color=colors[r.members] if r.n_members == 1 else "#6B7280", alpha=.7)
            if len(mean):
                frontier = coverage_frontier(mean)
                ax.plot(frontier.coverage, frontier.average_precision, "--", color="#1F2937", linewidth=1)
            ax.set(xlabel="Mean scored-pair coverage", ylabel="Mean average precision", title=subset.title())
        image = plots.finish(fig, axes, labels, report_root / "diagnostics" / "plots", "combination_frontier")
        overview_specs.append((image, "Exploratory predictor combinations", ["Diamonds: maximum within-tool scored percentile across a subset; circles: single tools.",
                              "FPS labels stay fixed. No fitted ensemble or statistical claim of superiority; at most eight tools."]))
    pdf_report(report_root / "diagnostics" / "REPORT.pdf", "Supplementary diagnostics", [f"Comparisons: {len(comps)}; predictors: {len(labels)}",
               "FPS is primary; IPS and own-scored sets are sensitivity analyses."], overview_specs, frame)
    index.append({"kind": "diagnostics", "name": "Supplementary diagnostics", "path": report_root / "diagnostics" / "REPORT.pdf"})
    for tool, label in labels.items():
        folder = report_root / "predictors" / component(tool)
        tool_frame = frame[frame.predictor == tool]
        export(tool_frame, folder / "metrics.tsv")
        links = [{"name": c.experiment_id, "path": report_root / "experiments" / component(c.experiment_id) / "predictors" / component(tool) / "REPORT.pdf"}
                 for c in comps.itertuples(index=False)]
        pdf_report(folder / "REPORT.pdf", "Predictor report", [f"Predictor: {label}", f"Comparisons: {len(comps)}",
                   "Each experiment has its own detailed classification, score and enrichment report."], [], tool_frame, links=links)
        index.append({"kind": "predictor", "name": label, "path": folder / "REPORT.pdf"})
    manifest = [{"kind": e["kind"], "name": e["name"], "report": str(e["path"].relative_to(out))} for e in index + pair_records]
    export(pd.DataFrame(manifest), report_root / "index.tsv")
    (report_root / "README.md").write_text("# Detailed reports\n\n" + "\n".join(
        f"- [{e['name']}]({e['path'].relative_to(report_root).as_posix()})" for e in index) + "\n", encoding="utf-8")
    return index
