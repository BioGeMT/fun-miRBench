"""Per-experiment and per-experiment/predictor reports on v2 benchmark inputs."""
from __future__ import annotations

import hashlib
import itertools
import logging
import pathlib
import re

import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

from funmirbench import subreport_plots as plots
from funmirbench.plots import drawn_tools, tool_colors
from funmirbench.report_pdf import (
    BLOCK_GAP, BOTTOM, INVENTORY_WRAP_SCALE, LEFT, REPORT_TABLE_SIZE, ROW_HEIGHT, TOP, WIDTH,
    _bullet_block, _draw_basic_table, _header, _new_page, _paginate_table_rows,
    _save_page, _table_row_weight, draw_plot,
)
from funmirbench.subreport_metrics import (
    UNIVERSES, UNIVERSE_LABELS, metrics, threshold_curves, top_weights, weighted_cdf,
)

logger = logging.getLogger(__name__)
RESPONSIVE_HEATMAP_LIMIT = 20


def responsive_gene_rank_data(table: pd.DataFrame, percentiles: pd.DataFrame, limit: int = RESPONSIVE_HEATMAP_LIMIT):
    """Return predictor percentile ranks for the strongest eligible responsive genes."""
    if limit < 1:
        raise ValueError("limit must be positive")
    effects = table["oriented_logFC_shrunk"].to_numpy(float)
    keep = (table["eligible"].to_numpy(bool) & np.isfinite(effects)
            & table["state_d05"].eq("responsive").to_numpy())
    candidates = pd.DataFrame({"gene_id": table.index[keep].astype(str), "effect": effects[keep]})
    selected = candidates.sort_values(["effect", "gene_id"], ascending=[False, True], kind="mergesort").head(limit)
    gene_ids = selected["gene_id"].tolist()
    effect_by_gene = pd.Series(selected["effect"].to_numpy(), index=gene_ids)
    return percentiles.loc[gene_ids], effect_by_gene


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


def key_results(frame: pd.DataFrame, predictor_labels: dict[str, str] | None = None) -> list[str]:
    """Summarize primary FPS results using the main report's comparison style."""
    predictor_labels = predictor_labels or {}
    fps = frame[frame.universe == "fps"]
    if fps.empty:
        return ["No FPS metric rows are available."]

    metrics = (("average_precision", "AP"), ("auroc", "AUROC"), ("spearman", "Spearman"))
    predictors, experiments = fps.predictor.nunique(), fps.experiment_id.nunique()
    if predictors == 1 and experiments == 1:
        row = fps.iloc[0]
        lines = [f"FPS {label}: {number(row[column])}." for column, label in metrics]
        if pd.isna(row.average_precision):
            lines.append(f"Classification unavailable: {row.classification_note}.")
        if pd.isna(row.spearman):
            lines.append(f"Continuous metrics unavailable: {row.spearman_na_reason}.")
        return lines

    if predictors == 1:
        tool = fps.predictor.iloc[0]
        label = predictor_labels.get(tool, tool)
        lines = []
        for column, name in metrics:
            values = fps[column].dropna()
            if len(values):
                lines.append(f"Mean FPS {name}: {number(values.mean())} over {len(values)} comparisons ({label}).")
            else:
                reason = fps.classification_note.iloc[0] if column in {"average_precision", "auroc"} else fps.spearman_na_reason.iloc[0]
                lines.append(f"No FPS {name} is available: {reason}.")
        return lines

    means = fps.groupby("predictor")[[column for column, _ in metrics]].mean()
    lines = []
    for column, name in metrics:
        available = means[column].dropna()
        if len(available):
            winner = available.idxmax()
            lines.append(f"Highest mean FPS {name}: {predictor_labels.get(winner, winner)} ({number(available[winner])}).")
        else:
            reason_column = "classification_note" if column in {"average_precision", "auroc"} else "spearman_na_reason"
            reasons = [str(value) for value in fps[reason_column].dropna().unique() if str(value)]
            lines.append(f"No FPS {name} is available: {reasons[0] if reasons else 'no usable comparisons'}.")
    return lines


def pdf_report(path, title, settings, specs, frame, *, predictor_labels=None):
    """Use the main PDF's pages, bullets, tables, margins and image placement."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(path) as pdf:
        fig, ax = _new_page()
        y = _header(ax, title)
        setting_rows = [setting.split(": ", 1) if ": " in setting else ["Description", setting]
                        for setting in settings]
        setting_height = ROW_HEIGHT * (len(setting_rows) + 1)
        _draw_basic_table(ax, setting_rows, columns=["Run setting", "Value"], col_widths=[.27, .73],
                          bbox=[LEFT, y - setting_height, WIDTH, setting_height], font_size=REPORT_TABLE_SIZE)
        y -= setting_height + BLOCK_GAP
        y = _bullet_block(ax, "Evaluation design", [
            "Positives: responsive. Negatives: no_response_above_margin. Unresolved genes are excluded from PR and ROC.",
            "FPS uses the fixed gene universe and ranks missing predictions last; IPS and own-scored are supplementary.",
            "Continuous metrics use eligible, finite oriented shrunken log2FC. Enrichment is descriptive, without an NES or p-value.",
        ], y=y)
        y = _bullet_block(ax, "Key results", key_results(frame, predictor_labels), y=y)
        _bullet_block(ax, "Output tables", ["metrics.tsv and plot-data TSVs contain the plotted values and unavailable-metric reasons.",
                     "Average precision (AP) and trapezoidal PR-AUC are distinct; AP is the primary PR summary."], y=y)
        _save_page(pdf, fig)
        fig = None
        for image, label, lines in specs:
            if fig is None:
                fig, ax = _new_page()
                y, count = TOP, 0
            axes_before, texts_before = len(fig.axes), len(ax.texts)
            next_y = draw_plot(fig, ax, image.parent, image.stem, label, lines, y)
            if next_y < BOTTOM and count:
                # Undo this block before saving the preceding plots on their page.
                for added in list(fig.axes)[axes_before:]:
                    added.remove()
                for added in list(ax.texts)[texts_before:]:
                    added.remove()
                _save_page(pdf, fig)
                fig, ax = _new_page()
                next_y = draw_plot(fig, ax, image.parent, image.stem, label, lines, TOP)
                count = 0
            if next_y < BOTTOM:
                raise ValueError(f"Report page overflow: {path}: {label}")
            y, count = next_y, count + 1
            if count == 2:
                _save_page(pdf, fig)
                fig = None
        if fig is not None:
            _save_page(pdf, fig)
        columns = ["Experiment", "Predictor", "Pair set", "n (+ / -)", "AP", "AUROC", "Spearman"]
        pair_set = {"fps": "FPS", "ips": "IPS", "own_scored": "Own"}
        rows = [[str(r.experiment_id).replace("_", "_ ") if len(str(r.experiment_id)) > 22 else r.experiment_id,
                 (predictor_labels or {}).get(r.predictor, r.predictor), pair_set[r.universe],
                 f"{r.n_responsive} / {r.n_no_response}", number(r.average_precision),
                 number(r.auroc), number(r.spearman)] for r in frame.itertuples(index=False)]
        widths = [.30, .14, .09, .15, .09, .09, .14]
        for page in _paginate_table_rows(rows, col_widths=widths, wrap_scale=74, max_weight=18):
            fig, ax = _new_page()
            y = _header(ax, "Metric details")
            weight = 1 + sum(_table_row_weight(row, col_widths=widths, wrap_scale=74) for row in page)
            h = min(y - BOTTOM, ROW_HEIGHT * weight)
            _draw_basic_table(ax, page, columns=columns, col_widths=widths, bbox=[LEFT, y - h, WIDTH, h],
                              wrap_scale=74)
            _save_page(pdf, fig)
def _write_experiment_subreports(task):
    """Write one experiment's complete report tree; workers never share output paths."""
    comparison, table, raw, labels, colors, report_root, min_per_state, top_n = task
    from funmirbench.benchmark import rank_last, recall_at_budgets

    exp = comparison["experiment_id"]
    logger.info("Detailed reports: %s", exp)
    folder = report_root / "experiments" / component(exp)
    overview_plots = folder / "plots" / "overview"
    predictor_plots = folder / "plots" / "predictors"
    plot_data = folder / "plots" / "data"
    ranked = pd.DataFrame({tool: rank_last(raw[tool], table.index) for tool in labels})
    specs, exp_records = [], []
    for universe in UNIVERSES:
        curves, recs = {}, {}
        for tool in labels:
            rec, curve = metrics(table, raw, ranked[tool], tool, universe, min_per_state)
            rec.update(experiment_id=exp, action=comparison["action"],
                       dependency_group=comparison["dependency_group"])
            recs[tool], curves[tool] = rec, curve
            exp_records.append(rec)
            export(curve, plot_data / f"{universe}_{component(tool)}_curves.tsv")
        image = plots.classification(curves, recs, labels, overview_plots, f"{universe}_classification")
        specs.append((image, UNIVERSE_LABELS[universe], ["PR baselines are dashed; tied genes enter at one threshold.",
                      "Missing curves have too few labelled genes; counts and reasons are in metrics.tsv."]))
        image = plots.enrichment(curves, labels, overview_plots, f"{universe}_enrichment_burden")
        specs.append((image, "Enrichment and false-positive burden", [f"Pair set: {UNIVERSE_LABELS[universe]}.",
                      "Running ES uses score blocks, so gene order inside ties cannot create enrichment."]))
        if universe == "fps":
            for tool in labels:
                pair = folder / "predictors" / component(tool)
                tool_plots = predictor_plots / component(tool)
                pair_frame = pd.DataFrame([recs[tool]])
                single = {tool: labels[tool]}
                first = plots.classification({tool: curves[tool]}, {tool: recs[tool]}, single, tool_plots,
                                             "classification", palette=colors)
                second = plots.scatter_and_enrichment(table, ranked[tool], curves[tool], tool,
                                                       labels[tool], colors[tool], tool_plots)
                export(pd.DataFrame({"gene_id": table.index, "score": raw[tool], "ranked_score": ranked[tool],
                                     "effect": table.oriented_logFC_shrunk, "eligible": table.eligible,
                                     "state": table.state_d05}), plot_data / f"{component(tool)}_score_and_enrichment.tsv")
                pdf_report(pair / "REPORT.pdf", f"{labels[tool]} predictor-experiment report",
                           [f"Experiment: {exp}", f"miRNA: {comparison['mirna_name']}; predictor: {labels[tool]}",
                            f"QC action: {comparison['action']}; pair set: FPS",
                            f"Minimum genes per classification state: {max(1, min_per_state)}"],
                           [(first, "Classification", ["Unresolved genes are excluded from PR and ROC."]),
                            (second, "Score and enrichment", ["Coloured points: responsive; grey: no-response or unresolved."])],
                           pair_frame, predictor_labels=labels)
    cdfs = {"background": weighted_cdf(table.loc[table.eligible, "oriented_logFC_shrunk"].to_numpy(float))}
    recovery = []
    for tool in labels:
        weights = top_weights(ranked[tool], top_n)
        valid_weights = weights.where(table.eligible, 0).to_numpy(float)
        cdfs[tool] = weighted_cdf(table.oriented_logFC_shrunk.to_numpy(float), valid_weights)
        export(cdfs[tool], plot_data / f"top_effect_cdf_{component(tool)}.tsv")
        if min((table.state_d05 == "responsive").sum(), (table.state_d05 == "no_response_above_margin").sum()) >= max(1, min_per_state):
            recovery += [{"predictor": tool, "budget": budget, "recall": value}
                         for budget, value in recall_at_budgets(
                             ranked[tool].to_numpy(), (table.state_d05 == "responsive").to_numpy()).items()]
    recovery = pd.DataFrame(recovery, columns=["predictor", "budget", "recall"])
    export(recovery, plot_data / "recovery.tsv")
    export(cdfs["background"], plot_data / "effect_cdf_background.tsv")
    image = plots.effect_and_recovery(cdfs, recovery, labels, overview_plots)
    specs.append((image, "Effect distributions and practical recovery", [f"Top {top_n}: fractional inclusion at tied cutoffs; dashed CDF: all eligible genes.",
                  "CDF view uses the central effect range; TSVs retain the tails. Budgets include unresolved genes."]))
    local = ranked.rank(method="average", pct=True)
    usable = table.eligible.to_numpy(bool) & np.isfinite(table.oriented_logFC_shrunk)
    rank_matrix, selected_effects = responsive_gene_rank_data(table, local, RESPONSIVE_HEATMAP_LIMIT)
    rank_export = rank_matrix.copy()
    rank_export.insert(0, "oriented_shrunken_log2FC", selected_effects)
    rank_export.index.name = "gene_id"
    export(rank_export.reset_index(), plot_data / "responsive_gene_ranks.tsv")
    image = plots.responsive_gene_heatmap(rank_matrix, labels, overview_plots)
    specs.append((image, "Responsive genes: predictor rank percentiles", [
        f"Up to {RESPONSIVE_HEATMAP_LIMIT} eligible responsive genes, ordered by oriented shrunken log2FC (largest first).",
        "Rows are genes; columns are predictor percentiles across the full gene universe. This descriptive view is not a performance metric.",
    ]))
    correlations = [local.loc[usable].corr(method="spearman"), pd.DataFrame(np.nan, index=list(labels), columns=list(labels))]
    for a, b in itertools.product(labels, repeat=2):
        common = usable & local[a].ge(.9) & local[b].ge(.9)
        if common.sum() >= 3 and local.loc[common, a].nunique() > 1 and local.loc[common, b].nunique() > 1:
            correlations[1].loc[a, b] = local.loc[common, a].corr(local.loc[common, b], method="spearman")
    for name, matrix in zip(["all", "shared_top_decile"], correlations):
        export(matrix.rename_axis("predictor").reset_index(), plot_data / f"predictor_correlations_{name}.tsv")
    image = plots.heatmaps(correlations, ["All eligible genes", "Shared top-decile genes"], labels,
                           overview_plots, "predictor_correlations")
    specs.append((image, "Predictor rank agreement", ["Spearman correlations describe agreement, not target validity.",
                  "Shared-top cells are NA with fewer than three genes or constant ranks."]))
    exp_frame = pd.DataFrame(exp_records)
    export(exp_frame, folder / "metrics.tsv")
    pdf_report(folder / "REPORT.pdf", "Experiment report", [f"Experiment: {exp}",
               f"miRNA: {comparison['mirna_name']}; cell line: {comparison['cell_line']}",
               f"QC action: {comparison['action']}; dependency group: {comparison['dependency_group']}",
               f"Genes: {len(table):,}; top-effect budget: {top_n}",
               f"Minimum genes per classification state: {max(1, min_per_state)}"],
               specs, exp_frame, predictor_labels=labels)


def write_subreports(tables, comps, preds, tools, out, *, min_per_state=30, top_n=100, workers=2):
    # Imported at call time, after the benchmark module has loaded.
    from concurrent.futures import ProcessPoolExecutor

    labels, colors = drawn_tools(tools), tool_colors(list(preds))
    report_root = out / "reports"
    tasks = []
    for comparison in comps.itertuples(index=False):
        comp = comparison._asdict()
        table = tables[comp["experiment_id"]]
        raw = pd.DataFrame({tool: by_mirna.get(comp["mirna_name"], pd.Series(dtype=float)).reindex(table.index)
                            for tool, by_mirna in preds.items()})
        tasks.append((comp, table, raw, labels, colors, report_root, min_per_state, top_n))

    if workers == 1 or len(tasks) < 2:
        for task in tasks:
            _write_experiment_subreports(task)
    else:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            list(executor.map(_write_experiment_subreports, tasks, chunksize=1))

