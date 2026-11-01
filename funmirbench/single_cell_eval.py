"""Evaluation of the selected single-cell comparisons: tables, two plots and REPORT.pdf.

Main and supplement comparisons go through the bulk evaluator unchanged, one row per named miRNA. Machinery comparisons
have no single miRNA: predicted targets of the families a cell line expresses are scored against the same number of
families it barely expresses. Seed effects and the control check come from pipelines/single_cell/evaluation, because
they need the 3'UTR sequences and the cell stores; they are reported and never used to keep or drop a comparison.
"""
from __future__ import annotations

import argparse
import logging
import pathlib
from concurrent.futures import ProcessPoolExecutor

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import rankdata

from funmirbench.benchmark import (expression_tenths, filter_df, load_predictions, rank_last,
                                   resolve_predictor_output_path, run_benchmark, score_metrics, spearman)
from funmirbench.logger import parse_log_level, setup_logging
from funmirbench.plots import (BAR_ALPHA, LARGE, SMALL, draw_boxes, draw_dots, drawn_tools, save, style_box_axes,
                               tool_colors)
from funmirbench.report_pdf import (BLOCK_GAP, LEFT as PAGE_LEFT, ROW_HEIGHT, TOP as PAGE_TOP, WIDTH as PAGE_WIDTH,
                                    _bullet_block, _draw_basic_table, _header, _new_page, _save_page, draw_inventory,
                                    draw_plot, draw_tables)

logger = logging.getLogger(__name__)

MAKING = {"DROSHA", "DGCR8", "DICER1", "XPO5"}
LINE_OF_SCREEN = {"Replogle2022_K562_GW": "K562", "Replogle2022_K562_essential": "K562", "VIPerturb_K562": "K562",
                  "XAtlasOrion_HEK293T": "HEK293T", "GSE218566_HEK293": "HEK293T", "Nadig2025_HepG2": "HepG2",
                  "GSE345058_A549": "A549", "XAtlasOrion_HCT116": "HCT116", "GSE249595_Jurkat_activated": "Jurkat",
                  "GSE249595_Jurkat_untreated": "Jurkat", "Nadig2025_Jurkat": "Jurkat"}
EXPRESSED_SHARE, DECOY_SHARE, N_DECOYS = 0.005, 0.0001, 20
METRICS = (("auroc_d05", "AUROC", 0.5), ("spearman_shrunk", "Spearman rho", 0.0), ("decoy_z_shrunk", "Decoy z-score", 0.0))
PANEL_W, PANEL_H, PANEL_LEFT, PANEL_GAP, PANEL_TOP, PANEL_RIGHT = 4.6, 3.85, 1.0, 0.95, 1.05, 0.2  # inches
MEASURE_LINES = ["AUROC: chance that a responsive gene scores above a no-response gene; 0.5 is chance.",
                 "Spearman rho: agreement of the score's ranking with the shrunken log2FC; 0 is none.",
                 "Decoy z-score: Spearman for the tested miRNA against different-seed miRNAs, in SDs."]
OUTPUT_TABLES = ["per_comparison.tsv: every metric, for each main and supplement experiment and tool.",
                 "machinery_per_comparison.tsv: every metric, for each machinery experiment and tool.",
                 "machinery_per_screen.tsv: medians per screen, the points of the machinery plot.",
                 "machinery_summary.tsv: medians per tool, by gene type.",
                 "data_checks.tsv: cells, replicates, states, seed effect and target change per experiment."]


def read_comparisons(final_dir: pathlib.Path, selected_path: pathlib.Path) -> pd.DataFrame:
    """The selected comparisons that the final table holds, with their group; each must be fitted."""
    selected = pd.read_csv(selected_path, sep="\t", dtype=str, keep_default_na=False)[["comparison_id", "group"]]
    final = pd.read_csv(final_dir / "comparisons.tsv", sep="\t", dtype=str, keep_default_na=False).drop(columns="group", errors="ignore")
    rows = final.merge(selected, on="comparison_id")
    if not (rows.group == "main").any():
        raise ValueError(f"{final_dir}: no main comparison in the final table")
    unfitted = rows.comparison_id[rows.status != "fitted"].tolist()
    if unfitted:
        raise ValueError(f"selected comparisons not fitted: {unfitted}")
    return rows


def evaluator_input(rows: pd.DataFrame, final_dir: pathlib.Path, mirnas: pd.DataFrame, families: pd.DataFrame,
                    tools: pd.DataFrame, evaluation: dict, out: pathlib.Path) -> pd.DataFrame:
    """Input folder and config of the bulk evaluator, one row per named miRNA. All are written as OE because the
    evaluator scores only OE rows and its measures use oriented values; orientation keeps the real direction."""
    out.mkdir(parents=True, exist_ok=True)
    comparisons, experiments = [], []
    for row in rows[rows.group != "machinery"].itertuples(index=False):
        named = row.mirna_or_gene.split(";")
        for mirna in named:
            experiment = f"{row.comparison_id}__{mirna}" if len(named) > 1 else row.comparison_id
            comparisons.append({"experiment_id": experiment, "comparison_id": row.comparison_id, "study": row.dataset, "mirna_name": mirna,
                                "perturbation": "OE", "orientation": row.orientation, "cell_line": row.dataset, "split": "single_cell",
                                "design": row.design, "dependency_group": row.dataset, "action": "keep" if row.group == "main" else "supplement",
                                "flags": row.replicate_type})
            experiments.append({"experiment_id": experiment, "table": str(final_dir / row.final_dir / row.table)})
    comparisons = pd.DataFrame(comparisons)
    comparisons.to_csv(out / "comparisons.tsv", sep="\t", index=False)
    pd.DataFrame(experiments).to_csv(out / "experiments.tsv", sep="\t", index=False)
    missing = sorted(set(comparisons.mirna_name) - set(mirnas.mirna_name))
    sequences = families.set_index("mirna_name").mirna_sequence.reindex(missing)
    added = pd.DataFrame({"mirna_name": missing, "mirna_sequence": sequences.to_numpy(), "seed_2_8": sequences.str[1:8].to_numpy(),
                          "seed_family": "single_cell"})
    pd.concat([mirnas, added]).to_csv(out / "mirnas.tsv", sep="\t", index=False)
    tools.to_csv(out / "predictions_info.tsv", sep="\t", index=False)
    config = {"input_dir": str(out), "mirnas_tsv": str(out / "mirnas.tsv"), "predictions_tsv": str(out / "predictions_info.tsv"),
              "predictors": {"tool_id": tools.tool_id.tolist()}, "evaluation": evaluation, "out_dir": str(out / "results")}
    (out / "benchmark.yaml").write_text(yaml.safe_dump(config, sort_keys=False))
    return comparisons


def family_sets(profiles: pd.DataFrame) -> dict[str, dict]:
    """Per cell line: expressed families with weights, and decoy sets of as many barely expressed families."""
    rng = np.random.default_rng(0)
    sets = {}
    for line, group in sorted(profiles.groupby("line")):
        share = group.set_index("seed").share.astype(float)
        expressed = share[share >= EXPRESSED_SHARE].sort_values(ascending=False)
        low = sorted(share[share < DECOY_SHARE].index)
        if len(low) < len(expressed):
            raise ValueError(f"{line}: {len(low)} low families for {len(expressed)} expressed")
        sets[line] = {"weights": (expressed / expressed.sum()).to_numpy(), "expressed": expressed.index.tolist(),
                      "decoys": [list(rng.choice(low, size=len(expressed), replace=False)) for _ in range(N_DECOYS)]}
    return sets


def combine_families(seeds: list[str], weights: np.ndarray, members: dict[str, list[str]], ranks: dict[str, np.ndarray],
                     n_genes: int) -> np.ndarray:
    """Weighted mean over families of each gene's rank, a family being the mean of its scored miRNAs."""
    total, used = np.zeros(n_genes), 0.0
    for seed, weight in zip(seeds, weights):
        scored = [ranks[name] for name in members.get(seed, []) if name in ranks]
        if scored:
            total += weight * np.mean(scored, axis=0)
            used += weight
    return total / used


def score_tool(args) -> list[dict]:
    tool, path, direction, members, sets, genes, rows, min_per_state = args
    wanted = {name for s in sets.values() for seeds in [s["expressed"], *s["decoys"]] for seed in seeds for name in members.get(seed, [])}
    by_mirna = load_predictions(path, direction, wanted, set(genes))
    ranks = {name: rankdata(rank_last(score, genes).to_numpy(float)) / len(genes) for name, score in by_mirna.items()}
    scores = {line: {"cognate": combine_families(s["expressed"], s["weights"], members, ranks, len(genes)),
                     "decoys": [combine_families(seeds, s["weights"], members, ranks, len(genes)) for seeds in s["decoys"]]}
              for line, s in sets.items()}
    records = []
    for row in rows:
        table = pd.read_csv(row["path"], sep="\t", usecols=["gene_id", "eligible", "oriented_logFC_shrunk", "control_mean_normalized_count",
                                                            "state_d05"]).set_index("gene_id")
        table["eligible"] = table.eligible.astype(bool)
        table["expression_tenth"] = expression_tenths(table.control_mean_normalized_count.fillna(0.0))
        table = table.reindex(genes)
        metrics = score_metrics(pd.Series(scores[row["line"]]["cognate"], index=genes), table, min_per_state=min_per_state)
        est = table.eligible.to_numpy(bool)
        y = table.oriented_logFC_shrunk.to_numpy(float)[est]
        decoys = np.array([spearman(d[est], y) for d in scores[row["line"]]["decoys"]])
        records.append({key: row[key] for key in ("comparison_id", "screen", "line", "gene", "gene_type")} |
                       {"predictor": tool, **metrics, "decoy_z_shrunk": (metrics["spearman_shrunk"] - decoys.mean()) / decoys.std(ddof=1)})
    logger.info("scored %d machinery comparisons with %s", len(records), tool)
    return records


def score_machinery(rows: pd.DataFrame, final_dir: pathlib.Path, profiles: pd.DataFrame, families: pd.DataFrame, tools: pd.DataFrame,
                    root: pathlib.Path, min_per_state: int) -> pd.DataFrame:
    members = families.groupby("seed").mirna_name.apply(list).to_dict()
    sets = family_sets(profiles)
    machinery = [{"comparison_id": row.comparison_id, "screen": row.dataset, "line": LINE_OF_SCREEN[row.dataset], "gene": row.mirna_or_gene,
                  "gene_type": "miRNA-making" if row.mirna_or_gene in MAKING else "repression", "path": str(final_dir / row.final_dir / row.table)}
                 for row in rows[rows.group == "machinery"].itertuples(index=False)]
    genes = pd.Index(pd.read_csv(machinery[0]["path"], sep="\t", usecols=["gene_id"]).gene_id)
    jobs = [(tool.tool_id, resolve_predictor_output_path(tool.predictor_output_path, root), tool.score_direction, members, sets, genes, machinery,
             min_per_state) for tool in tools.itertuples(index=False)]
    with ProcessPoolExecutor(len(jobs)) as pool:
        return pd.DataFrame([record for records in pool.map(score_tool, jobs) for record in records])


def plot_panels(tools: pd.DataFrame, per: pd.DataFrame, out: pathlib.Path, stem: str, *, boxes: pd.DataFrame | None = None,
                note: bool = False) -> None:
    """AUROC, Spearman and decoy z-score in one row. Boxes: the boxes data, or per where it has at least 10 points."""
    labels = drawn_tools(tools)
    tool_ids = list(labels)
    colors = tool_colors(tool_ids)
    bottom = 1.75 if note else 1.05
    width, height = PANEL_LEFT + 3 * PANEL_W + 2 * PANEL_GAP + PANEL_RIGHT, PANEL_TOP + PANEL_H + bottom
    fig = plt.figure(figsize=(width, height))
    for i, (metric, title, reference) in enumerate(METRICS):
        ax = fig.add_axes([(PANEL_LEFT + i * (PANEL_W + PANEL_GAP)) / width, bottom / height, PANEL_W / width, PANEL_H / height])
        behind = boxes if boxes is not None else per
        if behind[metric].notna().groupby(behind.predictor).sum().min() >= 10:
            draw_boxes(ax, behind, tool_ids, colors, metric)
        draw_dots(ax, per, tool_ids, colors, metric)
        style_box_axes(ax, len(tool_ids))
        ax.axhline(reference, color="#9CA3AF", linewidth=1.0, zorder=0)
        ax.set_title(title, pad=14, fontsize=LARGE, fontweight="bold")
    span = (PANEL_LEFT / width, 0.0, (3 * PANEL_W + 2 * PANEL_GAP) / width)
    handles = [plt.Rectangle((0, 0), 1, 1, color=colors[tool], alpha=BAR_ALPHA) for tool in tool_ids]
    fig.legend(handles, labels.values(), loc="upper center", frameon=False, borderaxespad=0.0, borderpad=0.0, fontsize=SMALL, handlelength=1.6,
               mode="expand", ncol=len(handles), bbox_to_anchor=(*span, (bottom - 0.25) / height))
    if note:
        notes = [plt.Rectangle((0, 0), 1, 1, color="#555555", alpha=0.28, label=f"Box: {boxes.experiment_id.nunique()} bulk RNA-seq experiments"),
                 plt.Line2D([], [], marker="o", linestyle="", color="#555555", markersize=8, label=f"Dot: {per.experiment_id.nunique()} single-cell experiments")]
        fig.legend(handles=notes, loc="upper center", frameon=False, borderaxespad=0.0, borderpad=0.0, fontsize=SMALL - 2, handletextpad=0.4,
                   ncol=2, columnspacing=3.0, bbox_to_anchor=(0.5, (bottom - 0.8) / height))
    save(fig, out, stem)


def data_checks(rows: pd.DataFrame, seed: pd.DataFrame) -> pd.DataFrame:
    columns = ["comparison_id", "group", "dataset", "orientation", "n_treated_cells", "n_control_cells", "n_treated_units", "n_control_units",
               "replicate_type", "design", "d05_responsive", "d05_no_response_above_margin", "d05_unresolved", "target_log2fc", "target_log2fc_shrunk", "caveats"]
    return rows[columns].merge(seed[["comparison_id", "mirna", "effect_7mer_8mer", "ci_low", "ci_high", "rank", "families", "passes"]], on="comparison_id", how="left")


def number(value, digits: int = 3) -> str:
    return "NA" if pd.isna(value) else f"{float(value):.{digits}f}"


def measure_rows(per: pd.DataFrame, names: pd.Series, bulk_means: pd.DataFrame, tools: pd.DataFrame, metric: str) -> list[list[str]]:
    wide = per[per.experiment_id.isin(names.index)].pivot(index="experiment_id", columns="predictor", values=metric).reindex(names.index)
    rows = [[names[experiment], *(number(wide.loc[experiment, tool]) for tool in tools.tool_id)] for experiment in names.index]
    return rows + [["Bulk mean", *(number(bulk_means.loc[tool, metric]) for tool in tools.tool_id)]]


def write_report(path: pathlib.Path, *, rows: pd.DataFrame, per: pd.DataFrame, comparisons: pd.DataFrame, bulk: pd.DataFrame, machinery: pd.DataFrame,
                 screens: pd.DataFrame, summary: pd.DataFrame, checks: pd.DataFrame, seed_bulk: pd.DataFrame, control: pd.DataFrame, labels: pd.DataFrame,
                 tools: pd.DataFrame, n_genes: int, min_per_state: int, plots: pathlib.Path) -> None:
    tool_labels = drawn_tools(tools)
    main = rows[rows.group == "main"]
    main_seed = checks[checks.group == "main"]
    names = labels[labels.comparison_id.isin(main.comparison_id)].set_index("comparison_id").label
    order = {comparison: position for position, comparison in enumerate(labels.comparison_id)}
    bulk_means = bulk.groupby("predictor")[[metric for metric, _, _ in METRICS]].mean()
    making = machinery[machinery.gene_type == "miRNA-making"]
    making_screens = screens[screens.gene_type == "miRNA-making"]
    medians = making_screens.groupby("predictor").decoy_z_shrunk.median()
    real = (main.replicate_type == "biological").sum()
    key = [f"Main experiments with real replicates: {real} of {len(main)}; the others are pseudo-replicates of one sample per group.",
           f"Seed check: {int(main_seed.passes.fillna(False).astype(bool).sum())} of {len(main)} main experiments pass, and {int(seed_bulk.passes.sum())} of {len(seed_bulk)} bulk experiments.",
           f"Control check: up to {int(control.responsive_overexpression_direction.max()):,} genes responsive between two libraries with nothing perturbed.",
           "Median decoy z-score, miRNA-making screens: " + ", ".join(f"{tool_labels[t]} {number(medians.get(t), 1)}" for t in tool_labels) + ".",
           "Too few main experiments to rank predictors, so no experiments-won count and no mean over them."]
    pages_tools = [tool_labels[t] for t in tools.tool_id]
    tool_widths = [max(0.12, (len(name) + 1) / 70) for name in pages_tools]  # draw_tables wraps at 70 characters per page width
    widths = [1 - sum(tool_widths), *tool_widths]
    wide = [(title, ["Experiment", *pages_tools], widths, measure_rows(per, names, bulk_means, tools, metric)) for metric, title, _ in METRICS]
    summary_rows = [[row.gene_type, tool_labels[row.predictor], str(row.comparisons), number(row.median_spearman), number(row.median_decoy_z, 2),
                     f"{number(row.median_auroc)} ({row.auroc_comparisons})"] for row in summary.itertuples(index=False)]
    seed_rows = [[labels.set_index("comparison_id").label.get(row.comparison_id, row.comparison_id), row.mirna.removeprefix("hsa-"), number(row.effect_7mer_8mer), f"{int(row['rank']):,}", "yes" if row.passes else "no"]
                 for _, row in checks[checks.group != "machinery"].sort_values("comparison_id", key=lambda ids: ids.map(order), kind="stable").iterrows()]
    target_rows = []
    for gene, group in checks[checks.comparison_id.isin(making.comparison_id)].assign(gene=lambda d: d.comparison_id.map(machinery.drop_duplicates("comparison_id").set_index("comparison_id").gene)).groupby("gene"):
        target_rows.append([gene, str(len(group)), number(pd.to_numeric(group.target_log2fc).median(), 2), number(pd.to_numeric(group.target_log2fc_shrunk).median(), 2)])
    control_rows = [[row.control_pair, f"{row.cells_side_1:,} and {row.cells_side_2:,}", f"{row.responsive_knockout_direction:,}", f"{row.responsive_overexpression_direction:,}"] for row in control.itertuples(index=False)]
    labelled = labels.set_index("comparison_id")
    inventory = []
    for row in rows.sort_values("group", key=lambda groups: groups.map({"main": 0, "supplement": 1, "machinery": 2}), kind="stable").itertuples(index=False):
        if row.group == "machinery":
            inventory.append([f"{row.mirna_or_gene} · {LINE_OF_SCREEN[row.dataset]} ({'plotted' if row.mirna_or_gene in MAKING else 'tables only'})", row.comparison_id])
        else:
            inventory.append([f"{labelled.label[row.comparison_id]} against {labelled.control[row.comparison_id]} ({row.group}). {row.caveats[0].upper()}{row.caveats[1:]}.", row.comparison_id])
    with PdfPages(path) as pdf:
        fig, ax = _new_page()
        y = _header(ax, "Single-cell Evaluation Report")
        settings = [["Main experiments", str(len(main))], ["Machinery experiments", f"{(rows.group == 'machinery').sum()} in {machinery.screen.nunique()} screens"],
                    ["Supplement experiments", str((rows.group == "supplement").sum())], ["Predictors evaluated", str(len(tool_labels))],
                    ["Responsive genes", "oriented shrunken log2FC > 0; FSOS s-value < 0.05 against d = 0.5"],
                    ["No-response genes", "whole shrunken 95% interval of the oriented change below d = 0.5"]]
        height = ROW_HEIGHT * (len(settings) + 1)
        _draw_basic_table(ax, settings, columns=["Run setting", "Value"], col_widths=[0.27, 0.73], bbox=[PAGE_LEFT, y - height, PAGE_WIDTH, height])
        y = _bullet_block(ax, "Evaluation design",
                          [f"Cells are summed per replicate unit and fitted with the bulk DESeq2 and labelling scripts, on the same {n_genes:,} genes.",
                           "Main experiments go through the bulk evaluator unchanged; a gene a tool does not score ranks last.",
                           f"AUROC needs at least {min_per_state} responsive and {min_per_state} no-response genes.",
                           "AUROC of miR-874 and miR-365 mimic is unreliable: treated and control cells are separate libraries.",
                           "Machinery: expressed families (at least 0.5% of reads) against 20 sets of families below 0.01%.",
                           "Seed effects and the control check are reported and never keep or drop an experiment."], y=y - height - BLOCK_GAP)
        _save_page(pdf, fig)

        fig, ax = _new_page()
        y = _bullet_block(ax, "Key results", key, y=PAGE_TOP)
        _bullet_block(ax, "Output tables", OUTPUT_TABLES, y=y)
        _save_page(pdf, fig)

        with plt.rc_context({"savefig.dpi": 300}):
            fig, ax = _new_page()
            y = draw_plot(fig, ax, plots, "single_cell", "Main experiments", ["Each box is the 21 bulk experiments, each dot one main experiment.", *MEASURE_LINES], PAGE_TOP)
            draw_plot(fig, ax, plots, "machinery", "Machinery experiments",
                      ["Each dot is one screen, the median of its miRNA-making experiments, because they share control cells.",
                       "Repression genes are left off the plot, a choice made after seeing their weak signal; they are in the tables."], y)
            _save_page(pdf, fig)

        draw_tables(pdf, [*wide, ("Machinery, medians over experiments", ["Gene type", "Predictor", "Experiments", "Spearman", "Decoy z-score", "AUROC (n)"],
                                  [0.2, 0.2, 0.17, 0.13, 0.14, 0.16], summary_rows)])
        draw_tables(pdf, [("Seed check", ["Experiment", "miRNA", "Effect (log2FC)", "Rank of 2,022", "Passes"], [0.42, 0.2, 0.14, 0.14, 0.1], seed_rows),
                          ("Target change, miRNA-making genes", ["Gene", "Experiments", "Median log2FC", "Median shrunken log2FC"], [0.25, 0.2, 0.25, 0.3], target_rows),
                          ("Control check: responsive genes with nothing perturbed", ["Two sides", "Cells per side", "Knockout direction", "Overexpression direction"], [0.42, 0.18, 0.18, 0.22], control_rows)])
        draw_inventory(pdf, inventory)


def run_evaluation(config_path: pathlib.Path) -> pathlib.Path:
    config_path = config_path.expanduser().resolve()
    config = yaml.safe_load(config_path.read_text())
    root = config_path.parent
    resolve = lambda key: resolve_predictor_output_path(config[key], root)
    final_dir, tables_dir, out = resolve("final_dir"), resolve("tables_dir"), resolve("out_dir")
    evaluation = config.get("evaluation") or {}
    min_per_state = int(evaluation.get("min_per_state", 30))
    out.mkdir(parents=True, exist_ok=True)

    rows = read_comparisons(final_dir, resolve("selected"))
    tools = filter_df(pd.read_csv(resolve("predictions_tsv"), sep="\t"), config["predictors"])
    tools["predictor_output_path"] = [str(resolve_predictor_output_path(path, root)) for path in tools.predictor_output_path]
    mirnas = pd.read_csv(resolve("mirnas_tsv"), sep="\t")
    families = pd.read_csv(tables_dir / "mirna_families.tsv", sep="\t")
    comparisons = evaluator_input(rows, final_dir, mirnas, families, tools, evaluation, out / "evaluator")
    results = run_benchmark(out / "evaluator" / "benchmark.yaml")
    per = pd.read_csv(results / "per_comparison.tsv", sep="\t").merge(comparisons[["experiment_id", "comparison_id", "orientation"]], on="experiment_id")
    per.to_csv(out / "per_comparison.tsv", sep="\t", index=False)
    bulk = pd.read_csv(resolve("bulk_per_comparison"), sep="\t")

    profiles = pd.read_csv(tables_dir / "mirna_profiles.tsv", sep="\t")
    machinery = score_machinery(rows, final_dir, profiles, families, tools, root, min_per_state)
    machinery.to_csv(out / "machinery_per_comparison.tsv", sep="\t", index=False)
    metric_columns = [metric for metric, _, _ in METRICS]
    screens = machinery.groupby(["screen", "line", "gene_type", "predictor"], as_index=False)[metric_columns].median()
    screens.to_csv(out / "machinery_per_screen.tsv", sep="\t", index=False)
    summary = machinery.groupby(["gene_type", "predictor"]).agg(comparisons=("comparison_id", "nunique"), screens=("screen", "nunique"), median_spearman=("spearman_shrunk", "median"),
                                                                 median_decoy_z=("decoy_z_shrunk", "median"), median_auroc=("auroc_d05", "median"),
                                                                 auroc_comparisons=("auroc_d05", "count")).reset_index()
    summary.to_csv(out / "machinery_summary.tsv", sep="\t", index=False)

    seed = pd.read_csv(tables_dir / "seed_check.tsv", sep="\t")
    checks = data_checks(rows, seed)
    checks.to_csv(out / "data_checks.tsv", sep="\t", index=False)
    plots = out / "plots"
    main = per[per.experiment_id.isin(comparisons.experiment_id[comparisons.action == "keep"])]
    plot_panels(tools, main, plots, "single_cell", boxes=bulk, note=True)
    plot_panels(tools, screens[screens.gene_type == "miRNA-making"], plots, "machinery")
    write_report(out / "REPORT.pdf", rows=rows, per=per, comparisons=comparisons, bulk=bulk, machinery=machinery, screens=screens, summary=summary, checks=checks,
                 seed_bulk=pd.read_csv(tables_dir / "seed_check_bulk.tsv", sep="\t"), control=pd.read_csv(tables_dir / "control_check.tsv", sep="\t"),
                 labels=pd.read_csv(tables_dir / "experiment_labels.tsv", sep="\t"), tools=tools, n_genes=int(per.n_genes.iloc[0]), min_per_state=min_per_state, plots=plots)
    logger.info("wrote %s", out)
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description="Evaluate the predictors on the selected single-cell comparisons.")
    parser.add_argument("--config", type=pathlib.Path, required=True)
    parser.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"])
    args = parser.parse_args()
    setup_logging(parse_log_level(args.log_level))
    run_evaluation(args.config)


if __name__ == "__main__":
    main()
