"""Evaluation of the selected single-cell comparisons: tables, two figures and REPORT.pdf.

Main and supplement comparisons go through the bulk evaluator unchanged, one row per named miRNA. Machinery comparisons
have no single miRNA: predicted targets of the families a cell line expresses are scored against the same number of
families it barely expresses. Single cell is supporting evidence, so the figures draw no ranking and no mean over
experiments. Seed effects come from pipelines/single_cell/evaluation, because they need the 3'UTR sequences and the
cell stores. Seed checks decide which comparisons the figures draw: the named miRNA's seed effect for
main comparisons, and for machinery the site genes of the cell's profiled miRNA families against those of rare families.
Every comparison stays in the output tables.
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

from funmirbench.benchmark import (decoy_z, expression_tenths, filter_df, load_predictions, rank_last,
                                   resolve_predictor_output_path, run_benchmark, score_metrics)
from funmirbench.logger import parse_log_level, setup_logging
from funmirbench.plots import DARK_ALPHA, LARGE, PAIR_SIZE, SMALL, drawn_tools, save, tool_colors
from funmirbench.report_pdf import (BLOCK_GAP, LEFT as PAGE_LEFT, ROW_HEIGHT, TOP as PAGE_TOP, WIDTH as PAGE_WIDTH,
                                    _bullet_block, _draw_basic_table, _header, _new_page, _save_page, draw_inventory,
                                    draw_plot, draw_tables)

logger = logging.getLogger(__name__)

MAKING = {"DROSHA", "DGCR8", "DICER1", "XPO5"}
# screen: name in the report, starting with its cell line, and its miRNA profile; HEK293 has no profile of its own and uses HEK293T's
SCREENS = {"Replogle2022_K562_GW": ("K562 genome-wide", "K562"), "Replogle2022_K562_essential": ("K562 essential-genes", "K562"), "VIPerturb_K562": ("K562 VIPerturb", "K562"),
           "XAtlasOrion_HEK293T": ("HEK293T X-Atlas", "HEK293T"), "GSE218566_HEK293": ("HEK293", "HEK293T"), "Nadig2025_HepG2": ("HepG2", "HepG2"),
           "GSE345058_A549": ("A549", "A549"), "XAtlasOrion_HCT116": ("HCT116 X-Atlas", "HCT116"), "GSE249595_Jurkat_activated": ("Jurkat activated", "Jurkat"),
           "GSE249595_Jurkat_untreated": ("Jurkat untreated", "Jurkat"), "Nadig2025_Jurkat": ("Jurkat", "Jurkat")}
# low guide load, at least 3 replicate units per side, not flagged by the control check (K562 essential not tested itself, K562 GW of the same study was), CRISPRi: the screens the figure draws
FIGURE_SCREENS = ("Replogle2022_K562_GW", "Replogle2022_K562_essential", "GSE218566_HEK293", "Nadig2025_HepG2", "Nadig2025_Jurkat")
EXPRESSED_SHARE, DECOY_SHARE, N_DECOYS = 0.005, 0.0001, 20
# AUROC stays in the output tables but is not drawn: unreliable or below the gene minimum in most comparisons (see the report)
METRICS = (("spearman_shrunk", "Spearman rho", 0.0), ("decoy_z_shrunk", "Decoy z-score", 0.0))
Z95 = 1.959964
SPEARMAN_LINE = "Spearman rho: how well the score's ranking matches the shrunken log2FC. 0 is none."
MAIN_CAPTION = ["One group per experiment.", SPEARMAN_LINE,
                "Decoy z-score: Spearman for the tested miRNA against miRNAs with a different seed, in standard deviations."]
MACHINERY_CAPTION = ["One group per knocked-down gene, the median over the screens that knocked it down.", SPEARMAN_LINE,
                     f"Decoy z-score: Spearman for the profiled families against {N_DECOYS} sets of rare families, in standard deviations."]


def read_comparisons(final_dir: pathlib.Path, selected_path: pathlib.Path) -> pd.DataFrame:
    """The selected comparisons with their group; the final table must hold all of them, fitted."""
    selected = pd.read_csv(selected_path, sep="\t", dtype=str, keep_default_na=False)[["comparison_id", "group"]]
    final = pd.read_csv(final_dir / "comparisons.tsv", sep="\t", dtype=str, keep_default_na=False).drop(columns="group", errors="ignore")
    rows = final.merge(selected, on="comparison_id")
    missing = sorted(set(selected.comparison_id) - set(rows.comparison_id))
    if missing:
        raise ValueError(f"{final_dir}: {len(missing)} selected comparisons are not in the final table, for example {missing[:3]}; run run.sh for all groups")
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


def scored_weight(seeds: list[str], weights: np.ndarray, members: dict[str, list[str]], ranks: dict[str, np.ndarray]) -> float:
    return float(sum(weight for seed, weight in zip(seeds, weights) if any(name in ranks for name in members.get(seed, []))))


def score_tool(args) -> tuple[list[dict], list[dict]]:
    tool, path, direction, members, sets, genes, rows, min_per_state = args
    wanted = {name for s in sets.values() for seeds in [s["expressed"], *s["decoys"]] for seed in seeds for name in members.get(seed, [])}
    by_mirna = load_predictions(path, direction, wanted, set(genes))
    ranks = {name: rankdata(rank_last(score, genes).to_numpy(float)) / len(genes) for name, score in by_mirna.items()}
    scores = {line: {"cognate": combine_families(s["expressed"], s["weights"], members, ranks, len(genes)),
                     "decoys": [combine_families(seeds, s["weights"], members, ranks, len(genes)) for seeds in s["decoys"]]}
              for line, s in sets.items()}
    coverage = [{"profile": line, "predictor": tool, "expressed_weight_scored": scored_weight(s["expressed"], s["weights"], members, ranks),
                 "lowest_decoy_weight_scored": min(scored_weight(seeds, s["weights"], members, ranks) for seeds in s["decoys"])} for line, s in sets.items()]
    records = []
    for row in rows:
        table = pd.read_csv(row["path"], sep="\t", usecols=["gene_id", "eligible", "oriented_logFC_shrunk", "control_mean_normalized_count",
                                                            "state_d05"]).set_index("gene_id")
        table["eligible"] = table.eligible.astype(bool)
        table["expression_tenth"] = expression_tenths(table.control_mean_normalized_count.fillna(0.0))
        table = table.reindex(genes)
        cognate = scores[row["profile"]]["cognate"]
        metrics = score_metrics(pd.Series(cognate, index=genes), table, min_per_state=min_per_state)
        est = table.eligible.to_numpy(bool)
        y = table.oriented_logFC_shrunk.to_numpy(float)[est]
        records.append({key: row[key] for key in ("comparison_id", "screen", "profile", "gene", "gene_type")} |
                       {"predictor": tool, **metrics, "decoy_z_shrunk": decoy_z(cognate[est], [d[est] for d in scores[row["profile"]]["decoys"]], y)})
    logger.info("scored %d machinery comparisons with %s", len(records), tool)
    return records, coverage


def score_machinery(rows: pd.DataFrame, final_dir: pathlib.Path, profiles: pd.DataFrame, families: pd.DataFrame, tools: pd.DataFrame,
                    root: pathlib.Path, min_per_state: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    members = families.groupby("seed").mirna_name.apply(list).to_dict()
    sets = family_sets(profiles)
    machinery = [{"comparison_id": row.comparison_id, "screen": row.dataset, "profile": SCREENS[row.dataset][1], "gene": row.mirna_or_gene,
                  "gene_type": "miRNA-making" if row.mirna_or_gene in MAKING else "repression", "path": str(final_dir / row.final_dir / row.table)}
                 for row in rows[rows.group == "machinery"].itertuples(index=False)]
    genes = pd.Index(pd.read_csv(machinery[0]["path"], sep="\t", usecols=["gene_id"]).gene_id)
    jobs = [(tool.tool_id, resolve_predictor_output_path(tool.predictor_output_path, root), tool.score_direction, members, sets, genes, machinery,
             min_per_state) for tool in tools.itertuples(index=False)]
    with ProcessPoolExecutor(len(jobs)) as pool:
        done = list(pool.map(score_tool, jobs))
    return pd.DataFrame([record for records, _ in done for record in records]), pd.DataFrame([row for _, coverage in done for row in coverage])


def target_drop(rows: pd.DataFrame, final_dir: pathlib.Path) -> pd.DataFrame:
    """The targeted gene's unshrunken log2FC and upper 95% limit. Shrunken values are not used: shrinkage can erase a real drop."""
    records = []
    for row in rows[rows.group == "machinery"].itertuples(index=False):
        table = pd.read_csv(final_dir / row.final_dir / row.table, sep="\t", usecols=["gene_id", "logFC", "lfcSE"]).set_index("gene_id")
        logfc, se = table.logFC.get(row.target_used_gene_id, np.nan), table.lfcSE.get(row.target_used_gene_id, np.nan)
        records.append({"comparison_id": row.comparison_id, "target_logFC": logfc, "target_lfcSE": se, "target_ci95_high": logfc + Z95 * se})
    return pd.DataFrame(records)


def figure_main(rows: pd.DataFrame, seed: pd.DataFrame) -> pd.DataFrame:
    """Main comparisons the first figure draws: those where the miRNA acted, the 95% interval of the seed effect above 0."""
    main = rows[rows.group == "main"].merge(seed[["comparison_id", "ci_low"]], on="comparison_id", how="left")
    return main[main.ci_low > 0]


def figure_machinery(rows: pd.DataFrame, drop: pd.DataFrame, machinery_seed: pd.DataFrame) -> pd.DataFrame:
    """Comparisons the machinery figure draws: the figure's screens, a clear drop of the targeted transcript, and a passed
    machinery seed check (site genes of the profiled families rose more than those of the rare sets)."""
    machinery = rows[rows.group == "machinery"].merge(drop, on="comparison_id")
    passed = machinery_seed.comparison_id[machinery_seed.passes.astype(bool)]
    return machinery[machinery.dataset.isin(FIGURE_SCREENS) & (machinery.target_ci95_high < 0) & machinery.comparison_id.isin(passed)]


def axis_label(text: str) -> str:
    return text.replace(" ", "\n")


def single_cell_tables(per: pd.DataFrame, labels: pd.DataFrame, drawn: pd.DataFrame, tools: pd.DataFrame) -> dict[str, pd.DataFrame]:
    main = labels[labels.comparison_id.isin(drawn.comparison_id)].set_index("comparison_id").label
    tables = {}
    for metric, _, _ in METRICS:
        wide = per[per.experiment_id.isin(main.index)].pivot(index="experiment_id", columns="predictor", values=metric).reindex(main.index).rename(index=main)
        tables[metric] = wide[tools.tool_id.tolist()]
    return tables


def machinery_tables(machinery: pd.DataFrame, kept: pd.DataFrame, tools: pd.DataFrame) -> tuple[dict[str, pd.DataFrame], int]:
    """Per gene type and knocked-down gene, the median over the screens that knocked it down. Each screen has one experiment per gene,
    so no two experiments of a gene share control cells. Also the number of miRNA-making rows."""
    drawn = machinery[machinery.comparison_id.isin(kept.comparison_id)]
    if drawn.duplicated(["gene", "screen", "predictor"]).any():
        raise ValueError("a screen has more than one drawn experiment of one gene; its experiments share control cells")
    tables, n_making = {}, 0
    for metric, _, _ in METRICS:
        parts = [drawn[drawn.gene_type == gene_type].groupby(["gene", "predictor"])[metric].median().unstack() for gene_type in ("miRNA-making", "repression")]
        n_making = len(parts[0])
        tables[metric] = pd.concat(parts)[tools.tool_id.tolist()]
    return tables, n_making


def plot_bars(tables: dict[str, pd.DataFrame], tools: pd.DataFrame, out: pathlib.Path, stem: str, *, divider: float | None = None, groups: tuple[tuple[float, str], ...] = ()) -> None:
    """One bar per tool in each group of the x axis, starting at the value that means no signal; one metric per row, at the
    width, font sizes and bar shading of the bulk plots, so that the text comes out the same size in the report."""
    labels = drawn_tools(tools)
    colors = tool_colors(list(labels))
    fig, axes = plt.subplots(len(METRICS), 1, figsize=(PAIR_SIZE[0], 10.5))
    for ax, (metric, title, start) in zip(axes, METRICS):
        table = tables[metric]
        width = 0.84 / len(labels)
        for i, tool in enumerate(labels):
            x = np.arange(len(table)) + (i - (len(labels) - 1) / 2) * width
            ax.bar(x, table[tool] - start, bottom=start, width=width * 0.92, color=colors[tool], alpha=DARK_ALPHA, zorder=3)
        ax.axhline(start, color="#9CA3AF", linewidth=1.0, zorder=0)
        if divider is not None:
            ax.axvline(divider, color="#9CA3AF", linewidth=1.0)
        ax.set_xlim(-0.6, len(table) - 0.4)
        ax.set_xticks(np.arange(len(table)), [axis_label(name) for name in table.index], fontsize=SMALL)
        ax.tick_params(axis="y", labelsize=SMALL)
        ax.tick_params(axis="x", length=0)
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", alpha=0.25, linewidth=0.8)
        ax.set_axisbelow(True)
        ax.set_title(title, pad=14, fontsize=LARGE, fontweight="bold")
        for center, text in groups:
            ax.text(center, -0.12, text, transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=SMALL, fontweight="bold")
    fig.tight_layout(rect=(0, 0.06, 1, 1), h_pad=2.0)
    span = axes[-1].get_position()
    handles = [plt.Rectangle((0, 0), 1, 1, color=colors[tool], alpha=DARK_ALPHA) for tool in labels]
    fig.legend(handles, labels.values(), loc="lower center", mode="expand", ncol=len(handles), frameon=False, borderaxespad=0.0,
               fontsize=SMALL, handlelength=1.6, bbox_to_anchor=(span.x0, 0.0, span.width, 0.04))
    save(fig, out, stem)


def data_checks(rows: pd.DataFrame, seed: pd.DataFrame, machinery_seed: pd.DataFrame, drop: pd.DataFrame, drawn: pd.DataFrame, labels: pd.DataFrame) -> pd.DataFrame:
    columns = ["comparison_id", "group", "dataset", "orientation", "n_treated_cells", "n_control_cells", "n_treated_units", "n_control_units",
               "replicate_type", "design", "d05_responsive", "d05_no_response_above_margin", "d05_unresolved", "caveats"]
    checks = rows[columns].merge(drop, on="comparison_id", how="left").assign(in_figure=lambda d: d.comparison_id.isin(drawn.comparison_id))
    written = labels.set_index("comparison_id").caveats
    checks["caveats"] = [written.get(comparison, caveats) for comparison, caveats in zip(checks.comparison_id, checks.caveats)]
    machinery_columns = {"difference": "machinery_seed_effect", "ci_low": "machinery_seed_ci_low", "ci_high": "machinery_seed_ci_high", "passes": "machinery_seed_passes"}
    return checks.merge(seed[["comparison_id", "mirna", "effect_7mer_8mer", "ci_low", "ci_high", "rank", "families", "passes"]], on="comparison_id", how="left").merge(
        machinery_seed[["comparison_id", *machinery_columns]].rename(columns=machinery_columns), on="comparison_id", how="left")


def number(value, digits: int = 3) -> str:
    return "NA" if pd.isna(value) else f"{float(value):.{digits}f}"


def report_tables(tables: dict[str, pd.DataFrame], tools: pd.DataFrame, first_column: str, prefix: str, n_making: int | None = None) -> list[tuple]:
    """The numbers each figure draws, one table per panel. n_making: how many rows are miRNA-making genes, the rest repression."""
    names = [drawn_tools(tools)[tool] for tool in tools.tool_id]
    widths = [max(0.12, (len(name) + 1) / 70) for name in names]  # draw_tables wraps at 70 characters per page width
    groups = lambda table: list(table.index) if n_making is None else [f"{'miRNA-making' if i < n_making else 'Repression'}, {name}" for i, name in enumerate(table.index)]
    return [(f"{prefix}: {title}", [first_column, *names], [1 - sum(widths), *widths],
             [[group, *(number(value) for value in values)] for group, values in zip(groups(tables[metric]), tables[metric].to_numpy())]) for metric, title, _ in METRICS]


def write_report(path: pathlib.Path, *, drawn_main: pd.DataFrame, single_cell: dict[str, pd.DataFrame], machinery: dict[str, pd.DataFrame], n_making: int, kept: pd.DataFrame,
                 checks: pd.DataFrame, labels: pd.DataFrame, tools: pd.DataFrame, n_genes: int, min_per_state: int, plots: pathlib.Path) -> None:
    """The report describes only the experiments the figures draw; the others stay in the output tables and EXPERIMENTS.md."""
    tool_labels = drawn_tools(tools)
    labelled = labels.set_index("comparison_id")
    main_checks = checks[checks.group == "main"].set_index("comparison_id")
    making = kept[kept.mirna_or_gene.isin(MAKING)]
    few_responsive = int((drawn_main.d05_responsive.astype(int) < min_per_state).sum() + (kept.d05_responsive.astype(int) < min_per_state).sum())
    span = lambda values: f"{number(np.nanmin(values), 1)} to {number(np.nanmax(values), 1)}"
    z, gene_z = single_cell["decoy_z_shrunk"], machinery["decoy_z_shrunk"]
    above = int((z > 0).sum().sum())
    n_repression = len(gene_z) - n_making
    main_line = (f"Main: every predictor's decoy z-score is above 0 in all {len(drawn_main)} experiments." if above == z.size
                 else f"Main: {above} of {z.size} decoy z-scores are above 0.")
    settings = [["Main experiments", str(len(drawn_main))],
                ["Machinery experiments", f"{len(kept)} ({len(making)} miRNA-making, {len(kept) - len(making)} repression) in {kept.dataset.nunique()} screens"],
                ["Predictors evaluated", str(len(tool_labels))],
                ["Responsive genes", "oriented shrunken log2FC > 0, FSOS s-value < 0.05 against d = 0.5"],
                ["No-response genes", "whole shrunken 95% interval of the oriented change below d = 0.5"]]
    design = [f"Supporting evidence only: too few experiments to rank predictors. "
              f"{int((main_checks.replicate_type[drawn_main.comparison_id] == 'technical').sum())} of the {len(drawn_main)} main experiments use pseudo-replicates.",
              f"Counts are summed per replicate unit and fitted with the bulk DESeq2 and labelling scripts. All predictors rank the same {n_genes:,} genes, "
              "with unscored genes last.",
              "A main experiment is included only if the miRNA acted: genes with a 7mer or 8mer site changed in the expected direction more than matched "
              "genes without one (95% interval above 0).",
              f"Machinery: the cell's profiled miRNA families (at least {EXPRESSED_SHARE:.1%} of reads) are weighted by abundance and compared with {N_DECOYS} "
              f"equal-size sets of rare families (below {DECOY_SHARE:.2%}).",
              "A knockdown is included only if site genes of the profiled families rose more than those of the rare families (95% interval above 0). "
              "Screens, experiments and both inclusion rules were chosen after seeing results.",
              f"AUROC is only in the output tables: {few_responsive} of {len(drawn_main) + len(kept)} experiments have fewer than {min_per_state} responsive genes."]
    key = [main_line,
           f"Machinery: per-gene median decoy z-scores are {span(gene_z.iloc[:n_making])} for the {n_making} miRNA-making genes and "
           f"{span(gene_z.iloc[n_making:])} for the {n_repression} repression gene{'s' if n_repression != 1 else ''}."]
    outputs = ["per_comparison.tsv: every metric, AUROC included, per experiment and predictor.",
               "machinery_per_comparison.tsv: the same for machinery experiments.",
               "machinery_family_coverage.tsv: family weight each predictor scores, per profile.",
               "data_checks.tsv: replicate units, gene counts, transcript drop and seed checks."]
    inventory = [[f"{labelled.label[comparison]} against {labelled.control[comparison]}. {labelled.caveats[comparison]}", comparison]
                 for comparison in labels.comparison_id if comparison in set(drawn_main.comparison_id)]
    # machinery in the order of the figure: miRNA-making genes, then repression genes, each gene in screen order
    order = kept.assign(repression=~kept.mirna_or_gene.isin(MAKING), screen_order=kept.dataset.map(FIGURE_SCREENS.index))
    for row in order.sort_values(["repression", "mirna_or_gene", "screen_order"]).itertuples(index=False):
        name, profile = SCREENS[row.dataset]
        note = "" if name.split()[0] == profile else f", scored with the {profile} miRNA profile"
        inventory.append([f"{row.mirna_or_gene} knockdown in the {name} screen{note}.", row.comparison_id])
    with PdfPages(path) as pdf:
        fig, ax = _new_page()
        y = _header(ax, "Single-cell Evaluation Report")
        height = ROW_HEIGHT * (len(settings) + 1)
        _draw_basic_table(ax, settings, columns=["Run setting", "Value"], col_widths=[0.27, 0.73], bbox=[PAGE_LEFT, y - height, PAGE_WIDTH, height])
        y = _bullet_block(ax, "Evaluation design", design, y=y - height - BLOCK_GAP)
        _bullet_block(ax, "Key results", key, y=y)
        _save_page(pdf, fig)

        with plt.rc_context({"savefig.dpi": 300}):
            fig, ax = _new_page()
            draw_plot(fig, ax, plots, "single_cell", "Main experiments", MAIN_CAPTION, PAGE_TOP)
            _save_page(pdf, fig)
            fig, ax = _new_page()
            y = draw_plot(fig, ax, plots, "machinery", "Machinery experiments", MACHINERY_CAPTION, PAGE_TOP)
            _bullet_block(ax, "Output tables", outputs, y=y)  # page 1 has no room left for it
            _save_page(pdf, fig)

        draw_tables(pdf, report_tables(single_cell, tools, "Experiment", "Main"))
        draw_tables(pdf, report_tables(machinery, tools, "Gene", "Machinery", n_making))
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
    labels = pd.read_csv(tables_dir / "experiment_labels.tsv", sep="\t")
    comparisons = evaluator_input(rows, final_dir, mirnas, families, tools, evaluation, out / "evaluator")
    results = run_benchmark(out / "evaluator" / "benchmark.yaml")
    per = pd.read_csv(results / "per_comparison.tsv", sep="\t").merge(comparisons[["experiment_id", "comparison_id", "orientation"]], on="experiment_id")
    per.to_csv(out / "per_comparison.tsv", sep="\t", index=False)

    profiles = pd.read_csv(tables_dir / "mirna_profiles.tsv", sep="\t")
    machinery, coverage = score_machinery(rows, final_dir, profiles, families, tools, root, min_per_state)
    machinery.to_csv(out / "machinery_per_comparison.tsv", sep="\t", index=False)
    coverage.to_csv(out / "machinery_family_coverage.tsv", sep="\t", index=False)
    drop = target_drop(rows, final_dir)
    machinery_seed = pd.read_csv(tables_dir / "machinery_seed_check.tsv", sep="\t")
    kept = figure_machinery(rows, drop, machinery_seed)
    seed = pd.read_csv(tables_dir / "seed_check.tsv", sep="\t")
    drawn_main = figure_main(rows, seed)
    checks = data_checks(rows, seed, machinery_seed, drop, pd.concat([drawn_main, kept]), labels)
    checks.to_csv(out / "data_checks.tsv", sep="\t", index=False)

    single_cell = single_cell_tables(per, labels, drawn_main, tools)
    screens, n_making = machinery_tables(machinery, kept, tools)
    plots = out / "plots"
    plot_bars(single_cell, tools, plots, "single_cell")
    n_groups = len(screens["spearman_shrunk"])
    plot_bars(screens, tools, plots, "machinery", divider=n_making - 0.5,
              groups=(((n_making - 1) / 2, "miRNA-making genes"), (n_making + (n_groups - n_making - 1) / 2, "Repression genes")))
    write_report(out / "REPORT.pdf", drawn_main=drawn_main, single_cell=single_cell, machinery=screens, n_making=n_making, kept=kept, checks=checks,
                 labels=labels, tools=tools, n_genes=int(per.n_genes.iloc[0]), min_per_state=min_per_state, plots=plots)
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
