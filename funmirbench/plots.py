"""Plots of the benchmark, OE comparisons, written by funmirbench to plots/ as PNG and SVG.

coverage: genes each tool scores and the responsive and no-response pairs it covers.
oe_all_genes: AUROC of responsive against no-response genes at d = 0.5, Spearman with the shrunken log2FC and the
experiments each tool wins, every tool on all genes (a gene a tool does not score ranks last).
recovery: recall by prediction budget and each experiment's best precision at 50% recall, relative to random.
"""
from __future__ import annotations

import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter, NullLocator

CURVE_COLORS = ["#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#17BECF", "#E377C2", "#BCBD22", "#8C564B"]

METRICS = (("auroc_d05", "AUROC", (0.2, 1.0), (0.2, 0.4, 0.6, 0.8, 1.0)),
           ("spearman_shrunk", "Spearman rho", (-0.5, 0.5), (-0.5, -0.25, 0.0, 0.25, 0.5)))
SLOT, BOX_HEIGHT, WON_WIDTH = 0.75, 2.6, 6.0  # inches per box, boxplot height, width of the experiments-won panel
GAP, WON_GAP = 0.95, 0.3  # inches between the boxplots, and before the experiments-won panel
LEFT, RIGHT, TOP, BOTTOM, LEGEND_GAP = 1.35, 0.2, 1.05, 0.95, 0.3  # inches around the panels
SMALL, LARGE = 17, 22  # font sizes: ticks and counts; panel titles and legend
BAR_HEIGHT, BAR_LABEL_SPACE = 0.45, 0.4  # experiments-won bars, in units of the distance between the two bars
BAR_ALPHA = 0.58
PAIR_SIZE, PAIR_LEFT, PAIR_GAP, PAIR_TOP, PAIR_BOTTOM = (14.0, 6.4), 1.0, 1.3, 1.05, 1.5  # two-panel figures, inches
MARKERS = ("o", "s", "^", "D", "P", "X", "v", "<", ">", "*")
DARK_ALPHA, LIGHT_ALPHA = 0.75, 0.40  # coverage bars


def drawn_tools(tools: pd.DataFrame) -> dict[str, str]:
    """tool_id -> label in file order."""
    return dict(zip(tools.tool_id, tools.official_name))


def tool_colors(tool_ids: list[str]) -> dict[str, str]:
    return {t: CURVE_COLORS[i % len(CURVE_COLORS)] for i, t in enumerate(sorted(tool_ids))}


def experiments_won(per: pd.DataFrame, tool_ids: list[str], metric: str) -> pd.Series:
    """Experiments in which each tool has the highest value; the first in tool order on a tie."""
    wide = per[per.predictor.isin(tool_ids)].pivot(index="experiment_id", columns="predictor", values=metric)
    wins = wide[tool_ids].dropna(how="all").idxmax(axis=1).value_counts()
    return wins.reindex(tool_ids, fill_value=0)


def draw_box(ax, per: pd.DataFrame, tool_ids: list[str], colors: dict, metric: str) -> None:
    values = [per.loc[per.predictor == t, metric].dropna().to_numpy(float) for t in tool_ids]
    positions = np.arange(1, len(tool_ids) + 1)
    box = ax.boxplot(values, positions=positions, widths=0.56, patch_artist=True, showfliers=False,
                     medianprops={"color": "#2D3748", "linewidth": 1.3},
                     whiskerprops={"color": "#6B7280", "linewidth": 1.0},
                     capprops={"color": "#6B7280", "linewidth": 1.0})
    for patch, tool in zip(box["boxes"], tool_ids):
        patch.set(facecolor=colors[tool], edgecolor=colors[tool], alpha=0.28)
    for x, tool, v in zip(positions, tool_ids, values):
        jitter = np.linspace(-0.09, 0.09, v.size) if v.size > 1 else np.zeros(v.size)
        ax.scatter(x + jitter, v, s=42, color=colors[tool], edgecolor="white", linewidth=0.45, alpha=0.75, zorder=3)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", alpha=0.25, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.set_xlim(0.5, len(tool_ids) + 0.5)
    ax.set_xticks([])
    ax.tick_params(axis="y", labelsize=SMALL)


def draw_won(ax, per: pd.DataFrame, tool_ids: list[str], colors: dict) -> None:
    for y, (metric, label, _, _) in enumerate(METRICS):
        counts = experiments_won(per, tool_ids, metric)
        left, total = 0.0, counts.sum()
        for tool, count in counts.items():
            if count:
                ax.barh(y, count / total, left=left, height=BAR_HEIGHT, color=colors[tool], alpha=BAR_ALPHA,
                        edgecolor="white", linewidth=2.0)
                ax.text(left + count / total / 2, y, str(count), ha="center", va="center", fontsize=SMALL,
                        fontweight="bold", color="#111111")
            left += count / total
        ax.text(0.0, y - BAR_HEIGHT / 2 - 0.06, label, ha="left", va="bottom", fontsize=SMALL)
    ax.set_ylim(len(METRICS) - 1 + BAR_HEIGHT / 2, -BAR_HEIGHT / 2 - BAR_LABEL_SPACE)
    ax.set_xlim(0, 1)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines[:].set_visible(False)


def plot_all_genes(per: pd.DataFrame, tools: pd.DataFrame, out: pathlib.Path) -> None:
    labels = drawn_tools(tools)
    tool_ids = list(labels)
    colors = tool_colors(tool_ids)
    box_width = SLOT * len(tool_ids)
    widths = [box_width, box_width, WON_WIDTH]
    gaps = [GAP, WON_GAP]
    fig_w, fig_h = LEFT + sum(widths) + sum(gaps) + RIGHT, TOP + BOX_HEIGHT + BOTTOM
    fig = plt.figure(figsize=(fig_w, fig_h))
    axes, left = [], LEFT
    for i, width in enumerate(widths):
        axes.append(fig.add_axes([left / fig_w, BOTTOM / fig_h, width / fig_w, BOX_HEIGHT / fig_h]))
        left += width + (gaps[i] if i < len(gaps) else 0)
    for ax, (metric, label, ylim, ticks) in zip(axes, METRICS):
        draw_box(ax, per, tool_ids, colors, metric)
        drawn = per.loc[per.predictor.isin(tool_ids), metric]
        ax.set_ylim(min(ylim[0], drawn.min()), max(ylim[1], drawn.max()))  # no value falls outside the panel
        ax.set_yticks(ticks)
        ax.set_title(label, pad=14, fontsize=LARGE, fontweight="bold")
    draw_won(axes[2], per, tool_ids, colors)
    axes[2].set_title("Experiments won", pad=14, fontsize=LARGE, fontweight="bold")
    fig.text(0.3 / fig_w, (BOTTOM + BOX_HEIGHT / 2) / fig_h, "OE", rotation=90, ha="center", va="center",
             fontsize=LARGE, fontweight="bold")
    handles = [plt.Rectangle((0, 0), 1, 1, color=colors[t], alpha=BAR_ALPHA) for t in tool_ids]
    fig.legend(handles, labels.values(), loc="upper center", frameon=False, borderaxespad=0.0, fontsize=LARGE,
               handlelength=1.6, handleheight=1.0, mode="expand", ncol=len(handles),
               bbox_to_anchor=(LEFT / fig_w, 0.0, (sum(widths) + sum(gaps)) / fig_w, (BOTTOM - LEGEND_GAP) / fig_h))
    save(fig, out, "oe_all_genes")


def save(fig, out: pathlib.Path, stem: str) -> None:
    out.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "svg"):
        fig.savefig(out / f"{stem}.{ext}", dpi=300, facecolor="white", bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)


def two_panels():
    w, h = PAIR_SIZE
    panel_w = (w - PAIR_LEFT - PAIR_GAP - RIGHT) / 2
    fig = plt.figure(figsize=PAIR_SIZE)
    axes = [fig.add_axes([left / w, PAIR_BOTTOM / h, panel_w / w, (h - PAIR_TOP - PAIR_BOTTOM) / h])
            for left in (PAIR_LEFT, PAIR_LEFT + panel_w + PAIR_GAP)]
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(labelsize=SMALL)
    return fig, axes


def tool_legend(fig, handles, labels) -> None:
    w, h = PAIR_SIZE
    fig.legend(handles, labels, loc="upper center", frameon=False, borderaxespad=0.0, fontsize=SMALL, mode="expand",
               ncol=len(handles), handlelength=1.6, bbox_to_anchor=(PAIR_LEFT / w, 0.0, (w - PAIR_LEFT - RIGHT) / w,
                                                                    0.65 / h))


def plot_coverage(coverage: pd.DataFrame, tools: pd.DataFrame, out: pathlib.Path) -> None:
    labels = drawn_tools(tools)
    table = coverage.set_index("predictor").loc[list(labels)]
    by_tool = tool_colors(list(labels))
    colors = [by_tool[t] for t in labels]
    x = np.arange(len(labels))
    fig, (genes, pairs) = two_panels()
    share = 100.0 * table.genes_covered / table.genes_measured
    genes.bar(x, share, color=colors, edgecolor=colors, alpha=DARK_ALPHA)
    for xi, value in zip(x, share):
        genes.text(xi, value, f"{value:.1f}%", ha="center", va="bottom", fontsize=SMALL, fontweight="bold")
    genes.set_xticks(x, [f"{r}/{n}" for r, n in zip(table.experiments_retained, table.experiments_total)])
    genes.set_xlabel("Experiments retained", fontsize=SMALL, fontweight="bold")
    genes.set_ylabel("Genes covered (%)", fontsize=SMALL, fontweight="bold")
    genes.set_title("Gene coverage", pad=46, fontsize=LARGE, fontweight="bold")
    low = 100.0 * table.no_response_scored / table.no_response_total
    high = 100.0 * table.responsive_scored / table.responsive_total
    pairs.bar(x, low, color=colors, edgecolor=colors, alpha=DARK_ALPHA)
    pairs.bar(x, (high - low).clip(lower=0), bottom=low, color=colors, edgecolor=colors, alpha=LIGHT_ALPHA)
    for xi, lo, hi in zip(x, low, high):
        pairs.text(xi, max(lo, hi), f"{hi:.1f}%", ha="center", va="bottom", fontsize=SMALL, fontweight="bold")
        pairs.text(xi, lo / 2, f"{lo:.1f}%", ha="center", va="center", fontsize=SMALL, fontweight="bold")
    pairs.set_xticks([])
    pairs.set_ylabel("Pairs covered (%)", fontsize=SMALL, fontweight="bold")
    pairs.set_title("miRNA-gene pair coverage", pad=46, fontsize=LARGE, fontweight="bold")
    shades = [Patch(facecolor="#4B5563", alpha=DARK_ALPHA), Patch(facecolor="#4B5563", alpha=LIGHT_ALPHA)]
    pairs.legend(shades, ["No-response pairs", "Responsive pairs"], loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2,
                 frameon=False, fontsize=SMALL)
    for ax in (genes, pairs):
        ax.set_ylim(0, 112)
        ax.grid(axis="y", alpha=0.25, linewidth=0.8)
        ax.set_axisbelow(True)
    tool_legend(fig, [Patch(facecolor=c, alpha=DARK_ALPHA) for c in colors], list(labels.values()))
    save(fig, out, "coverage")


def plot_recovery(recovery: pd.DataFrame, precision: pd.DataFrame, tools: pd.DataFrame, out: pathlib.Path) -> None:
    labels = drawn_tools(tools)
    colors = tool_colors(list(labels))
    markers = {t: MARKERS[i % len(MARKERS)] for i, t in enumerate(labels)}
    fig, (budget, best) = two_panels()
    means = recovery[recovery.predictor.isin(labels)].groupby(["predictor", "budget"]).recall.mean()
    handles = []
    for tool in labels:
        handles += budget.plot(means[tool].index, means[tool].to_numpy(), color=colors[tool], marker=markers[tool],
                               linewidth=2.5, markersize=9, markeredgecolor="#1F2937", markeredgewidth=0.7)
    budget.set_xscale("log")
    ticks = sorted(recovery.budget.unique())
    budget.set_xticks(ticks[::2])
    budget.xaxis.set_major_formatter(FuncFormatter(lambda value, _pos: f"{int(value):,}"))
    budget.xaxis.set_minor_locator(NullLocator())
    budget.set_ylim(0, 1.02)
    budget.set_xlabel("Prediction budget per experiment", fontsize=SMALL, fontweight="bold")
    budget.set_ylabel("Mean recall of responsive genes", fontsize=SMALL, fontweight="bold")
    budget.set_title("Recovery by prediction budget", pad=14, fontsize=LARGE, fontweight="bold")
    budget.grid(axis="both", alpha=0.22, linewidth=0.8)
    drawn = precision[precision.predictor.isin(labels)].assign(order=lambda d: d.predictor.map(list(labels).index))
    top = drawn.sort_values(["experiment_id", "precision", "false_positives_per_true_positive", "order"],
                            ascending=[True, False, True, True]).groupby("experiment_id").head(1)
    top = top.sort_values("precision_over_random", ascending=False).reset_index(drop=True)
    best.axhline(1.0, color="#1F2937", linewidth=1.0, linestyle=(0, (4, 3)))
    for tool, group in top.groupby("predictor"):
        best.scatter(group.index, group.precision_over_random, s=120, marker=markers[tool], color=colors[tool],
                     edgecolor="#1F2937", linewidth=0.75, zorder=3)
    best.set_yscale("log")
    best.set_ylim(min(0.8, top.precision_over_random.min() * 0.9), max(30, top.precision_over_random.max() * 1.1))
    best.set_yticks([1, 2, 5, 10, 20])
    best.yaxis.set_major_formatter(FuncFormatter(lambda value, _pos: f"{value:g}"))
    best.yaxis.set_minor_locator(NullLocator())
    best.set_xticks([])
    best.set_xlabel("Experiments", fontsize=SMALL, fontweight="bold")
    best.set_ylabel("Best precision / random", fontsize=SMALL, fontweight="bold")
    best.set_title("Best precision at 50% recall", pad=14, fontsize=LARGE, fontweight="bold")
    best.grid(axis="y", alpha=0.22, linewidth=0.8)
    tool_legend(fig, handles, list(labels.values()))
    save(fig, out, "recovery")
