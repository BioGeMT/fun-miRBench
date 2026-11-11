"""Detailed plots sharing the main report's palette, dimensions and typography."""
from __future__ import annotations

import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from funmirbench.plots import LARGE, SMALL, save, tool_colors, tool_legend, two_panels
from funmirbench.subreport_metrics import UNIVERSE_LABELS, weighted_cdf


def finish(fig, axes, tools: dict[str, str], out: pathlib.Path, stem: str, *, legend=True, palette=None):
    for ax in axes:
        ax.title.set_fontsize(LARGE)
        ax.xaxis.label.set_fontsize(SMALL)
        ax.yaxis.label.set_fontsize(SMALL)
        ax.grid(alpha=0.22, linewidth=0.8)
        ax.set_axisbelow(True)
    if legend:
        colors = palette or tool_colors(list(tools))
        handles = [plt.Line2D([], [], color=colors[t], linewidth=2.5) for t in tools]
        tool_legend(fig, handles, list(tools.values()))
    with plt.rc_context({"font.family": "DejaVu Sans"}):
        save(fig, out, stem)
    return out / f"{stem}.png"


def classification(curves, records, tools, out, stem, universe, palette=None):
    fig, axes = two_panels()
    colors = palette or tool_colors(list(tools))
    for tool in tools:
        c = curves[tool]
        if c.empty:
            continue
        axes[0].plot(np.r_[0, c.recall], np.r_[1, c.precision], color=colors[tool], linewidth=2.5)
        axes[0].axhline(records[tool]["prevalence"], color=colors[tool], linestyle="--", alpha=0.5)
        axes[1].plot(np.r_[0, c.fpr], np.r_[0, c.recall], color=colors[tool], linewidth=2.5)
    axes[1].plot([0, 1], [0, 1], "--", color="#1F2937", linewidth=1)
    for ax in axes:
        ax.set(xlim=(0, 1), ylim=(0, 1.02))
    axes[0].set(xlabel="Recall", ylabel="Precision", title="Precision-recall")
    axes[1].set(xlabel="False positive rate", ylabel="True positive rate", title="ROC")
    if all(c.empty for c in curves.values()):
        for ax in axes:
            ax.text(.5, .5, "Insufficient labelled genes", transform=ax.transAxes, ha="center", fontsize=SMALL)
    return finish(fig, axes, tools, out, stem, palette=palette)


def enrichment(curves, tools, out, stem):
    fig, axes = two_panels()
    colors = tool_colors(list(tools))
    for tool, c in curves.items():
        if c.empty:
            continue
        axes[0].plot(np.r_[0, c.selected / c.selected.iloc[-1]], np.r_[0, c.running_es], color=colors[tool], linewidth=2.5)
        axes[1].plot(c.recall, c.fp_per_tp, color=colors[tool], linewidth=2.5)
    axes[0].axhline(0, color="#1F2937", linestyle="--", linewidth=1)
    axes[0].set(xlabel="Fraction of labelled ranking", ylabel="Running enrichment", title="GSEA-style enrichment", xlim=(0, 1))
    axes[1].set(xlabel="Recall", ylabel="False positives / true positive", title="False-positive burden", xlim=(0, 1))
    if all(c.empty for c in curves.values()):
        for ax in axes:
            ax.text(.5, .5, "Insufficient labelled genes", transform=ax.transAxes, ha="center", fontsize=SMALL)
    return finish(fig, axes, tools, out, stem)


def effect_and_recovery(cdfs, recovery, tools, out):
    fig, axes = two_panels()
    colors = tool_colors(list(tools))
    for tool, c in cdfs.items():
        if not c.empty:
            axes[0].step(c.effect, c.cdf, where="post", color=colors.get(tool, "#1F2937"),
                         linestyle="--" if tool == "background" else "-", linewidth=2)
    axes[0].axvline(0, color="#1F2937", linestyle=":", linewidth=1)
    axes[0].set(xlabel="Oriented shrunken log2FC", ylabel="Cumulative fraction", title="Top-target effect CDF", ylim=(0, 1.02))
    # Keep all observations in the export; the plot follows the central range.
    finite = np.concatenate([c.effect.to_numpy() for c in cdfs.values() if len(c)]) if any(len(c) for c in cdfs.values()) else np.array([])
    if len(finite) and np.ptp(finite):
        lo, hi = np.quantile(finite, [.01, .99])
        pad = max((hi - lo) * .08, .1)
        axes[0].set_xlim(lo - pad, hi + pad)
    for tool in tools:
        c = recovery[recovery.predictor == tool]
        if len(c):
            axes[1].plot(c.budget, c.recall, color=colors[tool], linewidth=2.5, marker="o")
    axes[1].set(xscale="log", xlabel="Prediction budget", ylabel="Recall of responsive genes",
                title="Recovery by prediction budget", ylim=(0, 1.02))
    if recovery.empty:
        axes[1].text(.5, .5, "Insufficient labelled genes", transform=axes[1].transAxes, ha="center", fontsize=SMALL)
    return finish(fig, axes, tools, out, "effects_and_recovery")


def heatmaps(matrices, titles, tools, out, stem):
    fig, axes = two_panels()
    for panel, (ax, matrix, title) in enumerate(zip(axes, matrices, titles)):
        ax.set_title(title)
        ax.grid(False)
        if matrix.empty:
            ax.text(.5, .5, "No usable genes", transform=ax.transAxes, ha="center", fontsize=SMALL)
            continue
        correlation = stem == "predictor_correlations"
        image = ax.imshow(np.ma.masked_invalid(matrix.to_numpy(float)), aspect="auto", interpolation="nearest",
                          cmap="RdBu_r" if correlation else "viridis", vmin=-1 if correlation else 0, vmax=1)
        # Both matrices use the same predictor order; repeating long row labels
        # in the right panel would collide with the left colour bar.
        ax.set_yticks(range(len(matrix)), [tools.get(t, t) for t in matrix.index] if panel == 0 else [], fontsize=SMALL)
        ticks = np.arange(len(matrix.columns)) if correlation else np.arange(0, len(matrix.columns), max(1, len(matrix.columns) // 4))
        ax.set_xticks(ticks, [tools.get(str(matrix.columns[i]), str(matrix.columns[i])) for i in ticks], rotation=90, fontsize=SMALL)
        bar = fig.colorbar(image, ax=ax, shrink=.7)
        bar.ax.tick_params(labelsize=SMALL)
    # Heatmap colours encode values; predictor colours remain consistent in line plots.
    for ax in axes:
        ax.title.set_fontsize(LARGE)
    save(fig, out, stem)
    return out / f"{stem}.png"


def scatter_and_enrichment(table, raw, ranking, curve, tool, label, color, out):
    fig, axes = two_panels()
    usable = table.eligible & np.isfinite(table.oriented_logFC_shrunk)
    ranks = ranking.rank(method="average", pct=True)
    for state, marker, alpha in [("unresolved", ".", .2), ("no_response_above_margin", ".", .35), ("responsive", "o", .65)]:
        m = usable & (table.state_d05 == state)
        axes[0].scatter(ranks[m], table.loc[m, "oriented_logFC_shrunk"], s=12, alpha=alpha,
                        color=color if state == "responsive" else "#6B7280", marker=marker, rasterized=True)
    axes[0].axhline(0, color="#1F2937", linewidth=1, linestyle="--")
    axes[0].set(xlabel="Within-experiment score percentile", ylabel="Oriented shrunken log2FC", title="Score versus response", xlim=(0, 1))
    if len(curve):
        axes[1].plot(np.r_[0, curve.selected / curve.selected.iloc[-1]], np.r_[0, curve.running_es], color=color, linewidth=2.5)
    axes[1].axhline(0, color="#1F2937", linewidth=1, linestyle="--")
    axes[1].set(xlabel="Fraction of labelled ranking", ylabel="Running enrichment", title="GSEA-style enrichment", xlim=(0, 1))
    return finish(fig, axes, {tool: label}, out, "score_and_enrichment", palette={tool: color})
