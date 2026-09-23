"""Evaluate joined GT/prediction tables: metrics, plots, reports."""

import math
import os
import pathlib
import textwrap

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, TwoSlopeNorm
from matplotlib.patches import Rectangle
from sklearn.metrics import (
    auc,
    average_precision_score,
    precision_recall_curve,
    roc_auc_score,
    roc_curve,
)

from funmirbench.de_table import FDR_PLOT_FLOOR, add_fdr_derived_columns

SCORE_PREFIX = "score_"
GLOBAL_RANK_PREFIX = "global_rank_"
LOCAL_RANK_PREFIX = "local_rank_"
FDR_AUXILIARY_COLUMNS = (
    "PValue",
    "plot_FDR",
    "benchmark_FDR",
)
FIGURE_DPI = 300
REPORT_PAGE_SIZE = (8.27, 11.69)
PLOT_TITLE_SIZE = 17.0
PLOT_SUBTITLE_SIZE = 12.2
PLOT_AXIS_LABEL_SIZE = 12.0
PLOT_LEGEND_SIZE = 10.5
NEGATIVE_COLOR = "#B8C4D6"
POSITIVE_COLOR = "#D04E4E"
NEUTRAL_COLOR = "#5B6577"
GRID_COLOR = "#D8DEE9"
SCORE_CMAP = ListedColormap(
    ["#F6F7FB", "#C5D7EE", "#7FA8D8", "#2F5D8C", "#17324D"]
)
PREDICTOR_HEATMAP_CMAP = LinearSegmentedColormap.from_list(
    "predictor_heatmap",
    ["#4464D8", "#7EA6E7", "#8CC9B0", "#4AA36B", "#1E6A45"],
)
GT_CMAP = ListedColormap(["#EEF1F6", "#243B53"])
MISSING_COLOR = "#D7DEE8"
CURVE_COLORS = [
    "#1F77B4",
    "#FF7F0E",
    "#2CA02C",
    "#D62728",
    "#9467BD",
    "#17BECF",
    "#E377C2",
    "#BCBD22",
    "#8C564B",
]
TOP_PREDICTION_CDF_N = 100
TOOL_LABELS = {}
TOOL_COLORS = {}
DATASET_LABELS = {}


def _is_publication_tool(tool_id):
    return True


def _publication_tool_ids(tool_ids):
    return [str(tool_id) for tool_id in tool_ids if _is_publication_tool(tool_id)]


def _metric_plot_limits(metric_name):
    if metric_name == "spearman":
        return -1.02, 1.02
    return 0.0, 1.02


def _format_threshold_value(value):
    return str(float(value))


VALID_PERTURBATIONS = ("Overexpression", "Knockout", "Knockdown")
_PERTURBATION_EFFECT_SIGN = {
    "Overexpression": -1.0,
    "Knockout": 1.0,
    "Knockdown": 1.0,
}


def describe_gt_rule(fdr_threshold, effect_threshold, *, markdown=False):
    effect_text = f"perturbation-aware effect > {_format_threshold_value(effect_threshold)}"
    suffix = (
        "(`-logFC` for Overexpression, `+logFC` for Knockout/Knockdown)"
        if markdown
        else "(-logFC for Overexpression, +logFC for Knockout/Knockdown)"
    )
    if fdr_threshold is None:
        if markdown:
            return f"{effect_text.replace('>', '`> ', 1) + '`'} {suffix}"
        return f"{effect_text} {suffix}"
    fdr_text = f"FDR < {_format_threshold_value(fdr_threshold)}"
    if markdown:
        return (
            f"`{fdr_text}` and {effect_text.replace('>', '`> ', 1) + '`'} "
            f"{suffix}"
        )
    return f"{fdr_text} and {effect_text} {suffix}"


def _positive_mask(df, *, fdr_threshold, effect_threshold):
    effect_mask = df["expected_effect"] > float(effect_threshold)
    if fdr_threshold is None:
        return effect_mask
    fdr = df["benchmark_FDR"] if "benchmark_FDR" in df.columns else df["FDR"]
    return effect_mask & (fdr < float(fdr_threshold))


def _usable_gt_row_mask(df, *, fdr_threshold):
    mask = pd.to_numeric(df["logFC"], errors="coerce").notna()
    if fdr_threshold is not None:
        fdr = pd.to_numeric(
            df["benchmark_FDR"] if "benchmark_FDR" in df.columns else df["FDR"],
            errors="coerce",
        )
        mask = mask & fdr.notna() & (fdr >= 0.0) & (fdr <= 1.0)
    return mask


def _filter_usable_gt_rows(df, *, fdr_threshold):
    out = df.copy()
    out["logFC"] = pd.to_numeric(out["logFC"], errors="coerce")
    if "FDR" in out.columns:
        out["FDR"] = pd.to_numeric(out["FDR"], errors="coerce")
        out = add_fdr_derived_columns(out)
    return out.loc[_usable_gt_row_mask(out, fdr_threshold=fdr_threshold)].copy()


def _selection_caption(fdr_threshold):
    if fdr_threshold is None:
        return "selected by expected effect"
    return "selected by FDR and expected effect"


def _sort_heatmap_rows_by_logfc(work):
    sort_cols = ["expected_effect"]
    ascending = [False]
    if "FDR" in work.columns:
        sort_cols.append("FDR")
        ascending.append(True)
    if "gene_id" in work.columns:
        sort_cols.append("gene_id")
        ascending.append(True)
    return work.sort_values(sort_cols, ascending=ascending, kind="mergesort").reset_index(drop=True)


def _emit_log(logger, message):
    if logger is not None:
        logger(message)


def _set_tool_labels(tool_labels=None):
    global TOOL_LABELS
    TOOL_LABELS = {
        str(tool_id): str(label).strip()
        for tool_id, label in (tool_labels or {}).items()
        if str(label).strip()
    }


def _set_dataset_label(dataset_id, label=None):
    global DATASET_LABELS
    dataset_id = str(dataset_id)
    if label:
        DATASET_LABELS[dataset_id] = str(label).strip()
    else:
        DATASET_LABELS.pop(dataset_id, None)


def _set_tool_colors(tool_ids=None):
    global TOOL_COLORS
    TOOL_COLORS = {
        str(tool_id): CURVE_COLORS[index % len(CURVE_COLORS)]
        for index, tool_id in enumerate(tool_ids or [])
    }


def _tool_label(tool_id):
    resolved = TOOL_LABELS.get(str(tool_id))
    if resolved:
        return resolved
    return str(tool_id).replace("_", " ")


def _tool_color(tool_id, *, fallback=POSITIVE_COLOR):
    return TOOL_COLORS.get(str(tool_id), fallback)


def _positive_count_caption(scored_positives, positives_total):
    return f"{int(scored_positives):,}/{int(positives_total):,}"


def _dataset_heading(dataset_id, *, suffix=None):
    label = DATASET_LABELS.get(str(dataset_id), str(dataset_id))
    if suffix:
        return f"{label} | {suffix}"
    return label


def _dataset_caption(dataset_id):
    return DATASET_LABELS.get(str(dataset_id), str(dataset_id).replace("_", " "))


def _wrap_axis_label(text, *, width=14):
    wrapped = textwrap.wrap(str(text), width=width) or [str(text)]
    return "\n".join(wrapped)


def _nice_symmetric_limit(values, *, floor=1.0):
    raw = max(float(np.nanmax(np.abs(values))), float(floor))
    if raw <= 1.0:
        step = 0.25
    elif raw <= 2.0:
        step = 0.5
    elif raw <= 5.0:
        step = 1.0
    elif raw <= 10.0:
        step = 2.0
    else:
        step = 5.0
    return math.ceil(raw / step) * step


def _add_figure_heading(fig, *, title, subtitle, x=0.08, title_y=0.975, subtitle_y=0.93):
    del fig, title, subtitle, x, title_y, subtitle_y


def _add_horizontal_colorbar(fig, *, mappable, anchor_ax, label, ticks=None, height=0.014, pad=0.05):
    pos = anchor_ax.get_position()
    cax = fig.add_axes([pos.x0, pos.y0 - pad, pos.width, height])
    cbar = fig.colorbar(
        mappable,
        cax=cax,
        orientation="horizontal",
        ticks=ticks,
    )
    cbar.set_label(label)
    return cbar


def _style_axes(ax, *, grid_axis="y"):
    ax.set_facecolor("white")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#A0AEC0")
    ax.spines["bottom"].set_color("#A0AEC0")
    ax.tick_params(colors="#3C4858", labelsize=12.0)
    if grid_axis:
        ax.grid(True, axis=grid_axis, color=GRID_COLOR, linewidth=0.8, alpha=0.9)
    ax.set_axisbelow(True)


def _save_figure(fig, out_path):
    out_path = pathlib.Path(out_path)
    fig.savefig(out_path, dpi=FIGURE_DPI, bbox_inches="tight", facecolor="white")
    if out_path.suffix.lower() != ".svg":
        fig.savefig(out_path.with_suffix(".svg"), format="svg", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def _save_pdf_page(pdf, fig):
    fig.set_size_inches(*REPORT_PAGE_SIZE, forward=True)
    fig.patch.set_facecolor("white")
    pdf.savefig(fig, facecolor="white")
    plt.close(fig)


def _ecdf(values):
    ordered = np.sort(np.asarray(values, dtype=float))
    if ordered.size == 0:
        return ordered, ordered
    cumulative = np.arange(1, ordered.size + 1, dtype=float) / float(ordered.size)
    return ordered, cumulative


def _safe_neglog10(series):
    clipped = pd.to_numeric(series, errors="coerce").clip(lower=FDR_PLOT_FLOOR)
    return -clipped.map(math.log10)


def _normalize_perturbation(value):
    text = str(value or "").strip()
    if text.upper() in {"", "NAN", "NONE", "<NA>"}:
        return ""
    canonical = {
        "OVEREXPRESSION": "Overexpression",
        "KNOCKOUT": "Knockout",
        "KNOCKDOWN": "Knockdown",
    }.get(text.upper())
    return canonical or text


def _resolve_perturbation_series(df, *, perturbation=None):
    fallback = _normalize_perturbation(perturbation)
    if "perturbation" in df.columns:
        resolved = df["perturbation"].astype(str).map(_normalize_perturbation)
    else:
        resolved = pd.Series("", index=df.index, dtype="object")
    if fallback:
        resolved = resolved.where(resolved != "", fallback)
    return resolved.fillna("")


def _expected_effect_from_logfc(logfc, perturbations):
    invalid = sorted(
        str(value)
        for value in perturbations.dropna().unique()
        if value not in _PERTURBATION_EFFECT_SIGN
    )
    if invalid:
        expected = ", ".join(VALID_PERTURBATIONS)
        raise ValueError(
            "Missing or invalid perturbation value(s): "
            f"{invalid}. Expected one of: {expected}."
        )
    signs = perturbations.map(_PERTURBATION_EFFECT_SIGN).astype(float)
    return logfc.astype(float) * signs


def _annotate_ground_truth(df, *, perturbation=None):
    out = df.copy()
    out["logFC"] = pd.to_numeric(out["logFC"], errors="coerce")
    if "FDR" in out.columns:
        out = add_fdr_derived_columns(out)
        out["FDR"] = pd.to_numeric(out["FDR"], errors="coerce")
        out["neglog10_FDR"] = _safe_neglog10(out["plot_FDR"])
    else:
        out["FDR"] = np.nan
        out["neglog10_FDR"] = np.nan
    out["resolved_perturbation"] = _resolve_perturbation_series(out, perturbation=perturbation)
    out["expected_effect"] = _expected_effect_from_logfc(
        out["logFC"],
        out["resolved_perturbation"],
    )
    return out


def _rank_scale_scores(series):
    values = series.astype(float)
    ranks = values.rank(method="dense", ascending=True)
    max_rank = ranks.max(skipna=True)
    if pd.isna(max_rank):
        return pd.Series(float("nan"), index=series.index)
    if float(max_rank) <= 1.0:
        return pd.Series(1.0, index=series.index, dtype=float)
    return (ranks - 1.0) / (float(max_rank) - 1.0)


def _tool_id_from_score_col(score_col):
    if score_col.startswith(SCORE_PREFIX):
        return score_col[len(SCORE_PREFIX):]
    return score_col


def _rank_col_for_tool(tool_id, *, prefix=GLOBAL_RANK_PREFIX):
    return f"{prefix}{tool_id}"


def _rank_distribution_specs(frame: pd.DataFrame, *, rank_types=None):
    specs = []
    rank_types = tuple(rank_types) if rank_types is not None else ("local", "global", "score")
    seen = set()
    handlers = {
        "local": LOCAL_RANK_PREFIX,
        "global": GLOBAL_RANK_PREFIX,
        "score": SCORE_PREFIX,
    }
    for rank_type in rank_types:
        prefix = handlers.get(rank_type)
        if prefix is None:
            raise ValueError(f"Unsupported rank distribution type: {rank_type}")
        for column in frame.columns:
            if not column.startswith(prefix):
                continue
            tool_id = column[len(prefix):]
            if tool_id in seen:
                continue
            specs.append((tool_id, column, rank_type))
            seen.add(tool_id)
    return specs


def _top_fraction_mask(series, fraction, *, tie_breaker=None):
    valid = series.notna()
    if not bool(valid.any()):
        return pd.Series(False, index=series.index)
    fraction = float(fraction)
    if fraction <= 0.0:
        return pd.Series(False, index=series.index)
    valid_count = int(valid.sum())
    keep_count = min(valid_count, max(1, int(math.ceil(valid_count * fraction))))
    work = pd.DataFrame({"score": series[valid].astype(float)})
    if tie_breaker is not None:
        work["tie_breaker"] = tie_breaker[valid].astype(str)
        work = work.sort_values(
            ["score", "tie_breaker"], ascending=[False, True], kind="mergesort"
        )
    else:
        work = work.sort_values("score", ascending=False, kind="mergesort")
    selected = pd.Series(False, index=series.index)
    selected.loc[work.index[:keep_count]] = True
    return selected


def _prepare_scored_frame(
    joined, *, score_col, fdr_threshold, effect_threshold, perturbation=None,
):
    required_cols = {"gene_id", "logFC", "FDR", score_col}
    missing = [col for col in required_cols if col not in joined.columns]
    if missing:
        raise ValueError(f"Joined table missing required columns: {missing}")

    keep_cols = ["gene_id", "logFC", "FDR", score_col]
    for optional in ("dataset_id", "mirna", "perturbation", *FDR_AUXILIARY_COLUMNS):
        if optional in joined.columns:
            keep_cols.append(optional)

    keep = _filter_usable_gt_rows(joined[keep_cols], fdr_threshold=fdr_threshold)
    if keep.empty:
        raise ValueError(f"No usable rows remain for {score_col}.")

    total_rows = int(len(keep))
    missing_score_count = int(keep[score_col].isna().sum())
    keep = _annotate_ground_truth(keep, perturbation=perturbation)
    keep["is_positive"] = _positive_mask(
        keep,
        fdr_threshold=fdr_threshold,
        effect_threshold=effect_threshold,
    ).astype(int)
    positives_total = int(keep["is_positive"].sum())

    keep = keep[keep[score_col].notna()].copy()
    if keep.empty:
        raise ValueError(f"No scored rows remain for {score_col}.")

    rows_scored = int(len(keep))
    coverage = float(rows_scored / total_rows) if total_rows else float("nan")
    keep[score_col] = keep[score_col].astype(float)
    positives = int(keep["is_positive"].sum())
    negatives = int(len(keep) - positives)
    if positives == 0:
        raise ValueError(f"No positives remain for {score_col}.")
    if negatives == 0:
        raise ValueError(f"No negatives remain for {score_col}.")
    positive_coverage = float(positives / positives_total) if positives_total else float("nan")
    coverage_info = {
        "rows_total": total_rows,
        "rows_scored": rows_scored,
        "rows_missing_score": missing_score_count,
        "coverage": coverage,
        "positives_total": positives_total,
        "positives_scored": positives,
        "positive_coverage": positive_coverage,
    }
    return keep, coverage_info


__all__ = [name for name in globals() if not name.startswith("__")]
