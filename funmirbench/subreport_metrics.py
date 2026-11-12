"""Per-experiment metrics on the v2 labels; score ties enter as whole blocks.

FPS is the full fixed gene list, IPS the per-comparison intersection of scored
pairs, and own_scored each tool's scored pairs. Unresolved genes are never negatives.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.integrate import trapezoid

from funmirbench.metric_helpers import auroc, spearman

UNIVERSES = ("fps", "ips", "own_scored")
UNIVERSE_LABELS = {"fps": "Full pair set", "ips": "Intersection pair set", "own_scored": "Own scored pairs"}


def masks(scores: pd.DataFrame, tool: str, universe: str) -> np.ndarray:
    if universe == "fps":
        return np.ones(len(scores), dtype=bool)
    if universe == "ips":
        return scores.notna().all(axis=1).to_numpy()
    if universe == "own_scored":
        return scores[tool].notna().to_numpy()
    raise ValueError(f"Unknown pair universe: {universe}")


def threshold_curves(score: np.ndarray, positive: np.ndarray) -> pd.DataFrame:
    """One point per distinct score, descending; empty unless both states exist."""
    if len(score) == 0 or not positive.any() or positive.all():
        return pd.DataFrame(columns=["threshold", "selected", "tp", "fp", "recall", "precision", "fpr", "fp_per_tp", "running_es"])
    grouped = pd.DataFrame({"score": score, "positive": positive}).groupby("score").positive.agg(["size", "sum"])
    grouped = grouped.sort_index(ascending=False)
    selected, tp = grouped["size"].cumsum(), grouped["sum"].cumsum()
    fp = selected - tp
    recall, fpr = tp / positive.sum(), fp / (~positive).sum()
    return pd.DataFrame({"threshold": grouped.index, "selected": selected.to_numpy(), "tp": tp.to_numpy(),
                         "fp": fp.to_numpy(), "recall": recall.to_numpy(), "precision": (tp / selected).to_numpy(),
                         "fpr": fpr.to_numpy(), "fp_per_tp": (fp / tp.replace(0, np.nan)).to_numpy(),
                         "running_es": (recall - fpr).to_numpy()})


def metrics(table: pd.DataFrame, raw: pd.DataFrame, ranking: pd.Series, tool: str, universe: str,
            min_per_state: int) -> tuple[dict, pd.DataFrame]:
    selected = masks(raw, tool, universe)
    states = table.state_d05.to_numpy()
    labelled = selected & np.isin(states, ["responsive", "no_response_above_margin"])
    positive = states[labelled] == "responsive"
    n_pos, n_neg = int(positive.sum()), int((~positive).sum())
    effects = table.oriented_logFC_shrunk.to_numpy(float)
    continuous = selected & table.eligible.to_numpy(bool) & np.isfinite(effects)
    usable = continuous & np.isfinite(ranking.to_numpy(float))
    x, y = ranking.to_numpy(float)[usable], effects[usable]
    rho = spearman(x, y)
    variable = np.isfinite(rho)
    if len(x) < 3:
        continuous_note = "Needs at least 3 eligible genes with finite scores and effects"
    elif not variable:
        continuous_note = "Score or oriented effect is constant"
    else:
        continuous_note = ""
    n_scored = int((selected & raw[tool].notna().to_numpy()).sum())
    rec = {"predictor": tool, "universe": universe, "n_genes": int(selected.sum()), "n_usable": int(usable.sum()),
           "n_scored": n_scored,
           "n_responsive": n_pos, "n_no_response": n_neg,
           "n_unresolved": int((selected & (states == "unresolved")).sum()),
           "prevalence": n_pos / (n_pos + n_neg) if n_pos + n_neg else np.nan,
           "spearman": rho,
           "spearman_na_reason": continuous_note,
           "average_precision": np.nan, "pr_auc": np.nan, "auroc": np.nan, "enrichment_score": np.nan,
           "classification_note": ""}
    if min(n_pos, n_neg) < max(1, min_per_state):
        rec["classification_note"] = f"Needs at least {max(1, min_per_state)} genes in each labelled state"
        return rec, threshold_curves(np.array([]), np.array([], dtype=bool))
    curve = threshold_curves(ranking.to_numpy(float)[labelled], positive)
    recall, precision = curve.recall.to_numpy(), curve.precision.to_numpy()
    rec["average_precision"] = float(np.sum(np.diff(np.r_[0, recall]) * precision))
    rec["pr_auc"] = float(trapezoid(np.r_[1, precision], np.r_[0, recall]))
    rec["auroc"] = auroc(ranking.to_numpy(float)[labelled], positive)
    rec["enrichment_score"] = float(curve.running_es.iloc[np.abs(curve.running_es.to_numpy()).argmax()])
    return rec, curve


def top_weights(ranking: pd.Series, budget: int) -> pd.Series:
    """Expected inclusion in the top budget over all orders within a tied block."""
    n = min(budget, len(ranking))
    weights = pd.Series(0.0, index=ranking.index)
    if n == 0:
        return weights
    cutoff = ranking.nlargest(n).iloc[-1]
    above, tied = ranking > cutoff, ranking == cutoff
    weights.loc[above] = 1.0
    weights.loc[tied] = (n - int(above.sum())) / int(tied.sum())
    return weights


def weighted_cdf(values: np.ndarray, weights: np.ndarray | None = None) -> pd.DataFrame:
    weights = np.ones(len(values)) if weights is None else weights
    keep = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    if not keep.any():
        return pd.DataFrame(columns=["effect", "cdf", "weight"])
    grouped = pd.DataFrame({"effect": values[keep], "weight": weights[keep]}).groupby("effect").weight.sum()
    return pd.DataFrame({"effect": grouped.index, "cdf": (grouped.cumsum() / grouped.sum()).to_numpy(),
                         "weight": grouped.to_numpy()})
