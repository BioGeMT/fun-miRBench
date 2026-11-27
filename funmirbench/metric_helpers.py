"""Shared ranking metrics used by the benchmark and detailed reports."""
from __future__ import annotations

import numpy as np
from scipy.stats import rankdata, spearmanr


def auroc(score: np.ndarray, positive: np.ndarray) -> float:
    """Mann–Whitney AUROC, assigning half credit to tied positive/negative pairs."""
    score = np.asarray(score, dtype=float)
    positive = np.asarray(positive, dtype=bool)
    n_pos = int(positive.sum())
    n_neg = len(positive) - n_pos
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    ranks = rankdata(score, method="average")
    return float((ranks[positive].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    """Spearman correlation, returning NA for fewer than three or constant values."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) < 3 or len(y) != len(x) or not np.isfinite(x).all() or not np.isfinite(y).all():
        return float("nan")
    if np.ptp(x) == 0 or np.ptp(y) == 0:
        return float("nan")
    return float(spearmanr(x, y).statistic)
