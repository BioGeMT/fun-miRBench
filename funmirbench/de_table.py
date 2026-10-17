"""DE table reading, gene ID detection, and benchmark FDR derivations."""

import re
from pathlib import Path

import numpy as np
import pandas as pd

_ENSEMBL = re.compile(r"^ENS[A-Z]*G\d+", re.IGNORECASE)
FDR_PLOT_FLOOR = np.nextafter(0, 1)
FDR_DERIVED_COLUMNS = (
    "plot_FDR",
    "benchmark_FDR",
)
# How a missing FDR with a reported PValue (e.g. DESeq2 independent filtering) enters
# the ground truth: "negative" (legacy) treats it as FDR=1.0, so the gene can become a
# negative; "exclude" leaves benchmark_FDR missing so the gene is in neither class.
FDR_NA_POLICIES = ("negative", "exclude")
DEFAULT_FDR_NA_POLICY = "negative"
# Optional DE covariates carried from the DE table into joined tables for stratified
# metrics and diagnostics. They are never predictor scores.
DE_COVARIATE_COLUMNS = (
    "baseMean",
    "logCPM",
    "lfcSE",
    "stat",
    "F",
    "control_mean_normalized_count",
    "treated_mean_normalized_count",
    "normalized_control_mean",
    "normalized_treated_mean",
    "raw_control_mean",
    "raw_treated_mean",
)
DE_STATUS_COLUMN_SUFFIX = "status"
_NON_COVARIATE_PREFIXES = ("score_", "global_rank_", "local_rank_")
DE_COLUMN_ALIASES = {
    "log2foldchange": "logFC",
    "log2_fold_change": "logFC",
    "log_fold_change": "logFC",
    "padj": "FDR",
    "adj_p_val": "FDR",
    "adj_p_value": "FDR",
    "adjusted_p_value": "FDR",
    "qvalue": "FDR",
    "q_value": "FDR",
    "pvalue": "PValue",
    "p_value": "PValue",
    "p_val": "PValue",
}


def _column_key(column) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(column).strip().lower()).strip("_")


def normalize_de_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize common DE-result column aliases to the benchmark schema."""
    rename = {}
    present = {str(column) for column in df.columns}
    for column in df.columns:
        canonical = DE_COLUMN_ALIASES.get(_column_key(column))
        if canonical and canonical not in present:
            rename[column] = canonical
            present.add(canonical)
    if not rename:
        return df
    return df.rename(columns=rename)


def find_gene_id_column(df):
    """Identify which column (or index) holds gene identifiers."""
    for col in df.columns:
        low = str(col).lower().replace(" ", "_")
        if low in ("gene_id", "geneid", "gene_name", "genename"):
            return col

    best_col, best_frac = None, 0.0
    for col in df.columns:
        s = df[col].dropna().astype(str)
        if len(s) == 0:
            continue
        frac = float(s.str.match(_ENSEMBL).mean())
        if frac > best_frac:
            best_frac = frac
            best_col = str(col)
    if best_col is not None and best_frac >= 0.5:
        return best_col

    idx = df.index.astype(str)
    if len(idx) > 0 and float(idx.str.match(_ENSEMBL).mean()) >= 0.5:
        return "__index__"

    return df.columns[0]


def extract_gene_ids(df):
    """Return gene ID values as a list of strings."""
    col = find_gene_id_column(df)
    if col == "__index__":
        return list(df.index.astype(str))
    return df[col].dropna().astype(str).tolist()


def de_covariate_columns(columns) -> list[str]:
    """Return the DE covariate and status columns present in ``columns``."""
    covariates = []
    for column in columns:
        name = str(column)
        if name.startswith(_NON_COVARIATE_PREFIXES):
            continue
        if name in DE_COVARIATE_COLUMNS or _column_key(name).endswith(DE_STATUS_COLUMN_SUFFIX):
            covariates.append(name)
    return covariates


def resolve_fdr_na_policy(value=None) -> str:
    policy = DEFAULT_FDR_NA_POLICY if value is None else str(value).strip().lower()
    if policy not in FDR_NA_POLICIES:
        raise ValueError(
            f"Unsupported fdr_na_policy {value!r}; expected one of {list(FDR_NA_POLICIES)}."
        )
    return policy


def add_fdr_derived_columns(
    df: pd.DataFrame,
    *,
    fdr_na_policy: str = DEFAULT_FDR_NA_POLICY,
) -> pd.DataFrame:
    """Add benchmark-safe FDR columns without changing original FDR/PValue."""
    fdr_na_policy = resolve_fdr_na_policy(fdr_na_policy)
    out = df.copy()
    if "FDR" not in out.columns:
        return out

    fdr = pd.to_numeric(out["FDR"], errors="coerce")
    if "PValue" in out.columns:
        pvalue = pd.to_numeric(out["PValue"], errors="coerce")
    else:
        pvalue = pd.Series(float("nan"), index=out.index)

    fdr_for_plot = fdr.copy()
    fdr_for_plot.loc[fdr.isna()] = 1.0
    fdr_for_plot.loc[fdr.eq(0.0)] = FDR_PLOT_FLOOR

    fdr_for_eval = fdr.copy()
    if fdr_na_policy == "negative":
        fdr_for_eval.loc[fdr.isna() & pvalue.notna()] = 1.0

    out["plot_FDR"] = fdr_for_plot
    out["benchmark_FDR"] = fdr_for_eval
    return out


def read_de_table(path: Path, *, fdr_na_policy: str = DEFAULT_FDR_NA_POLICY) -> pd.DataFrame:
    """Read a DE table TSV, normalizing malformed gene ID headers."""
    df = pd.read_csv(path, sep="\t")
    cols = [str(c).strip() for c in df.columns]
    if cols[0] == "" or cols[0].startswith("Unnamed:"):
        cols[0] = "gene_id"
    df.columns = cols
    if (
        "gene_id" not in df.columns
        and len(df.index) > 0
        and not isinstance(df.index, pd.RangeIndex)
        and float(pd.Series(df.index.astype(str)).str.match(_ENSEMBL).mean()) >= 0.5
    ):
        df = df.reset_index().rename(columns={"index": "gene_id"})
    df = normalize_de_columns(df)
    return add_fdr_derived_columns(df, fdr_na_policy=fdr_na_policy)
