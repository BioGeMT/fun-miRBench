"""Run the fun-miRBench benchmark from a single YAML config: TSV outputs, plots and REPORT.pdf.

For every overexpression (OE) comparison and every tool or reference baseline:
  * all scores are on the same fixed gene list; a gene the tool does not score is ranked last;
  * Spearman between the score and the oriented shrunken log2FC, over genes with a usable estimate, overall and
    averaged within tenths of control expression;
  * AUROC of responsive against no-response genes at d = 0.5 (unresolved left out), plain and within tenths of control
    expression, only where each state has at least min_per_state genes;
  * for tools, a decoy check: the same Spearman for the predictor's scores of unrelated
    benchmark miRNAs (different nt 2-8 seed), giving a z-score of the cognate miRNA against its decoys.
Comparisons are summarised on all comparisons and on unflagged ones, as the mean over comparisons with 95% intervals
from a bootstrap over dependency groups.

Usage:
  uv run funmirbench --config benchmark.yaml
"""
from __future__ import annotations

import argparse
import logging
import pathlib
import shutil

import numpy as np
import pandas as pd
import yaml
from scipy.stats import rankdata, spearmanr

from funmirbench.benchmark_config import (
    build_run_dir_name,
    filter_df,
    resolve_predictor_output_path,
    validate_predictor_output_files,
)
from funmirbench.gene_ids import strip_ensembl_version
from funmirbench.logger import parse_log_level, setup_logging
from funmirbench.plots import plot_all_genes, plot_coverage, plot_recovery
from funmirbench.report_pdf import write_report

logger = logging.getLogger(__name__)
BUDGETS = (10, 25, 50, 100, 250, 500, 1000, 2500, 5000)
RECALL_TARGET = 0.5
TABLE_COLUMNS = ["gene_id", "eligible", "oriented_logFC_shrunk", "control_mean_normalized_count", "state_d05"]


# ---------------------------------------------------------------- metrics
def rank_last(score: pd.Series, genes: pd.Index) -> pd.Series:
    """Scores on the universe, higher = stronger; unscored genes tie below every scored gene."""
    s = score.reindex(genes)
    floor = (s.min() - 1.0) if s.notna().any() else 0.0
    return s.fillna(floor)


def auroc(score: np.ndarray, positive: np.ndarray) -> float:
    n_pos, n_neg = int(positive.sum()), int((~positive).sum())
    if n_pos == 0 or n_neg == 0:
        return np.nan
    ranks = rankdata(score)
    return (ranks[positive].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)


def binned_auroc(score: np.ndarray, positive: np.ndarray, bins: np.ndarray) -> float:
    """Mann-Whitney AUROC over positive/negative pairs within the same bin (weights n_pos * n_neg)."""
    num = den = 0.0
    for b in np.unique(bins):
        m = bins == b
        p, n = int(positive[m].sum()), int((~positive[m]).sum())
        if p and n:
            num += auroc(score[m], positive[m]) * p * n
            den += p * n
    return num / den if den else np.nan


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 3 or np.all(x == x[0]) or np.all(y == y[0]):
        return np.nan
    return float(spearmanr(x, y).statistic)


def binned_spearman(x: np.ndarray, y: np.ndarray, bins: np.ndarray, min_size: int = 20) -> float:
    vals, weights = [], []
    for b in np.unique(bins):
        m = bins == b
        if m.sum() >= min_size:
            r = spearman(x[m], y[m])
            if not np.isnan(r):
                vals.append(r)
                weights.append(m.sum())
    return float(np.average(vals, weights=weights)) if vals else np.nan


def expression_tenths(expression: pd.Series) -> pd.Series:
    ranks = expression.rank(method="first")
    return np.ceil(10 * ranks / len(ranks)).clip(1, 10).astype(int)


def score_metrics(score: pd.Series, table: pd.DataFrame, *, min_per_state: int) -> dict:
    """All metrics of one ranking on one comparison. table is indexed by gene_id (the universe)."""
    s = score.reindex(table.index).to_numpy(float)
    est = table.eligible.to_numpy(bool)
    bins = table.expression_tenth.to_numpy()
    shrunk = table.oriented_logFC_shrunk.to_numpy(float)[est]
    rec = {"spearman_shrunk": spearman(s[est], shrunk),
           "spearman_within_tenths_shrunk": binned_spearman(s[est], shrunk, bins[est])}
    state = table.state_d05.to_numpy()
    labelled = (state == "responsive") | (state == "no_response_above_margin")
    pos = state[labelled] == "responsive"
    n_pos, n_neg = int(pos.sum()), int((~pos).sum())
    rec["n_responsive_d05"], rec["n_no_response_d05"] = n_pos, n_neg
    if n_pos >= min_per_state and n_neg >= min_per_state:
        rec["auroc_d05"] = auroc(s[labelled], pos)
        rec["auroc_within_tenths_d05"] = binned_auroc(s[labelled], pos, bins[labelled])
        rec["auroc_na_reason_d05"] = ""
    else:
        rec["auroc_d05"] = rec["auroc_within_tenths_d05"] = np.nan
        rec["auroc_na_reason_d05"] = f"fewer than {min_per_state} genes in a state"
    return rec


def recall_at_budgets(score: np.ndarray, positive: np.ndarray, genes: np.ndarray) -> dict[int, float]:
    """Share of the positives among the top genes for each budget; ties in gene order."""
    found = np.cumsum(positive[np.lexsort((genes, -score))])
    return {b: found[min(b, len(found)) - 1] / positive.sum() for b in BUDGETS}


def precision_at_recall(score: np.ndarray, positive: np.ndarray, target: float = RECALL_TARGET) -> dict:
    """Precision at the highest score whose genes reach the target recall; tied genes enter together."""
    by_score = pd.DataFrame({"score": score, "positive": positive}).groupby("score").positive.agg(["size", "sum"])
    by_score = by_score.sort_index(ascending=False).cumsum()
    row = by_score[by_score["sum"] / positive.sum() >= target].iloc[0]
    return {"precision": row["sum"] / row["size"], "random_precision": positive.mean(),
            "precision_over_random": row["sum"] / row["size"] / positive.mean(),
            "false_positives_per_true_positive": (row["size"] - row["sum"]) / row["sum"]}


def decoy_z(cognate: np.ndarray, decoys: list[np.ndarray], y: np.ndarray) -> float:
    """z of the Spearman of the cognate scores with y against the Spearman of each decoy's scores; needs two decoys."""
    dvals = np.array([v for v in (spearman(d, y) for d in decoys) if not np.isnan(v)])
    return (spearman(cognate, y) - dvals.mean()) / dvals.std(ddof=1) if len(dvals) > 1 else np.nan


# ---------------------------------------------------------------- inputs
def load_predictions(path: pathlib.Path, direction: str, mirnas: set[str], genes: set[str]) -> dict[str, pd.Series]:
    """miRNA -> gene score (higher = stronger) from a standardized predictor file, strongest score per gene."""
    sign = -1.0 if direction == "lower_is_stronger" else 1.0
    parts = []
    for chunk in pd.read_csv(path, sep="\t", usecols=["Ensembl_ID", "miRNA_Name", "Score"], chunksize=2_000_000):
        chunk = chunk[chunk.miRNA_Name.isin(mirnas)]
        if len(chunk):
            chunk = chunk.assign(gene_id=chunk.Ensembl_ID.astype(str).map(strip_ensembl_version),
                                 score=sign * pd.to_numeric(chunk.Score, errors="coerce"))
            parts.append(chunk.loc[chunk.gene_id.isin(genes) & chunk.score.notna(), ["miRNA_Name", "gene_id", "score"]])
    if not parts:
        return {}
    allp = pd.concat(parts).groupby(["miRNA_Name", "gene_id"]).score.max()
    return {m: allp.loc[m] for m in allp.index.get_level_values(0).unique()}


def group_bootstrap(values: pd.Series, groups: pd.Series, *, n: int, rng: np.random.Generator) -> tuple[float, float]:
    """95% percentile interval of the mean over comparisons, resampling dependency groups."""
    ok = values.notna()
    values, groups = values[ok], groups[ok]
    labels = groups.unique()
    if len(labels) < 2:
        return np.nan, np.nan
    by_group = {g: values[groups == g].to_numpy() for g in labels}
    means = []
    for _ in range(n):
        pick = rng.choice(labels, size=len(labels), replace=True)
        means.append(np.concatenate([by_group[g] for g in pick]).mean())
    return float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5))


# ---------------------------------------------------------------- main
def run_benchmark(config_path: pathlib.Path) -> pathlib.Path:
    config_path = config_path.expanduser().resolve()
    with config_path.open(encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    root = config_path.parent
    input_dir = resolve_predictor_output_path(config["input_dir"], root)
    if not (input_dir / "comparisons.tsv").is_file():
        raise FileNotFoundError(f"no benchmark tables in {input_dir}; run `uv run funmirbench-download-data`")
    evaluation = config.get("evaluation") or {}
    min_per_state = int(evaluation.get("min_per_state", 30))
    bootstrap = int(evaluation.get("bootstrap", 2000))
    seed = int(evaluation.get("seed", 0))
    comps = pd.read_csv(input_dir / "comparisons.tsv", sep="\t", keep_default_na=False)
    skipped = comps.experiment_id[comps.perturbation != "OE"].tolist()
    if skipped:
        logger.warning("Only OE comparisons are scored; skipping %d others: %s", len(skipped), ", ".join(skipped))
    comps = comps[comps.perturbation == "OE"]
    if config.get("experiments"):
        comps = filter_df(comps, config["experiments"])
    set_aside = comps[comps.action == "set_aside"]
    comps = comps[comps.action != "set_aside"]
    mirnas = pd.read_csv(resolve_predictor_output_path(config["mirnas_tsv"], root), sep="\t").set_index("mirna_name")
    tools = pd.read_csv(resolve_predictor_output_path(config["predictions_tsv"], root), sep="\t")
    if config.get("predictors"):
        tools = filter_df(tools, config["predictors"])
    validate_predictor_output_files({r["tool_id"]: r for r in tools.to_dict("records")}, root)
    out = resolve_predictor_output_path(config.get("out_dir", "results"), root) / build_run_dir_name(
        experiments=None, tool_ids=None, eval_cfg=None)
    out.mkdir(parents=True, exist_ok=True)
    shutil.copy2(config_path, out / "benchmark_config.yaml")
    final = pd.read_csv(input_dir / "experiments.tsv", sep="\t").set_index("experiment_id")
    missing = set(comps.experiment_id) - set(final.index)
    if missing:
        raise ValueError(f"comparisons without a benchmark table: {sorted(missing)}")
    if (set(comps.mirna_name) - set(mirnas.index)):
        raise ValueError("every comparison miRNA needs a sequence in mirnas_tsv")

    tables = {}
    for exp in comps.experiment_id:
        t = pd.read_csv(input_dir / final.loc[exp, "table"], sep="\t", usecols=TABLE_COLUMNS).set_index("gene_id")
        t["eligible"] = t.eligible.astype(bool)
        t["expression_tenth"] = expression_tenths(t.control_mean_normalized_count.fillna(0.0))
        tables[exp] = t
    universe = next(iter(tables.values())).index
    if any(not t.index.equals(universe) for t in tables.values()):
        raise ValueError("benchmark tables must share the same gene list in the same order")
    logger.info("%d comparisons, %d genes", len(tables), len(universe))

    seeds = mirnas.mirna_sequence.str.upper().str.replace("T", "U").str[1:8]
    preds = {}
    for r in tools.itertuples(index=False):
        path = resolve_predictor_output_path(r.predictor_output_path, root)
        logger.info("loading %s", r.tool_id)
        preds[r.tool_id] = load_predictions(path, r.score_direction, set(mirnas.index), set(universe))

    rng = np.random.default_rng(seed)
    rows, decoy_rows, pair_rows, recovery_rows, precision_rows = [], [], [], [], []
    covered, measured = {tool: set() for tool in preds}, set()
    for c in comps.itertuples(index=False):
        t = tables[c.experiment_id]
        expr = t.control_mean_normalized_count.fillna(0.0)
        rankings = {
            "random": pd.Series(rng.random(len(universe)), index=universe),
            "expression_high": expr,
            "expression_low": -expr,
        }
        kinds = {k: "reference" for k in rankings}
        decoys = [m for m in mirnas.index if m != c.mirna_name and seeds[m] != seeds[c.mirna_name]]
        decoy_scores = {}
        for tool, by_mirna in preds.items():
            scored = by_mirna.get(c.mirna_name, pd.Series(dtype=float))
            rankings[tool] = rank_last(scored, universe)
            kinds[tool] = "tool"
            decoy_scores[tool] = [rank_last(by_mirna[d], universe) for d in decoys if d in by_mirna]
        est = t.eligible.to_numpy(bool)
        y = t.oriented_logFC_shrunk.to_numpy(float)[est]
        for name, scores in decoy_scores.items():
            z = decoy_z(rankings[name].to_numpy(float)[est], [d.to_numpy(float)[est] for d in scores], y)
            decoy_rows.append({"experiment_id": c.experiment_id, "predictor": name, "decoy_z_shrunk": z})
        for name, score in rankings.items():
            rec = {"experiment_id": c.experiment_id, "action": c.action,
                   "dependency_group": c.dependency_group, "predictor": name, "kind": kinds[name],
                   "n_genes": len(universe),
                   "coverage": float(preds[name].get(c.mirna_name, pd.Series(dtype=float)).size / len(universe))
                   if name in preds else 1.0}
            rec.update(score_metrics(score, t, min_per_state=min_per_state))
            rows.append(rec)
        state = t.state_d05.to_numpy()
        responsive, no_response = state == "responsive", state == "no_response_above_margin"
        labelled = responsive | no_response
        enough = min(responsive.sum(), no_response.sum()) >= min_per_state
        measured.update(universe[est])
        for tool, by_mirna in preds.items():
            scored = universe.isin(by_mirna.get(c.mirna_name, pd.Series(dtype=float)).index)
            covered[tool].update(universe[scored & est])
            pair_rows.append({"experiment_id": c.experiment_id, "predictor": tool,
                              "responsive_scored": int((scored & responsive).sum()),
                              "responsive_total": int(responsive.sum()),
                              "no_response_scored": int((scored & no_response).sum()),
                              "no_response_total": int(no_response.sum())})
            if enough:
                s = rankings[tool].to_numpy(float)[labelled]
                recalls = recall_at_budgets(s, responsive[labelled], universe[labelled].to_numpy())
                recovery_rows += [{"experiment_id": c.experiment_id, "predictor": tool, "budget": b, "recall": r}
                                  for b, r in recalls.items()]
                precision_rows.append({"experiment_id": c.experiment_id, "predictor": tool,
                                       **precision_at_recall(s, responsive[labelled])})
        logger.info("scored %s", c.experiment_id)

    per = pd.DataFrame(rows).merge(pd.DataFrame(decoy_rows), on=["experiment_id", "predictor"], how="left")
    per.to_csv(out / "per_comparison.tsv", sep="\t", index=False)

    metrics = ["spearman_shrunk", "spearman_within_tenths_shrunk", "auroc_d05", "auroc_within_tenths_d05",
               "decoy_z_shrunk", "coverage"]
    summary = []
    boot_rng = np.random.default_rng(seed + 1)
    for subset, sub in (("all", per), ("unflagged", per[per.action == "keep"])):
        for name, g in sub.groupby("predictor"):
            gi = g.set_index("experiment_id")
            for metric in metrics:
                v = gi[metric]
                lo, hi = group_bootstrap(v, gi.dependency_group, n=bootstrap, rng=boot_rng)
                summary.append({"subset": subset, "predictor": name, "kind": gi.kind.iloc[0],
                                "metric": metric, "mean": v.mean(), "ci95_low": lo, "ci95_high": hi,
                                "n_comparisons": int(v.notna().sum()),
                                "n_groups": int(gi.dependency_group[v.notna()].nunique())})
    summary = pd.DataFrame(summary)
    summary.to_csv(out / "summary.tsv", sep="\t", index=False)
    pairs = pd.DataFrame(pair_rows)
    retained = pairs[(pairs.responsive_scored > 0) & (pairs.no_response_scored > 0)].groupby("predictor").size()
    coverage = pairs.groupby("predictor")[["responsive_scored", "responsive_total", "no_response_scored",
                                           "no_response_total"]].sum()
    coverage.insert(0, "genes_covered", pd.Series({tool: len(genes) for tool, genes in covered.items()}))
    coverage.insert(1, "genes_measured", len(measured))
    coverage.insert(2, "experiments_retained", retained.reindex(coverage.index, fill_value=0))
    coverage.insert(3, "experiments_total", len(comps))
    coverage = coverage.reset_index()
    recovery, precision = pd.DataFrame(recovery_rows), pd.DataFrame(precision_rows)
    coverage.to_csv(out / "coverage.tsv", sep="\t", index=False)
    recovery.to_csv(out / "recovery.tsv", sep="\t", index=False)
    precision.to_csv(out / "precision_at_recall.tsv", sep="\t", index=False)
    plots = out / "plots"
    plot_coverage(coverage, tools, plots)
    plot_all_genes(per, tools, plots)
    plot_recovery(recovery, precision, tools, plots)
    write_report(out / "REPORT.pdf", comps=comps, set_aside=set_aside, summary=summary, tools=tools,
                 n_genes=len(universe), min_per_state=min_per_state, plots=plots)
    logger.info("wrote %s", out)
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the fun-miRBench benchmark.")
    parser.add_argument("--config", type=pathlib.Path, required=True)
    parser.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"])
    args = parser.parse_args()
    setup_logging(parse_log_level(args.log_level))
    run_benchmark(args.config)


if __name__ == "__main__":
    main()
