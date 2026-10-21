"""FuNmiRBench v2 evaluation on the final three-state tables (no plots; TSV outputs).

For every comparison and every tool or reference baseline:
  * all scores are on the same fixed gene list; a gene the tool does not score is ranked last;
  * Spearman between the score and the oriented log2FC (raw and shrunken), over genes with a usable estimate,
    overall and averaged within tenths of control expression;
  * AUROC of responsive against no-response genes (unresolved left out), plain and within tenths of control
    expression, only where each state has at least --min-per-state genes; for the states of the table (state_*)
    and for the raw states kept as a check (state_raw_*, metrics named *_raw_*);
  * for tools and the two seed baselines, a decoy check: the same Spearman for the predictor's scores of
    unrelated benchmark miRNAs (different nt 2-8 seed), giving a z-score of the cognate miRNA against its
    decoys, on the raw (decoy_z) and shrunken (decoy_z_shrunk) oriented log2FC.
Comparisons are summarised per track (OE, KO, and pooled over both), on all comparisons and on unflagged ones, as
the mean over comparisons with 95% intervals from a bootstrap over dependency groups; tools are also compared
with the seed-plus-expression reference on the same comparisons.

Usage:
  uv run funmirbench-v2-evaluate --final-dir <final-v2> --utr-fasta <MANE 3'UTR fasta> --out <dir>
"""
from __future__ import annotations

import argparse
import json
import logging
import pathlib
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

from funmirbench.gene_ids import strip_ensembl_version

logger = logging.getLogger(__name__)
REPO = pathlib.Path(__file__).resolve().parents[1]
MARGINS = ("d05", "d10")
STATE_SETS = {"": "state_", "raw_": "state_raw_"}  # metric tag prefix -> state column prefix
TRACKS = ("OE", "KO", "pooled")
REFERENCE = "seed_plus_expression"
TABLE_COLUMNS = ["gene_id", "eligible", "oriented_logFC", "oriented_logFC_shrunk", "control_mean_normalized_count",
                 *[f"{column}{m}" for column in STATE_SETS.values() for m in MARGINS]]
COMPLEMENT = str.maketrans("ACGU", "UGCA")


# ---------------------------------------------------------------- seed sites
def read_fasta_by_gene(path: pathlib.Path) -> dict[str, str]:
    """FASTA with headers '>ENSG...|...'; returns gene_id -> sequence (RNA alphabet, upper case)."""
    seqs, gene, parts = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if gene is not None:
                    seqs[gene] = "".join(parts).upper().replace("T", "U")
                gene, parts = strip_ensembl_version(line[1:].split("|")[0]), []
            elif line:
                parts.append(line)
    if gene is not None:
        seqs[gene] = "".join(parts).upper().replace("T", "U")
    return seqs


def seed_site_counts(utr: str, mirna: str) -> tuple[int, int, int, int]:
    """Canonical 3'UTR sites: (8mer, 7mer-m8, 7mer-A1, 6mer), each site counted once in that hierarchy."""
    mirna = mirna.upper().replace("T", "U")
    s6 = mirna[1:7].translate(COMPLEMENT)[::-1]          # pairs with nt 2-7
    m8 = mirna[7].translate(COMPLEMENT)                 # pairs with nt 8, 5' of s6 in the target
    counts = [0, 0, 0, 0]
    start = utr.find(s6)
    while start >= 0:
        has8 = start > 0 and utr[start - 1] == m8
        hasA = start + 6 < len(utr) and utr[start + 6] == "A"
        counts[0 if has8 and hasA else 1 if has8 else 2 if hasA else 3] += 1
        start = utr.find(s6, start + 1)
    return tuple(counts)


def seed_score(utrs: dict[str, str], genes: pd.Index, mirna_sequence: str) -> pd.Series:
    """3 x 8mer + 2 x 7mer-m8 + 1 x 7mer-A1 + 0.5 x 6mer; 0 for genes without a 3'UTR."""
    weights = np.array([3.0, 2.0, 1.0, 0.5])
    values = [float(np.dot(weights, seed_site_counts(utrs[g], mirna_sequence))) if g in utrs else 0.0 for g in genes]
    return pd.Series(values, index=genes)


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
    raw = table.oriented_logFC.to_numpy(float)[est]
    shrunk = table.oriented_logFC_shrunk.to_numpy(float)[est]
    rec = {"spearman": spearman(s[est], raw),
           "spearman_shrunk": spearman(s[est], shrunk),
           "spearman_within_tenths": binned_spearman(s[est], raw, bins[est]),
           "spearman_within_tenths_shrunk": binned_spearman(s[est], shrunk, bins[est])}
    for prefix, column in STATE_SETS.items():
        for m in MARGINS:
            tag = f"{prefix}{m}"
            state = table[f"{column}{m}"].to_numpy()
            labelled = (state == "responsive") | (state == "no_response_above_margin")
            pos = state[labelled] == "responsive"
            n_pos, n_neg = int(pos.sum()), int((~pos).sum())
            rec[f"n_responsive_{tag}"], rec[f"n_no_response_{tag}"] = n_pos, n_neg
            if n_pos >= min_per_state and n_neg >= min_per_state:
                rec[f"auroc_{tag}"] = auroc(s[labelled], pos)
                rec[f"auroc_within_tenths_{tag}"] = binned_auroc(s[labelled], pos, bins[labelled])
                rec[f"auroc_na_reason_{tag}"] = ""
            else:
                rec[f"auroc_{tag}"] = rec[f"auroc_within_tenths_{tag}"] = np.nan
                rec[f"auroc_na_reason_{tag}"] = f"fewer than {min_per_state} genes in a state"
    return rec


def decoy_record(cognate: np.ndarray, decoys: list[np.ndarray], y: np.ndarray) -> dict:
    """Spearman of the cognate scores and of each decoy's scores with y; z of the cognate against the decoys."""
    rho = spearman(cognate, y)
    dvals = np.array([v for v in (spearman(d, y) for d in decoys) if not np.isnan(v)])
    return {"spearman": rho, "n_decoys": len(dvals), "decoy_mean": dvals.mean() if len(dvals) else np.nan,
            "decoy_sd": dvals.std(ddof=1) if len(dvals) > 1 else np.nan,
            "decoy_z": ((rho - dvals.mean()) / dvals.std(ddof=1)) if len(dvals) > 1 else np.nan,
            "share_decoys_at_least_cognate": (dvals >= rho).mean() if len(dvals) else np.nan}


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
def evaluate(args: argparse.Namespace) -> None:
    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    comps = pd.read_csv(args.comparisons or args.final_dir / "comparisons.tsv", sep="\t", keep_default_na=False)
    comps = comps[comps.action != "set_aside"]
    mirnas = pd.read_csv(args.mirnas, sep="\t").set_index("mirna_name")
    tools = pd.read_csv(args.predictions, sep="\t")
    if args.tools:
        tools = tools[tools.tool_id.isin(args.tools)]
    final = pd.read_csv(args.final_dir / "experiments.tsv", sep="\t").set_index("experiment_id")
    missing = set(comps.experiment_id) - set(final.index)
    if missing:
        raise ValueError(f"comparisons without a final table: {sorted(missing)}")
    if (set(comps.mirna_name) - set(mirnas.index)):
        raise ValueError("every comparison miRNA needs a sequence in --mirnas")

    tables = {}
    for exp in comps.experiment_id:
        t = pd.read_csv(args.final_dir / final.loc[exp, "table"], sep="\t", usecols=TABLE_COLUMNS).set_index("gene_id")
        t["eligible"] = t.eligible.astype(bool)
        t["expression_tenth"] = expression_tenths(t.control_mean_normalized_count.fillna(0.0))
        tables[exp] = t
    universe = next(iter(tables.values())).index
    if any(not t.index.equals(universe) for t in tables.values()):
        raise ValueError("final tables must share the same gene list in the same order")
    utrs = read_fasta_by_gene(args.utr_fasta)
    logger.info("%d comparisons, %d genes, %d 3'UTRs", len(tables), len(universe), sum(g in utrs for g in universe))

    seeds = mirnas.mirna_sequence.str.upper().str.replace("T", "U").str[1:8]
    preds = {}
    for r in tools.itertuples(index=False):
        path = REPO / r.predictor_output_path
        if not path.exists():
            raise FileNotFoundError(f"{r.tool_id}: {path} (run funmirbench-download-data)")
        logger.info("loading %s", r.tool_id)
        preds[r.tool_id] = load_predictions(path, r.score_direction, set(mirnas.index), set(universe))

    sites = {}

    def site_score(mirna: str) -> pd.Series:
        if mirna not in sites:
            sites[mirna] = seed_score(utrs, universe, mirnas.loc[mirna, "mirna_sequence"])
        return sites[mirna]

    rng = np.random.default_rng(args.seed)
    rows, decoy_rows = [], []
    for c in comps.itertuples(index=False):
        t = tables[c.experiment_id]
        expr = t.control_mean_normalized_count.fillna(0.0)
        expr_pct = expr.rank(pct=True)
        site = site_score(c.mirna_name)
        rankings = {
            "random": pd.Series(rng.random(len(universe)), index=universe),
            "expression_high": expr,
            "expression_low": -expr,
            "seed_site": site,
            REFERENCE: site * 10 + expr_pct,
        }
        kinds = {k: "reference" for k in rankings}
        decoys = [m for m in mirnas.index if m != c.mirna_name and seeds[m] != seeds[c.mirna_name]]
        decoy_scores = {"seed_site": [site_score(m) for m in decoys],
                        REFERENCE: [site_score(m) * 10 + expr_pct for m in decoys]}
        for tool, by_mirna in preds.items():
            scored = by_mirna.get(c.mirna_name, pd.Series(dtype=float))
            rankings[tool] = rank_last(scored, universe)
            kinds[tool] = "tool"
            decoy_scores[tool] = [rank_last(by_mirna[d], universe) for d in decoys if d in by_mirna]
        est = t.eligible.to_numpy(bool)
        for target, column in (("raw", "oriented_logFC"), ("shrunk", "oriented_logFC_shrunk")):
            y = t[column].to_numpy(float)[est]
            for name, scores in decoy_scores.items():
                rec = decoy_record(rankings[name].to_numpy(float)[est], [d.to_numpy(float)[est] for d in scores], y)
                decoy_rows.append({"experiment_id": c.experiment_id, "predictor": name, "target": target, **rec})
        for name, score in rankings.items():
            rec = {"experiment_id": c.experiment_id, "perturbation": c.perturbation, "action": c.action,
                   "dependency_group": c.dependency_group, "predictor": name, "kind": kinds[name],
                   "n_genes": len(universe),
                   "coverage": float(preds[name].get(c.mirna_name, pd.Series(dtype=float)).size / len(universe))
                   if name in preds else 1.0}
            rec.update(score_metrics(score, t, min_per_state=args.min_per_state))
            rows.append(rec)
        logger.info("scored %s", c.experiment_id)

    per = pd.DataFrame(rows)
    decoy = pd.DataFrame(decoy_rows, columns=["experiment_id", "predictor", "target", "spearman", "n_decoys",
                                              "decoy_mean", "decoy_sd", "decoy_z", "share_decoys_at_least_cognate"])
    z = decoy.pivot(index=["experiment_id", "predictor"], columns="target", values="decoy_z")
    z = z.rename(columns={"raw": "decoy_z", "shrunk": "decoy_z_shrunk"}).reset_index()
    per = per.merge(z, on=["experiment_id", "predictor"], how="left")
    per.to_csv(out / "per_comparison.tsv", sep="\t", index=False)
    decoy.to_csv(out / "decoy.tsv", sep="\t", index=False)

    tags = [f"{prefix}{m}" for prefix in STATE_SETS for m in MARGINS]
    metrics = ["spearman", "spearman_shrunk", "spearman_within_tenths", "spearman_within_tenths_shrunk",
               *[f"auroc_{t}" for t in tags], *[f"auroc_within_tenths_{t}" for t in tags], "decoy_z",
               "decoy_z_shrunk", "coverage"]
    summary, paired = [], []
    boot_rng = np.random.default_rng(args.seed + 1)
    for subset, sub in (("all", per), ("unflagged", per[per.action == "keep"])):
        for track in TRACKS:
            in_track = sub if track == "pooled" else sub[sub.perturbation == track]
            for name, g in in_track.groupby("predictor"):
                ref = in_track[in_track.predictor == REFERENCE].set_index("experiment_id")
                gi = g.set_index("experiment_id")
                for metric in metrics:
                    v = gi[metric]
                    lo, hi = group_bootstrap(v, gi.dependency_group, n=args.bootstrap, rng=boot_rng)
                    summary.append({"subset": subset, "track": track, "predictor": name, "kind": gi.kind.iloc[0],
                                    "metric": metric, "mean": v.mean(), "ci95_low": lo, "ci95_high": hi,
                                    "n_comparisons": int(v.notna().sum()),
                                    "n_groups": int(gi.dependency_group[v.notna()].nunique())})
                    if name != REFERENCE and metric != "coverage" and metric in ref and ref[metric].notna().any():
                        d = (v - ref[metric].reindex(v.index)).dropna()
                        lo, hi = group_bootstrap(d, gi.dependency_group.reindex(d.index), n=args.bootstrap,
                                                 rng=boot_rng)
                        paired.append({"subset": subset, "track": track, "predictor": name, "metric": metric,
                                       "mean_difference_vs_" + REFERENCE: d.mean() if len(d) else np.nan,
                                       "ci95_low": lo, "ci95_high": hi, "n_comparisons": len(d),
                                       "n_groups": int(gi.dependency_group.reindex(d.index).nunique())})
    pd.DataFrame(summary).to_csv(out / "summary.tsv", sep="\t", index=False)
    pd.DataFrame(paired).to_csv(out / "paired_vs_reference.tsv", sep="\t", index=False)
    (out / "run_info.json").write_text(json.dumps({
        "built": datetime.now(timezone.utc).isoformat(timespec="seconds"), "final_dir": str(args.final_dir),
        "comparisons": int(len(comps)), "genes": int(len(universe)), "tools": list(preds),
        "min_per_state": args.min_per_state, "bootstrap": args.bootstrap, "seed": args.seed,
        "reference": REFERENCE, "tracks": list(TRACKS), "states": {"state_*": "main", "state_raw_*": "raw check"},
        "decoys": "benchmark miRNAs (--mirnas) with a different nt 2-8 seed"}, indent=2) + "\n")
    logger.info("wrote %s", out)


def parse_args(argv=None) -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--final-dir", type=pathlib.Path, required=True, help="output of build_final_data.py")
    ap.add_argument("--utr-fasta", type=pathlib.Path, required=True, help="MANE Select 3'UTR per gene")
    ap.add_argument("--out", type=pathlib.Path, required=True)
    ap.add_argument("--comparisons", type=pathlib.Path, help="default: <final-dir>/comparisons.tsv")
    ap.add_argument("--mirnas", type=pathlib.Path, default=REPO / "metadata/v2/benchmark_mirnas.tsv")
    ap.add_argument("--predictions", type=pathlib.Path, default=REPO / "metadata/v2/predictions.tsv")
    ap.add_argument("--tools", nargs="*")
    ap.add_argument("--min-per-state", type=int, default=30)
    ap.add_argument("--bootstrap", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=0)
    return ap.parse_args(argv)


def main(argv=None) -> None:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
    evaluate(parse_args(argv))


if __name__ == "__main__":
    main()
