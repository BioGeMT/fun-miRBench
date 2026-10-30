"""State counts and opposite-side genes across the margin grid of the FuNmiRBench v2 final tables (TSV outputs).

For every comparison and every margin d of the final tables (manifest.json "margins"):
  * the number of responsive, no-response and unresolved genes (state_<tag>);
  * target side: the responsive genes; opposite side: genes with a usable estimate, FSOS s-value < 0.05 against
    d and a shrunken change in the other direction, i.e. the responsive rule mirrored;
  * opposite / target ratio, the expected share of wrong responsive calls (mean over responsive genes of the
    probability that the change is smaller than d or of the other sign), and the median shrunken change and
    95% interval of the opposite-side genes;
  * the share of genes with an 8mer, 7mer-m8 or 7mer-A1 site of the comparison miRNA in the 3'UTR, on the target
    side, the opposite side and over all genes with a usable estimate.
margin_summary.tsv has one row per comparison at --main-margin plus, for each --ratio-limits value, the smallest
d where the ratio is at or below it with at least --min-target target-side genes (empty if none).

Usage:
  uv run python -m funmirbench.v2_margin_report --final-dir <final-v2-svalue> --utr-fasta <MANE 3'UTR fasta> \
      --out <dir>
"""
from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np
import pandas as pd
from scipy.stats import norm

from funmirbench.v2_evaluate import REPO, read_fasta_by_gene, seed_site_counts

COLUMNS = ["gene_id", "eligible", "logFC_shrunk", "lfcSE_shrunk", "oriented_logFC_shrunk",
           "oriented_ci95_low_shrunk", "oriented_ci95_high_shrunk"]


def has_site(utrs: dict[str, str], genes: pd.Series, mirna_sequence: str) -> np.ndarray:
    return np.array([g in utrs and sum(seed_site_counts(utrs[g], mirna_sequence)[:3]) > 0 for g in genes])


def margin_rows(table: pd.DataFrame, margins: dict[str, float], site: np.ndarray) -> list[dict]:
    est = table.eligible.astype(bool).to_numpy()
    o = table.oriented_logFC_shrunk.to_numpy(float)
    rows = []
    for tag, d in margins.items():
        state = table[f"state_{tag}"].to_numpy()
        target = state == "responsive"
        opposite = est & (table[f"svalue_{tag}"].to_numpy(float) < 0.05) & (o < 0)
        small = norm.cdf((d - np.abs(table.logFC_shrunk.to_numpy(float))) / table.lfcSE_shrunk.to_numpy(float))
        n_t, n_o = int(target.sum()), int(opposite.sum())
        rows.append({
            "margin": tag, "d": d,
            "n_responsive": n_t, "n_no_response": int((state == "no_response_above_margin").sum()),
            "n_unresolved": int((state == "unresolved").sum()), "n_opposite": n_o,
            "opposite_to_target_ratio": n_o / n_t if n_t else np.nan,
            "expected_wrong_share_responsive": float(small[target].mean()) if n_t else np.nan,
            "opposite_median_change": float(np.median(o[opposite])) if n_o else np.nan,
            "opposite_median_ci95_low": float(table.oriented_ci95_low_shrunk[opposite].median()) if n_o else np.nan,
            "opposite_median_ci95_high": float(table.oriented_ci95_high_shrunk[opposite].median()) if n_o else np.nan,
            "site_share_target": float(site[target].mean()) if n_t else np.nan,
            "site_share_opposite": float(site[opposite].mean()) if n_o else np.nan,
            "site_share_all": float(site[est].mean()),
        })
    return rows


def report(args: argparse.Namespace) -> None:
    margins = json.loads((args.final_dir / "manifest.json").read_text())["margins"]
    if args.main_margin not in margins:
        raise ValueError(f"--main-margin {args.main_margin} is not one of {list(margins)}")
    comps = pd.read_csv(args.final_dir / "comparisons.tsv", sep="\t", keep_default_na=False)
    comps = comps[(comps.action != "set_aside") & (comps.perturbation == "OE")]
    final = pd.read_csv(args.final_dir / "experiments.tsv", sep="\t").set_index("experiment_id")
    mirnas = pd.read_csv(args.mirnas, sep="\t").set_index("mirna_name")
    utrs = read_fasta_by_gene(args.utr_fasta)
    cols = COLUMNS + [f"{p}_{tag}" for tag in margins for p in ("state", "svalue")]
    rows = []
    for c in comps.itertuples(index=False):
        t = pd.read_csv(args.final_dir / final.loc[c.experiment_id, "table"], sep="\t", usecols=cols)
        site = has_site(utrs, t.gene_id, mirnas.loc[c.mirna_name, "mirna_sequence"])
        for r in margin_rows(t, margins, site):
            rows.append({"experiment_id": c.experiment_id, "action": c.action,
                         "dependency_group": c.dependency_group, **r})
    grid = pd.DataFrame(rows)
    args.out.mkdir(parents=True, exist_ok=True)
    grid.to_csv(args.out / "margin_grid.tsv", sep="\t", index=False)

    summary = grid[grid.margin == args.main_margin].set_index("experiment_id")
    for limit in args.ratio_limits:
        ok = grid[(grid.opposite_to_target_ratio <= limit) & (grid.n_responsive >= args.min_target)]
        summary[f"smallest_d_ratio_at_most_{limit:g}"] = ok.groupby("experiment_id").d.min().reindex(summary.index)
    summary.reset_index().to_csv(args.out / "margin_summary.tsv", sep="\t", index=False)


def parse_args(argv=None) -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--final-dir", type=pathlib.Path, required=True)
    ap.add_argument("--utr-fasta", type=pathlib.Path, required=True)
    ap.add_argument("--out", type=pathlib.Path, required=True)
    ap.add_argument("--mirnas", type=pathlib.Path, default=REPO / "metadata/v2/benchmark_mirnas.tsv")
    ap.add_argument("--main-margin", default="d05")
    ap.add_argument("--ratio-limits", type=float, nargs="+", default=[0.1, 0.5])
    ap.add_argument("--min-target", type=int, default=30)
    return ap.parse_args(argv)


def main(argv=None) -> None:
    report(parse_args(argv))


if __name__ == "__main__":
    main()
