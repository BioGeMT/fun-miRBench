"""Build the final v2 benchmark tables from saved DESeq2 analyses.

For each comparison, label_states.R writes one row per gene of the benchmark gene list: statistics, shrunken
log2FC, FSOS s-values and the three states for each d of the grid (from the shrunken log2FC), and the raw threshold
tests and raw states at d = 0.5 and 1.0 as a check. This driver runs it for every comparison and writes a summary
table and a manifest with SHA-256 hashes.

Input --fits is a TSV with columns: experiment_id, perturbation (OE or KO), fit_rds (complete_analysis.rds written
by pipelines/experiments/run_deseq2_counts.R).

The comparison metadata (metadata/v2/comparisons.tsv: groups, flags, KO evidence) is copied into the output.

Usage:
  uv run python pipelines/final_data/build_final_data.py --fits fits.tsv \
      --genes metadata/v2/benchmark_genes.tsv.gz --out <final_dir> [--rscript Rscript] [--workers 8]
"""
import argparse
import hashlib
import json
import pathlib
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone

import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
MARGINS = ("d05", "d10")
D_GRID = {"d00": 0.0, "d01": 0.1, "d025": 0.25, "d05": 0.5, "d075": 0.75, "d10": 1.0, "d15": 1.5, "d20": 2.0,
          "d30": 3.0}
STATES = ("responsive", "no_response_above_margin", "unresolved")


def sha256(path: pathlib.Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def run_one(row, *, genes: pathlib.Path, out: pathlib.Path, rscript: str) -> dict:
    prefix = out / "genes" / row.experiment_id
    log = out / "logs" / f"{row.experiment_id}.log"
    with open(log, "w") as fh:
        subprocess.run([rscript, str(HERE / "label_states.R"), str(row.fit_rds), str(genes), row.perturbation,
                        str(prefix)], stdout=fh, stderr=subprocess.STDOUT, check=True)
    for suffix in ("_events.tsv", "_metadata.txt"):
        pathlib.Path(f"{prefix}{suffix}").replace(out / "logs" / f"{row.experiment_id}{suffix}")
    table = pathlib.Path(f"{prefix}.tsv.gz")
    tags = [*MARGINS, *[f"raw_{m}" for m in MARGINS]]
    t = pd.read_csv(table, sep="\t", usecols=["eligible", "analysis_status", *[f"state_{m}" for m in tags]])
    rec = {"experiment_id": row.experiment_id, "perturbation": row.perturbation, "n_genes": len(t),
           "n_fitted": int((t.analysis_status != "not_fitted_all_zero").sum()), "n_eligible": int(t.eligible.sum())}
    for m in tags:
        counts = t[f"state_{m}"].value_counts()
        for s in STATES:
            rec[f"{m}_{s}"] = int(counts.get(s, 0))
    rec.update({"fit_rds": str(row.fit_rds), "fit_sha256": sha256(pathlib.Path(row.fit_rds)),
                "table": f"genes/{row.experiment_id}.tsv.gz", "table_sha256": sha256(table)})
    return rec


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--fits", type=pathlib.Path, required=True)
    ap.add_argument("--genes", type=pathlib.Path, required=True)
    ap.add_argument("--out", type=pathlib.Path, required=True)
    ap.add_argument("--rscript", default="Rscript")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--comparisons", type=pathlib.Path, default=pathlib.Path("metadata/v2/comparisons.tsv"),
                    help="copied into the output so the evaluation reads one folder")
    args = ap.parse_args()
    fits = pd.read_csv(args.fits, sep="\t")
    if fits.experiment_id.duplicated().any() or not fits.perturbation.isin(["OE", "KO"]).all():
        raise ValueError("--fits needs unique experiment_id values and perturbation OE or KO")
    missing = [p for p in fits.fit_rds if not pathlib.Path(p).exists()]
    if missing:
        raise FileNotFoundError(f"missing fits: {missing}")
    (args.out / "genes").mkdir(parents=True, exist_ok=False)
    (args.out / "logs").mkdir(parents=True, exist_ok=False)
    with ThreadPoolExecutor(args.workers) as pool:
        rows = list(pool.map(lambda r: run_one(r, genes=args.genes, out=args.out, rscript=args.rscript),
                             fits.itertuples(index=False)))
    summary = pd.DataFrame(rows)
    summary.to_csv(args.out / "experiments.tsv", sep="\t", index=False)
    comps = pd.read_csv(args.comparisons, sep="\t")
    if set(comps.experiment_id) != set(summary.experiment_id):
        raise ValueError("--comparisons must list exactly the fitted comparisons")
    shutil.copyfile(args.comparisons, args.out / "comparisons.tsv")
    manifest = {
        "built": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "gene_list": str(args.genes), "gene_list_sha256": sha256(args.genes),
        "label_states_R_sha256": sha256(HERE / "label_states.R"),
        "states": {"responsive": "oriented shrunken log2FC > 0 and FSOS s-value against d < 0.05 (over listed genes)",
                   "no_response_above_margin": "oriented shrunken log2FC + 1.959964 * lfcSE_shrunk < d",
                   "unresolved": "everything else, including genes with no usable estimate"},
        "raw_states": {"responsive": "threshold-test FDR < 0.05 (BH over listed genes after independent filtering)",
                       "no_response_above_margin": "oriented log2FC + 1.959964 * lfcSE < d"},
        "margins": D_GRID,
        "raw_margins": {"d05": 0.5, "d10": 1.0},
        "experiments": len(summary),
        "comparisons_sha256": sha256(args.out / "comparisons.tsv"),
        "tables": {r["experiment_id"]: r["table_sha256"] for r in rows},
    }
    (args.out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(summary[["experiment_id", "perturbation", "d05_responsive", "d05_no_response_above_margin",
                   "d05_unresolved"]].to_string(index=False))


if __name__ == "__main__":
    main()
