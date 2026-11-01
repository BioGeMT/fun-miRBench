"""Check one tier 2 dataset and write its row in final_build/outcomes/tier2.tsv (one row per dataset, replaced on rerun).
Usage: tier2_record.py <dataset> <category> <source_files> <status> <reason> [--final-dir DIR] [--store H5AD] [--solo PREFIX=DIR ...]
--solo checks, for the store cells named PREFIX:<barcode> (plain barcodes when PREFIX is empty), that each cell's total equals the alignment's matrix."""
import argparse
import pathlib
import subprocess

import numpy as np
import pandas as pd

from funmirbench.single_cell_store import read_10x_mtx, read_cell_store

OUTCOMES = pathlib.Path("/SCRATCH/dtzim01/fun-miRBind-work/single_cell_mirna/final_build/outcomes/tier2.tsv")
COLUMNS = ["dataset", "category", "source_files", "size_gb", "cells_in_store", "comparisons", "fitted", "not_fitted", "failed", "final_dir", "status", "reason"]
FEATURES_LEFT_OUT_TOLERANCE = 0.02  # features repeated within a library are left out of the store


def check_against_alignment(store, prefix: str, folder: pathlib.Path) -> str:
    counts, _, barcodes = read_10x_mtx(folder / "matrix.mtx", folder / "features.tsv", folder / "barcodes.tsv")
    names = pd.Index([f"{prefix}:{barcode}" if prefix else barcode for barcode in barcodes])
    cells = store.obs_names[store.obs_names.str.startswith(f"{prefix}:")] if prefix else store.obs_names
    position = names.get_indexer(cells)
    if len(cells) == 0 or (position < 0).any():
        raise SystemExit(f"{prefix}: store cells missing from {folder}")
    aligned = np.asarray(counts.tocsc()[:, position].sum(axis=0)).ravel()
    stored = np.asarray(store[cells].X.sum(axis=1)).ravel()
    gap = np.abs(aligned - stored) / np.maximum(aligned, 1)
    if gap.max() > FEATURES_LEFT_OUT_TOLERANCE:
        raise SystemExit(f"{prefix}: store totals differ from the alignment by up to {gap.max():.3f}")
    return f"{prefix}: {len(cells)} cells, largest total difference to the alignment {gap.max():.4f}"


def main() -> None:
    parser = argparse.ArgumentParser()
    for name in ("dataset", "category", "source_files", "status", "reason"):
        parser.add_argument(name)
    parser.add_argument("--final-dir", type=pathlib.Path)
    parser.add_argument("--store", type=pathlib.Path)
    parser.add_argument("--solo", action="append", default=[])
    args = parser.parse_args()
    row = dict.fromkeys(COLUMNS, 0)
    row.update(dataset=args.dataset, category=args.category, source_files=args.source_files, status=args.status, reason=args.reason, final_dir="")
    paths = [path for path in args.source_files.split(";") if path and pathlib.Path(path).exists()]
    if paths:
        row["size_gb"] = round(sum(int(subprocess.run(["du", "-sb", path], capture_output=True, text=True, check=True).stdout.split()[0]) for path in paths) / 1e9, 2)
    if args.store:
        store = read_cell_store(args.store)
        row["cells_in_store"] = store.n_obs
        print(f"{args.dataset}: {store.n_obs} cells, {int(store.obs.passed_qc.sum())} pass QC, {store.n_vars} features, {store.X.sum()} counts; store valid", flush=True)
        for pair in args.solo:
            prefix, folder = pair.split("=", 1)
            print(check_against_alignment(store, prefix, pathlib.Path(folder)), flush=True)
    if args.final_dir:
        row["final_dir"] = str(args.final_dir)
        listing = args.final_dir / "comparisons.tsv"
        groups = args.final_dir / "correlation/groups.tsv"
        if groups.exists():  # co-occurrence: groups, not comparisons; fitted means correlated
            table = pd.read_csv(groups, sep="\t", keep_default_na=False)
            row.update(cells_in_store=int(table.n_cells.sum()), comparisons=len(table), fitted=int((table.status == "correlated").sum()), not_fitted=int((table.status != "correlated").sum()))
        elif listing.exists():
            table = pd.read_csv(listing, sep="\t", keep_default_na=False)
            row.update(comparisons=len(table), fitted=int((table.status == "fitted").sum()),
                       not_fitted=int((table.status == "not_fitted").sum()), failed=int(table.status.isin(["fit_failed", "error"]).sum()))
            missing = [c for c in table.table if c and not (args.final_dir / c).exists()]
            if missing:
                raise SystemExit(f"tables listed but missing: {missing[:3]}")
            print(table[["comparison_id", "status", "reason"]].to_string(index=False), flush=True)
    table = pd.read_csv(OUTCOMES, sep="\t", keep_default_na=False) if OUTCOMES.exists() else pd.DataFrame(columns=COLUMNS)
    table = pd.concat([table[table.dataset != args.dataset], pd.DataFrame([row])], ignore_index=True)
    temporary = OUTCOMES.with_suffix(".tsv.tmp")
    table[COLUMNS].to_csv(temporary, sep="\t", index=False)
    temporary.replace(OUTCOMES)


if __name__ == "__main__":
    main()
