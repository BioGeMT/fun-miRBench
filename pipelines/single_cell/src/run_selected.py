"""Build the cell store of each selected dataset (builds.tsv) and fit only the selected comparisons
(selected_comparisons.tsv), then write the combined table. A finished dataset is not run again.
"""
from __future__ import annotations

import argparse
import os
import pathlib
import shlex
import subprocess
import sys

import pandas as pd

from funmirbench.single_cell_final import combine_comparisons, run_final

SINGLE_CELL = pathlib.Path(__file__).resolve().parents[1]
ADAPTERS = SINGLE_CELL / "adapters"


def adapter_environment(data: pathlib.Path, final: pathlib.Path, gtf: pathlib.Path) -> dict[str, str]:
    return {**os.environ, "SINGLE_CELL_DATA": str(data), "SINGLE_CELL_FINAL": str(final), "SINGLE_CELL_GTF": str(gtf)}


def run_dataset(row, selected: pd.DataFrame, *, data: pathlib.Path, final: pathlib.Path, gtf: pathlib.Path,
                gene_list: pathlib.Path, workers: int) -> str:
    """Empty if it worked, otherwise why not."""
    out = final / row.dataset
    if (out / "comparisons.tsv").exists():
        print(f"done    {row.dataset}", flush=True)
        return ""
    store = out / "store" / f"{row.dataset}.h5ad"
    if not store.exists():
        store.parent.mkdir(parents=True, exist_ok=True)
        command = row.adapter_command.format(store=store, data=data)
        print(f"store   {row.dataset}: {command}", flush=True)
        done = subprocess.run([sys.executable, *shlex.split(command)], cwd=ADAPTERS, check=False,
                              env=adapter_environment(data, final, gtf))
        if done.returncode != 0 or not store.exists():
            return f"adapter exit code {done.returncode}"
    wanted = set(selected.comparison_id[selected.dataset == row.dataset])
    listed = pd.read_csv(SINGLE_CELL / "comparisons" / f"{row.dataset}.tsv", sep="\t", dtype=str, keep_default_na=False)
    chosen = listed[listed.comparison_id.isin(wanted)]
    if set(chosen.comparison_id) != wanted:
        return f"not in the comparison list: {sorted(wanted - set(chosen.comparison_id))}"
    work = final / "work"
    work.mkdir(parents=True, exist_ok=True)
    list_path = work / f"{row.dataset}.tsv"
    chosen.to_csv(list_path, sep="\t", index=False)
    table = run_final(store, list_path, gene_list, out, workers=workers)
    errors = table[table.status == "error"]
    return f"{len(errors)} comparisons with errors" if len(errors) else ""


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--data", type=pathlib.Path, required=True)
    parser.add_argument("--final", type=pathlib.Path, required=True)
    parser.add_argument("--gtf", type=pathlib.Path, required=True)
    parser.add_argument("--genes", type=pathlib.Path, required=True)
    parser.add_argument("--groups", default="main,machinery,supplement")
    parser.add_argument("--dataset", action="append")
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--units", type=pathlib.Path)
    args = parser.parse_args()
    for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[variable] = "1"  # one R process per worker

    builds = pd.read_csv(SINGLE_CELL / "builds.tsv", sep="\t", dtype=str, keep_default_na=False)
    selected = pd.read_csv(SINGLE_CELL / "selected_comparisons.tsv", sep="\t", dtype=str, keep_default_na=False)
    groups = set(args.groups.split(","))
    builds = builds[builds.group.isin(groups)]
    if args.dataset:
        builds = builds[builds.dataset.isin(args.dataset)]
    if builds.empty:
        raise ValueError(f"no builds for groups {sorted(groups)} and datasets {args.dataset}")
    units, failed = [], 0
    for count, row in enumerate(builds.itertuples(index=False), start=1):
        reason = run_dataset(row, selected, data=args.data, final=args.final, gtf=args.gtf, gene_list=args.genes,
                             workers=args.workers)
        units.append((row.dataset, "failed" if reason else "ok", reason))
        failed += bool(reason)
        if reason:
            print(f"FAILED  {row.dataset}: {reason}", flush=True)
        print(f"PROGRESS {count}/{len(builds)} {row.dataset}", flush=True)
    if args.units:
        pd.DataFrame(units).to_csv(args.units, sep="\t", index=False, header=False)
    if failed:
        sys.exit(1)
    combined = combine_comparisons(args.final)
    merged = combined.merge(selected[["comparison_id", "group"]], on="comparison_id", how="left")
    merged.to_csv(args.final / "comparisons.tsv", sep="\t", index=False)
    print(merged.groupby(["group", "status"]).size().to_string())


if __name__ == "__main__":
    main()
