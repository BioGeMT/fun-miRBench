"""Relabel finished single-cell comparisons from their saved fits, into a new root: no store is read and nothing is refit.

Every dataset folder of --root with a comparisons.tsv and no SUPERSEDED file gets a folder under --out with comparisons.tsv
and genes/, and comparisons/<comparison_id>/ with the label logs. Stores and fits stay in --root. Only the s-value and
state columns of the d grid may change: any other difference makes the comparison an error. Finished tables under --out
are kept, so a second run adds what is new. The output has the same form as a --root, so
funmirbench-single-cell-combine works on it.
"""

from __future__ import annotations

import argparse
import json
import os
import pathlib
import shutil
import subprocess
import sys
import traceback
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone

import pandas as pd

from funmirbench.single_cell_final import LABEL_SCRIPT, STATES, combine_comparisons, last_lines, repo_root, sha256


def may_change(column: str) -> bool:
    """The FSOS columns of the d grid; state_raw_* come from the threshold tests and stay as they are."""
    return column.startswith("svalue_d") or column.startswith("state_d")


def check_only_labels_changed(old_path: pathlib.Path, new_path: pathlib.Path) -> None:
    old = pd.read_csv(old_path, sep="\t", dtype=str, keep_default_na=False)
    new = pd.read_csv(new_path, sep="\t", dtype=str, keep_default_na=False)
    if old.columns.tolist() != new.columns.tolist():
        raise ValueError(f"{new_path}: columns differ from {old_path}")
    changed = [column for column in old if not may_change(column) and not old[column].equals(new[column])]
    if changed:
        raise ValueError(f"{new_path}: columns other than the labels changed: {changed}")


def relabel_comparison(row: dict, old: pathlib.Path, new: pathlib.Path, gene_list: pathlib.Path,
                       repo: pathlib.Path) -> dict:
    """The comparison's row of the new comparisons.tsv. A comparison that was not fitted keeps its table and its row."""
    comparison_id = row["comparison_id"]
    old_table = old / "genes" / f"{comparison_id}.tsv.gz"
    new_table = new / "genes" / f"{comparison_id}.tsv.gz"
    if row["status"] != "fitted":
        if old_table.exists():
            shutil.copy2(old_table, new_table)
        return row
    log_directory = new / "comparisons" / comparison_id
    log_directory.mkdir(parents=True, exist_ok=True)
    if not new_table.exists():
        prefix = log_directory / "relabel"
        log = log_directory / "label_states.log"
        with open(log, "w") as handle:
            done = subprocess.run(["Rscript", str(repo / LABEL_SCRIPT),
                                   str(old / "comparisons" / comparison_id / "complete_analysis.rds"), str(gene_list),
                                   str(prefix), row["orientation"]],
                                  stdout=handle, stderr=subprocess.STDOUT, check=False)
        if done.returncode != 0:
            raise RuntimeError(f"label_states.R exit code {done.returncode}: {last_lines(log)}")
        for suffix in ("_events.tsv", "_metadata.txt"):
            pathlib.Path(f"{prefix}{suffix}").replace(log_directory / f"label_states{suffix}")
        check_only_labels_changed(old_table, pathlib.Path(f"{prefix}.tsv.gz"))
        pathlib.Path(f"{prefix}.tsv.gz").replace(new_table)
    states = pd.read_csv(new_table, sep="\t", usecols=["state_d05"]).state_d05.value_counts()
    return {**row, "table_sha256": sha256(new_table), **{f"d05_{state}": str(int(states.get(state, 0))) for state in STATES}}


def run_relabel(root: pathlib.Path, out: pathlib.Path, gene_list: pathlib.Path, *, workers: int,
                units_path: pathlib.Path | None = None, repo: pathlib.Path | None = None) -> pd.DataFrame:
    repo = (repo or repo_root()).resolve()
    tables = {}
    for path in sorted(root.glob("*/comparisons.tsv")):
        if not (path.parent / "SUPERSEDED").exists():
            tables[path.parent.name] = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    ids = pd.concat([table.comparison_id for table in tables.values()])
    if ids.duplicated().any():
        raise ValueError(f"{root}: duplicated comparison_id, for example {ids[ids.duplicated()].head(3).tolist()}")
    for dataset in tables:
        (out / dataset / "genes").mkdir(parents=True, exist_ok=True)

    def run_one(dataset: str, row: dict) -> dict:
        try:
            return relabel_comparison(row, root / dataset, out / dataset, gene_list, repo)
        except Exception as error:  # one comparison's failure must not stop the others; main exits non-zero
            traceback.print_exc()
            return {**row, "status": "error", "reason": f"{type(error).__name__}: {error}"[:400]}

    done = {}
    with ThreadPoolExecutor(workers) as pool:
        futures = {pool.submit(run_one, dataset, row): dataset for dataset, table in tables.items()
                   for row in table.to_dict("records")}
        for count, future in enumerate(as_completed(futures), start=1):
            row = future.result()
            done[row["comparison_id"]] = row
            print(f"PROGRESS {count}/{len(futures)} {row['comparison_id']} {row['status']}", flush=True)
    for dataset, table in tables.items():
        new_table = pd.DataFrame([done[comparison_id] for comparison_id in table.comparison_id])
        new_table.to_csv(out / dataset / "comparisons.tsv", sep="\t", index=False)
        manifest = {"built": datetime.now(timezone.utc).isoformat(timespec="seconds"), "fits_and_store": str((root / dataset).resolve()),
                    "gene_list_sha256": sha256(gene_list), "label_states_R_sha256": sha256(repo / LABEL_SCRIPT),
                    "status_counts": {key: int(value) for key, value in new_table.status.value_counts().items()}}
        (out / dataset / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    combined = combine_comparisons(out)
    if units_path:
        unit_status = combined.status.map({"fitted": "ok", "error": "failed"}).fillna("skipped")
        combined.assign(unit_status=unit_status)[["comparison_id", "unit_status", "reason"]].to_csv(
            units_path, sep="\t", index=False, header=False)
    return combined


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--root", type=pathlib.Path, required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    parser.add_argument("--genes", type=pathlib.Path, required=True)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--units", type=pathlib.Path)
    args = parser.parse_args()
    for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[variable] = "1"  # one R process per worker
    table = run_relabel(args.root, args.out, args.genes, workers=args.workers, units_path=args.units)
    print(table.status.value_counts().to_string())
    if (table.status == "error").any():
        sys.exit(1)


if __name__ == "__main__":
    main()
