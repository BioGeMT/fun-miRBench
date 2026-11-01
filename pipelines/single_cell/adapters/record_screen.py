"""Adds one dataset's row to the outcome table from its comparisons.tsv. Usage: record_screen.py <dataset> <source_files>
<size_gb> (or: record_screen.py <dataset> <source_files> <size_gb> blocked|failed <reason>)."""

import sys

import anndata as ad
import pandas as pd

from screen_common import FINAL, record_outcome


def main() -> None:
    dataset, source_files, size_gb = sys.argv[1:4]
    row = {"dataset": dataset, "source_files": source_files, "size_gb": size_gb}
    if len(sys.argv) > 4:
        row.update(status=sys.argv[4], reason=sys.argv[5])
        record_outcome(row)
        return
    table = pd.read_csv(FINAL / dataset / "comparisons.tsv", sep="\t", dtype=str, keep_default_na=False)
    base = dataset.removesuffix("_pc_hosts")
    store = ad.read_h5ad(FINAL / base / "store" / f"{base}.h5ad", backed="r")
    counts = table.status.value_counts()
    row.update(category=";".join(sorted(table.category.unique())), cells_in_store=store.n_obs,
               comparisons=len(table), fitted=int(counts.get("fitted", 0)), not_fitted=int(counts.get("not_fitted", 0)),
               failed=int(counts.get("fit_failed", 0) + counts.get("error", 0)), final_dir=str(FINAL / dataset),
               status="done")
    record_outcome(row)
    print(row)


if __name__ == "__main__":
    main()
