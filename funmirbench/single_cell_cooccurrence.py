"""Co-occurrence datasets (small RNA and mRNA of the same cells): per group of cells (cell line or cell type), the
Spearman correlation across cells of each expressed miRNA with each expressed gene, its 95% confidence interval, the
expressed miRNA list and the per-cell tables.

Before ranking, genes are divided by the cell's mRNA table total, and miRNAs by the whole-library total of the library
that holds them (--library-totals; cells without a miRNA count stay as zeros). Without those totals, miRNAs are divided
by the cell's miRNA table total, and cells without a miRNA count are left out. So sequencing depth does not correlate
miRNAs with genes. The interval is the Bonett-Wright one: Fisher z with standard error sqrt((1 + r^2 / 2) / (n - 3)).
"""

from __future__ import annotations

import argparse
import json
import pathlib
import re

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import rankdata

from funmirbench.single_cell_store import check_counts

Z95 = 1.959964
MIN_CELLS_FOR_INTERVAL = 4


def read_cell_table(path: pathlib.Path, *, allow_fractional: bool = False) -> tuple[np.ndarray, np.ndarray, sparse.csr_matrix]:
    """Cells, features and counts (cells x features) from an h5ad (features: var.gene_id where it is set, otherwise
    var_names) or a TSV (first column the cell, one column per feature). Counts must be raw integers unless
    allow_fractional is set; they must never be negative or missing."""
    if path.name.endswith(".h5ad"):
        table = ad.read_h5ad(path)
        cells, features, counts = table.obs_names.to_numpy(dtype=object), table.var_names.to_numpy(dtype=object), table.X
        if "gene_id" in table.var:
            keep = table.var.gene_id.astype(str).to_numpy() != ""
            features, counts = table.var.gene_id.astype(str).to_numpy(dtype=object)[keep], counts[:, keep]
    else:
        table = pd.read_csv(path, sep="\t", index_col=0)
        table.index = table.index.astype(str)
        cells, features, counts = table.index.to_numpy(dtype=object), table.columns.to_numpy(dtype=object), table.to_numpy()
    counts = sparse.csr_matrix(counts)
    if allow_fractional:
        if counts.nnz and not (counts.data >= 0).all():
            raise ValueError(f"{path}: negative or missing values")
    else:
        check_counts(counts)
    for name, values in (("cell", cells), ("feature", features)):
        if len(pd.unique(values)) != len(values):
            raise ValueError(f"{path}: duplicated {name} names")
    return cells, features, counts


def spearman_table(mirna: np.ndarray, genes: np.ndarray, mirna_names: np.ndarray, gene_names: np.ndarray) -> pd.DataFrame:
    """Spearman correlation and 95% interval of every miRNA with every gene (cells x features arrays)."""
    def standardised_ranks(values: np.ndarray) -> np.ndarray:
        ranks = rankdata(values, axis=0)
        with np.errstate(invalid="ignore", divide="ignore"):
            return (ranks - ranks.mean(axis=0)) / ranks.std(axis=0, ddof=1)

    n = mirna.shape[0]
    r = np.clip(standardised_ranks(mirna).T @ standardised_ranks(genes) / (n - 1), -1, 1)
    half_width = Z95 * np.sqrt((1 + r ** 2 / 2) / (n - 3))
    with np.errstate(divide="ignore"):
        z = np.arctanh(r)
    return pd.DataFrame({"mirna": np.repeat(mirna_names, len(gene_names)),
                         "gene": np.tile(gene_names, len(mirna_names)),
                         "spearman_r": r.ravel(),
                         "ci95_low": np.tanh(z - half_width).ravel(),
                         "ci95_high": np.tanh(z + half_width).ravel(),
                         "n_cells": n})


def run_group(name: str, mirna: sparse.csr_matrix, genes: sparse.csr_matrix, mirna_names: np.ndarray,
              gene_names: np.ndarray, cells: np.ndarray, *, min_cells: int, min_detected: float,
              out: pathlib.Path, library_totals: np.ndarray | None = None) -> tuple[dict, pd.DataFrame]:
    mirna_total = np.asarray(mirna.sum(axis=1)).ravel()
    gene_library = np.asarray(genes.sum(axis=1)).ravel()
    mirna_library = mirna_total if library_totals is None else library_totals
    has_counts = (mirna_library > 0) & (gene_library > 0)
    mirna, genes, cells = mirna[has_counts], genes[has_counts], cells[has_counts]
    mirna_library, gene_library = mirna_library[has_counts], gene_library[has_counts]
    n_cells = len(cells)
    row = {"group": name, "n_cells": n_cells, "n_cells_without_counts": int((~has_counts).sum()),
           "n_cells_zero_mirna": int((mirna_total[has_counts] == 0).sum())}

    mirna_detected = np.asarray((mirna > 0).sum(axis=0)).ravel()
    gene_detected = np.asarray((genes > 0).sum(axis=0)).ravel()
    expressed_mirna = (mirna_detected >= min_detected * n_cells) & (mirna_detected > 0)
    expressed_gene = (gene_detected >= min_detected * n_cells) & (gene_detected > 0)
    mirna_counts = mirna[:, np.flatnonzero(expressed_mirna)].toarray()
    gene_counts = genes[:, np.flatnonzero(expressed_gene)].toarray()
    row.update(n_expressed_mirnas=int(expressed_mirna.sum()), n_expressed_genes=int(expressed_gene.sum()))

    expressed = pd.DataFrame({"group": name, "mirna": mirna_names[expressed_mirna],
                              "cells_detected": mirna_detected[expressed_mirna],
                              "fraction_detected": mirna_detected[expressed_mirna] / max(n_cells, 1),
                              "total_counts": mirna_counts.sum(axis=0)})
    for label, values, names, library in (("mirna", mirna_counts, mirna_names[expressed_mirna], mirna_library),
                                          ("mrna", gene_counts, gene_names[expressed_gene], gene_library)):
        per_cell = pd.DataFrame(values, columns=names, index=pd.Index(cells, name="cell"))
        per_cell.insert(0, "library_size", library)
        per_cell.to_csv(out / "per_cell" / f"{name}_{label}.tsv.gz", sep="\t")

    if n_cells < min_cells:
        return {**row, "status": "not_correlated", "reason": f"{n_cells} cells with counts, fewer than {min_cells}"}, expressed
    if not len(mirna_counts.T) or not len(gene_counts.T):
        return {**row, "status": "not_correlated", "reason": "no expressed miRNA or no expressed gene"}, expressed
    table = spearman_table(mirna_counts / mirna_library[:, None], gene_counts / gene_library[:, None],
                           mirna_names[expressed_mirna], gene_names[expressed_gene])
    table.to_csv(out / "correlations" / f"{name}.tsv.gz", sep="\t", index=False, float_format="%.5g")
    return {**row, "status": "correlated", "reason": ""}, expressed


def read_library_totals(path: pathlib.Path, cells: pd.Index, mirna: sparse.csr_matrix) -> np.ndarray:
    """Whole-library total of each cell, in the order of cells; mirna holds the same cells in the same order."""
    table = pd.read_csv(path, sep="\t", dtype={"cell": str})
    if not {"cell", "library_total"} <= set(table):
        raise ValueError(f"{path}: needs columns cell and library_total")
    if table.cell.duplicated().any():
        raise ValueError(f"{path}: duplicated cells")
    totals = table.set_index("cell").library_total.reindex(cells).to_numpy(dtype=float)
    if not (totals >= 0).all():
        raise ValueError(f"{path}: a cell has no library total, or a negative one")
    if (np.asarray(mirna.sum(axis=1)).ravel() > totals * (1 + 1e-4) + 1e-4).any():
        raise ValueError(f"{path}: a cell's miRNA counts exceed its library total")
    return totals


def run_cooccurrence(mirna_path: pathlib.Path, mrna_path: pathlib.Path, cells_path: pathlib.Path, group_column: str,
                     out: pathlib.Path, *, min_cells: int = 10, min_detected: float = 0.1, allow_fractional: bool = False,
                     library_totals_path: pathlib.Path | None = None) -> pd.DataFrame:
    if min_cells < MIN_CELLS_FOR_INTERVAL:
        raise ValueError(f"--min-cells must be at least {MIN_CELLS_FOR_INTERVAL} for a confidence interval")
    mirna_cells, mirna_names, mirna = read_cell_table(mirna_path, allow_fractional=allow_fractional)
    gene_cells, gene_names, genes = read_cell_table(mrna_path, allow_fractional=allow_fractional)
    groups = pd.read_csv(cells_path, sep="\t", dtype=str)
    if not {"cell", group_column} <= set(groups):
        raise ValueError(f"{cells_path}: needs columns cell and {group_column}")
    if groups.cell.duplicated().any():
        raise ValueError(f"{cells_path}: duplicated cells")
    labelled = groups.dropna(subset=[group_column])
    labelled = labelled[labelled[group_column] != ""].set_index("cell")[group_column]
    shared = pd.Index(mirna_cells).intersection(pd.Index(gene_cells)).intersection(labelled.index)
    if shared.empty:
        raise ValueError("no cell is in the miRNA table, the mRNA table and the group table")
    mirna_position = pd.Series(np.arange(len(mirna_cells)), index=mirna_cells)
    gene_position = pd.Series(np.arange(len(gene_cells)), index=gene_cells)
    labels = labelled.reindex(shared).to_numpy()
    mirna_normalisation = "miRNA_table_total" if library_totals_path is None else "whole_library_total"
    totals = None
    if library_totals_path is not None:
        totals = read_library_totals(library_totals_path, shared, mirna[mirna_position[shared].to_numpy()])
    file_names = {group: re.sub(r"[^A-Za-z0-9._-]+", "_", group) for group in pd.unique(labels)}
    if len(set(file_names.values())) != len(file_names):
        raise ValueError(f"groups collide as file names: {file_names}")

    (out / "correlations").mkdir(parents=True, exist_ok=False)
    (out / "per_cell").mkdir()
    rows, expressed = [], []
    for group, file_name in sorted(file_names.items()):
        members = shared[labels == group]
        row, listed = run_group(file_name, mirna[mirna_position[members].to_numpy()],
                                genes[gene_position[members].to_numpy()], mirna_names, gene_names,
                                members.to_numpy(dtype=object), min_cells=min_cells, min_detected=min_detected, out=out,
                                library_totals=None if totals is None else totals[labels == group])
        rows.append({**row, "group": group, "mirna_normalisation": mirna_normalisation})
        expressed.append(listed.assign(group=group))
        print(f"PROGRESS {len(rows)}/{len(file_names)} {group} {row['status']}", flush=True)
    table = pd.DataFrame(rows)
    table.to_csv(out / "groups.tsv", sep="\t", index=False)
    pd.concat(expressed, ignore_index=True).to_csv(out / "expressed_mirnas.tsv", sep="\t", index=False)
    (out / "manifest.json").write_text(json.dumps(
        {"mirna": str(mirna_path.resolve()), "mrna": str(mrna_path.resolve()), "cells": str(cells_path.resolve()),
         "group_column": group_column, "min_cells": min_cells, "min_detected_fraction": min_detected,
         "allow_fractional": allow_fractional, "mirna_normalisation": mirna_normalisation,
         "library_totals": None if library_totals_path is None else str(library_totals_path.resolve()),
         "cells_in_miRNA_table": len(mirna_cells), "cells_in_mRNA_table": len(gene_cells),
         "cells_in_group_table": len(groups), "cells_used": len(shared)}, indent=2) + "\n")
    return table


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--mirna", type=pathlib.Path, required=True)
    parser.add_argument("--mrna", type=pathlib.Path, required=True)
    parser.add_argument("--cells", type=pathlib.Path, required=True)
    parser.add_argument("--group-column", required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    parser.add_argument("--min-cells", type=int, default=10)
    parser.add_argument("--min-detected-fraction", type=float, default=0.1)
    parser.add_argument("--allow-fractional", action="store_true")
    parser.add_argument("--library-totals", type=pathlib.Path)
    args = parser.parse_args()
    table = run_cooccurrence(args.mirna, args.mrna, args.cells, args.group_column, args.out,
                             min_cells=args.min_cells, min_detected=args.min_detected_fraction,
                             allow_fractional=args.allow_fractional, library_totals_path=args.library_totals)
    print(table.to_string(index=False))


if __name__ == "__main__":
    main()
