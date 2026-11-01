"""Shared steps of the pooled-screen adapters: which arms are targets, reading cells from deposits, writing the cell
store and the comparison list, and the outcome row. No expression or response value is used to choose anything here."""

from __future__ import annotations

import fcntl
import os
import pathlib
import re

import h5py
import numpy as np
import pandas as pd
from scipy import sparse, stats

from funmirbench.single_cell_store import write_cell_store

ADAPTERS = pathlib.Path(__file__).resolve().parent
REPO = ADAPTERS.parents[2]
TARGET_TABLE = ADAPTERS / "target_genes.tsv"
DATA = pathlib.Path(os.environ.get("SINGLE_CELL_DATA", REPO / "data/single_cell/raw"))
FINAL = pathlib.Path(os.environ.get("SINGLE_CELL_FINAL", REPO / "data/single_cell/final"))
GTF = pathlib.Path(os.environ.get("SINGLE_CELL_GTF", DATA / "refs/Homo_sapiens.GRCh38.115.gtf.gz"))
OUTCOMES = pathlib.Path("/SCRATCH/dtzim01/fun-miRBind-work/single_cell_mirna/final_build/outcomes/screens.tsv")
OUTCOME_COLUMNS = ("dataset", "category", "source_files", "size_gb", "cells_in_store", "comparisons", "fitted",
                   "not_fitted", "failed", "final_dir", "status", "reason")
MIN_CELLS_PER_UNIT = 3
PLAIN_LANE_SET_MIN_LANES = 16


class Targets:
    """target_genes.tsv: symbol, Ensembl IDs, category, host_type, miRNAs inside, and the aliases safe to match (not an
    approved symbol of another gene, unique in the list). RAX is not in it: it is its own gene, not PRKRA."""

    def __init__(self):
        table = pd.read_csv(TARGET_TABLE, sep="\t", dtype=str, keep_default_na=False)
        self.table = table.set_index("symbol", drop=False)
        self.names = {}
        for row in table.itertuples(index=False):
            for name in (*row.aliases.split(";"), *row.gene_ids.split(";"), row.symbol):
                if name:
                    self.names[name] = row.symbol

    def symbol_of(self, label: str) -> str | None:
        label = str(label).split(".")[0] if str(label).startswith("ENSG") else str(label)
        return self.names.get(label)


def lane_sets(lanes) -> tuple[np.ndarray, str]:
    """Lane set per cell: eight consecutive lanes when the screen has enough lanes for at least two sets of eight,
    otherwise every lane is its own unit. Returns the units and the rule used, for the caveats."""
    lanes = np.asarray(lanes)
    observed, position = np.unique(lanes, return_inverse=True)
    if len(observed) >= PLAIN_LANE_SET_MIN_LANES:
        n_sets = len(observed) // 8
        sets = np.minimum(position // 8, n_sets - 1)
        return np.array([f"lane_set_{value + 1}" for value in sets], dtype=object), \
            f"replicate units are sets of 8 consecutive lanes ({n_sets} sets of {len(observed)} lanes)"
    return np.array([f"lane_{value}" for value in lanes], dtype=object), \
        f"replicate units are single lanes ({len(observed)} lanes, too few for sets of 8)"


def natural_order(values) -> np.ndarray:
    """Integer rank of each text lane name in natural order (lane_2 before lane_10), for lane_sets."""
    names = sorted(set(values), key=lambda name: [int(part) if part.isdigit() else part
                                                  for part in re.split(r"(\d+)", str(name))])
    rank = {name: position for position, name in enumerate(names)}
    return np.array([rank[value] for value in values], dtype=np.int64)


def decode(values) -> np.ndarray:
    return np.array([value.decode() if isinstance(value, bytes) else value for value in values], dtype=object)


def h5_frame(group: h5py.Group) -> pd.DataFrame:
    """An h5ad obs or var group, old layout (__categories) or new (categorical subgroups). Numbers stay numbers."""
    index_name = group.attrs.get("_index", "_index")
    categories = group["__categories"] if "__categories" in group else {}
    out = {}
    for key in group.keys():
        if key == "__categories":
            continue
        item = group[key]
        if isinstance(item, h5py.Group):
            if "codes" not in item:
                continue
            codes, labels = item["codes"][()], decode(item["categories"][()])
            out[key] = np.where(codes >= 0, labels[np.maximum(codes, 0)], "") if len(labels) else np.full(len(codes), "", dtype=object)
            continue
        values = item[()]
        if key in categories:
            labels = decode(categories[key][()])
            values = np.where(values >= 0, labels[np.maximum(values, 0)], "")
        elif values.dtype.kind in "OS":
            values = decode(values)
        out[key] = values
    return pd.DataFrame(out).set_index(index_name, drop=False)


def h5_rows(handle: h5py.File, rows: np.ndarray, key: str = "X", block: int = 16384,
            integer: bool = True) -> sparse.csr_matrix:
    """Rows of a matrix in an h5ad (dense or CSR), by blocks of sorted rows so only the needed part of a large file is
    touched. integer: CSR int32 after checking the values are integers; otherwise CSR float64 as stored."""
    rows = np.asarray(rows)
    if len(rows) and not (np.diff(rows) > 0).all():
        raise ValueError("rows must be sorted and unique")
    node = handle[key]
    parts = []
    if isinstance(node, h5py.Dataset):
        for start in range(0, len(rows), block):
            chunk = rows[start:start + block]
            first, last = int(chunk.min()), int(chunk.max()) + 1
            parts.append(sparse.csr_matrix(node[first:last][chunk - first], dtype=np.float64))
    else:
        encoding = node.attrs.get("encoding-type", "csr_matrix")
        encoding = encoding.decode() if isinstance(encoding, bytes) else str(encoding)
        if encoding != "csr_matrix":
            raise ValueError(f"{key}: expected CSR or dense, found {encoding}")
        indptr = node["indptr"][()]
        n_columns = int(node.attrs.get("shape", node.attrs.get("h5sparse_shape"))[1])
        for start in range(0, len(rows), block):
            chunk = rows[start:start + block]
            first, last = int(chunk.min()), int(chunk.max()) + 1
            span_start, span_end = int(indptr[first]), int(indptr[last])
            whole = sparse.csr_matrix((np.asarray(node["data"][span_start:span_end], np.float64),
                                       node["indices"][span_start:span_end], indptr[first:last + 1] - span_start),
                                      shape=(last - first, n_columns))
            parts.append(whole[chunk - first])
    counts = sparse.vstack(parts).tocsr()
    if not integer:
        return counts
    if counts.data.size and not np.array_equal(counts.data, np.round(counts.data)):
        raise ValueError(f"{key} holds non-integer values; it is not a raw count matrix")
    return counts.astype(np.int32)


def plan_cells(obs: pd.DataFrame, targets: Targets, *, all_cells: bool = False) -> tuple[pd.DataFrame, np.ndarray]:
    """obs of all passing and failing cells of a deposit, with the store's columns and arm as deposited. Renames target
    arms to the listed symbol. Returns the rows to keep and their positions: target cells, plus control cells of the
    batches where those cells occur (every cell at high guide load, where controls are the other cells of the batch)."""
    obs = obs.reset_index(drop=True).copy()
    if all_cells:
        arms = obs.arm.astype(str).str.split(";")
        obs["arm"] = [";".join(targets.symbol_of(part) or part for part in value) for value in arms]
        return obs, np.arange(len(obs))
    renamed = obs.arm.astype(str).map(targets.symbol_of)
    is_target = renamed.notna().to_numpy() & ~obs.is_control.to_numpy()
    obs.loc[is_target, "arm"] = renamed[is_target]
    batches = set(obs.batch[is_target])
    keep = is_target | (obs.is_control.to_numpy() & obs.batch.isin(batches).to_numpy())
    return obs[keep].reset_index(drop=True), np.flatnonzero(keep)


def comparison_rows(dataset: str, arms, targets: Targets, *, control_arm: str, caveats: str, raised: bool = False,
                    other_pathway_arms: str = "") -> pd.DataFrame:
    """One row per target arm. raised: the guide or ORF raises the target (CRISPRa, ORF): orientation OE."""
    rows = []
    for symbol in sorted(set(arms)):
        target = targets.table.loc[symbol]
        notes = [caveats] if caveats else []
        if target.category == "host_gene":
            action = "activated" if raised else "repressed"
            notes.append(f"host gene or locus {action}, miRNAs inside: {target.mirnas_inside}; the mature miRNA was not "
                         f"measured, and the perturbation also changes the host transcript")
            if target.host_type == "protein_coding_host":
                notes.append(f"{symbol} is protein-coding: the perturbation also changes the {symbol} protein")
        if raised and symbol in ("LIN28A", "LIN28B"):
            notes.append(f"{symbol} raised, orientation OE records the gain of {symbol}; LIN28 represses let-7 "
                         f"maturation, so let-7 activity is expected to go down")
        elif raised and symbol == "TARBP2":
            notes.append("TARBP2 raised, orientation OE records the gain of TARBP2, a DICER1 cofactor")
        elif symbol in ("LIN28A", "LIN28B"):
            notes.append(f"{symbol} lowered, orientation KO records the loss of {symbol}; LIN28 represses let-7 "
                         f"maturation, so let-7 activity is expected to go up")
        rows.append({"comparison_id": f"{dataset}__{symbol}_{'OE' if raised else 'KO'}", "dataset": dataset,
                     "category": target.category, "mirna_or_gene": symbol, "orientation": "OE" if raised else "KO",
                     "treated_arm": symbol, "control_arm": control_arm, "min_cells": MIN_CELLS_PER_UNIT,
                     "caveats": "; ".join(notes), "target_gene_id": target.gene_ids,
                     "other_pathway_arms": other_pathway_arms})
    return pd.DataFrame(rows)


def write_screen(dataset: str, counts, obs: pd.DataFrame, features, *, source: str, comparisons: pd.DataFrame,
                 symbols=None, gene_ids=None, folder: str | None = None, fractional: bool = False) -> pathlib.Path:
    """The store under FINAL/<folder, default dataset>/store and the lists under the worktree's comparisons folder: <dataset>.tsv for
    pathway and non-coding host targets; <dataset>_pc_hosts.tsv for protein-coding hosts, run last and reported apart."""
    obs = obs.assign(dataset=dataset)
    store_path = FINAL / (folder or dataset) / "store" / f"{dataset}.h5ad"
    store_path.parent.mkdir(parents=True, exist_ok=True)
    write_cell_store(store_path, counts, obs, features, gtf=GTF, source=source, symbols=symbols, gene_ids=gene_ids,
                     fractional=fractional)
    table = Targets().table
    coding = comparisons.mirna_or_gene.isin(table.index[table.host_type == "protein_coding_host"])
    comparisons[~coding].to_csv(ADAPTERS.parent / "comparisons" / f"{dataset}.tsv", sep="\t", index=False)
    if coding.any():
        comparisons[coding].to_csv(ADAPTERS.parent / "comparisons" / f"{dataset}_pc_hosts.tsv", sep="\t", index=False)
    print(f"{dataset}: {counts.shape[0]} cells, {comparisons.shape[0]} comparisons, store {store_path}", flush=True)
    return store_path


def record_outcome(row: dict) -> None:
    """Add or replace the dataset's row of the outcome table; safe with several runs writing at once."""
    OUTCOMES.parent.mkdir(parents=True, exist_ok=True)
    with open(OUTCOMES.with_suffix(".lock"), "w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        table = pd.read_csv(OUTCOMES, sep="\t", dtype=str, keep_default_na=False) if OUTCOMES.exists() \
            else pd.DataFrame(columns=OUTCOME_COLUMNS)
        table = table[table.dataset != row["dataset"]]
        table = pd.concat([table, pd.DataFrame([{key: row.get(key, "") for key in OUTCOME_COLUMNS}])], ignore_index=True)
        table.sort_values("dataset").to_csv(OUTCOMES, sep="\t", index=False)


def read_mtx_rows(path: pathlib.Path, rows: np.ndarray, *, cells_on_rows: bool = True) -> sparse.csr_matrix:
    """Selected cells of a large gzipped Matrix Market file, read in chunks so memory follows the kept cells only.
    rows: 0-based positions of the cells to keep; the result is in that order (cells x genes)."""
    import gzip

    with gzip.open(path, "rt") as handle:
        for line in handle:
            if not line.startswith("%"):
                n_rows, n_columns, _ = (int(value) for value in line.split())
                break
        cells, genes = (n_rows, n_columns) if cells_on_rows else (n_columns, n_rows)
        lookup = np.full(cells, -1, dtype=np.int64)
        lookup[rows] = np.arange(len(rows))
        kept_cell, kept_gene, kept_value = [], [], []
        for chunk in pd.read_csv(handle, sep=" ", header=None, chunksize=20_000_000, dtype=np.int64, comment="%",
                                 engine="c"):
            first, second, value = chunk[0].to_numpy() - 1, chunk[1].to_numpy() - 1, chunk[2].to_numpy()
            cell, gene = (first, second) if cells_on_rows else (second, first)
            mask = lookup[cell] >= 0
            kept_cell.append(lookup[cell[mask]]), kept_gene.append(gene[mask].astype(np.int32))
            kept_value.append(value[mask].astype(np.int32))
    return sparse.csr_matrix((np.concatenate(kept_value), (np.concatenate(kept_cell), np.concatenate(kept_gene))),
                             shape=(len(rows), genes), dtype=np.int32)


def target_arms(arms, targets: Targets) -> list[str]:
    """Listed targets found in the arm column, including inside a ';'-separated arm at high guide load."""
    found = set()
    for arm in pd.unique(pd.Series(arms, dtype=str)):
        for part in str(arm).split(";"):
            if part in targets.table.index:
                found.add(part)
    return sorted(found)


def pathway_arms(arms, targets: Targets) -> str:
    """';'-separated pathway genes present, for other_pathway_arms at high guide load."""
    return ";".join(arm for arm in target_arms(arms, targets) if targets.table.loc[arm, "category"] == "pathway")


def read_h5_with_guides(path) -> tuple[sparse.csr_matrix, pd.DataFrame, np.ndarray, sparse.csr_matrix, np.ndarray]:
    """One 10x h5 with a CRISPR Guide Capture block: gene counts (cells x genes), gene table, barcodes, guide counts
    (cells x guides), guide names."""
    with h5py.File(path, "r") as handle:
        group = handle["matrix"]
        whole = sparse.csc_matrix((group["data"][()], group["indices"][()], group["indptr"][()]),
                                  shape=tuple(group["shape"][()])).T.tocsr()
        types = decode(group["features/feature_type"][()])
        ids, names = decode(group["features/id"][()]), decode(group["features/name"][()])
        barcodes = decode(group["barcodes"][()])
    genes = types == "Gene Expression"
    guides = types == "CRISPR Guide Capture"
    return (whole[:, np.flatnonzero(genes)].astype(np.int32), pd.DataFrame({"id": ids[genes], "name": names[genes]}),
            barcodes, whole[:, np.flatnonzero(guides)].astype(np.int32), names[guides])


GUIDE_FDR = 0.05
GUIDE_MIN_UMI = 3


def call_guides(guide_counts: sparse.csr_matrix, guide_targets: np.ndarray) -> tuple[list[str], np.ndarray, list[str]]:
    """Guide calls for deposits without the authors' calls, one lane at a time, by the binomial test GSE249595 deposits as its
    pipeline: a guide's background frequency is its share of the lane's guide UMIs; a guide is in a cell when its UMIs exceed that
    background given the cell's guide UMIs at Benjamini-Hochberg 0.05 over the lane, with at least GUIDE_MIN_UMI UMIs. Returns per cell the ';'-joined targets, the
    number of guides and the ';'-joined guide indices. Uses guide counts only."""
    counts = guide_counts.tocoo()
    per_guide = np.asarray(guide_counts.sum(axis=0)).ravel()
    per_cell = np.asarray(guide_counts.sum(axis=1)).ravel()
    p_values = stats.binom.sf(counts.data - 1, per_cell[counts.row], per_guide[counts.col] / max(per_guide.sum(), 1))
    order = np.argsort(p_values)
    adjusted = np.minimum.accumulate((p_values[order] * len(p_values) / np.arange(1, len(p_values) + 1))[::-1])[::-1]
    passed = np.zeros(len(p_values), dtype=bool)
    passed[order] = adjusted < GUIDE_FDR
    passed &= counts.data >= GUIDE_MIN_UMI
    guide_targets = np.asarray(guide_targets, dtype=object)
    called = [[] for _ in range(guide_counts.shape[0])]
    for cell, guide in zip(counts.row[passed], counts.col[passed]):
        called[cell].append(guide)
    arms = [";".join(sorted(set(guide_targets[guides]))) if guides else "" for guides in called]
    return arms, np.array([len(guides) for guides in called], dtype=np.int64), [";".join(map(str, guides)) for guides in called]
