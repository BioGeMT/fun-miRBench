"""Shared by the adapters of datasets whose arms are samples or hashed arms: one library is read, QC thresholds from
the authors or the deposit's own filter are applied, the libraries are aligned and one cell store is written."""

from __future__ import annotations

import argparse
import os
import pathlib
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import sparse

from funmirbench.single_cell_store import OBS_COLUMNS, pseudo_replicate_units, write_cell_store

REPO = pathlib.Path(__file__).resolve().parents[3]
DATA = pathlib.Path(os.environ.get("SINGLE_CELL_DATA", REPO / "data/single_cell/raw"))
GTF = pathlib.Path(os.environ.get("SINGLE_CELL_GTF", DATA / "refs/Homo_sapiens.GRCh38.115.gtf.gz"))
ENSEMBL_PREFIX = "ENSG"


@dataclass
class Library:
    """counts: genes x cells; features: id in the first column, symbol in the second if there is one. arm: one name per
    library, or one per cell for hashed arms (None drops the cell from the store). unit: replicate unit of real
    replicates, one per library; None asks for three pseudo-replicates."""
    sample: str
    batch: str
    arm: object
    counts: sparse.spmatrix
    features: pd.DataFrame
    barcodes: np.ndarray
    control_type: object = ""
    unit: str | None = None
    unit_type: str = "biological"
    pseudo_group: str | None = None  # pseudo-replicates split the cells of this group and arm; default the sample


def feature_keys(features: pd.DataFrame) -> pd.Series:
    keys = features[0].astype(str)
    if keys.str.startswith(ENSEMBL_PREFIX).mean() >= 0.9:
        keys = keys.str.replace(r"\.\d+$", "", regex=True)
    return keys


def symbols_of(features: pd.DataFrame) -> pd.Series:
    return features[1].astype(str) if 1 in features else features[0].astype(str)


def align(libraries: list[Library]) -> tuple[list[sparse.csc_matrix], pd.DataFrame]:
    """The features every library lists, in the order of the first; a feature repeated within a library is dropped,
    since it cannot be assigned to one gene."""
    kept = []
    for library in libraries:
        keys = feature_keys(library.features)
        kept.append(keys[~keys.duplicated(keep=False)])
    common = set(kept[0])
    for keys in kept[1:]:
        common &= set(keys)
    first = kept[0][kept[0].isin(common)]
    order = pd.Index(first.to_numpy())
    matrices = []
    for library, keys in zip(libraries, kept):
        position = pd.Series(keys.index.to_numpy(), index=keys.to_numpy()).reindex(order).to_numpy()
        matrices.append(sparse.csc_matrix(library.counts)[position, :])
    table = pd.DataFrame({"key": order, "symbol": symbols_of(libraries[0].features).iloc[first.index].to_numpy()})
    for library, keys in zip(libraries, kept):
        dropped = len(keys) - len(order)
        if dropped:
            print(f"{library.sample}: {dropped} features not in every library left out", flush=True)
    return matrices, table


def build_store(dataset: str, libraries: list[Library], *, qc: dict, qc_source: str, source: str, out: pathlib.Path,
                mito_prefix: str = "MT-", fractional: bool = False) -> None:
    """qc: min_genes, max_genes, min_counts, max_counts, max_mito (fraction); an absent key is not applied. Writes all
    cells of the listed arms; cells failing a threshold stay with passed_qc = False."""
    matrices, features = align(libraries)
    mito = features.symbol.str.upper().str.startswith(mito_prefix).to_numpy()
    if "max_mito" in qc and not mito.any():
        raise ValueError(f"{dataset}: max_mito given but no symbol starts with {mito_prefix}")
    frames, kept_matrices = [], []
    for library, matrix in zip(libraries, matrices):
        n_cells = matrix.shape[1]
        arm = np.full(n_cells, library.arm, dtype=object) if isinstance(library.arm, str) else np.asarray(library.arm, dtype=object)
        keep = np.array([value is not None for value in arm])
        control = np.full(n_cells, library.control_type, dtype=object) if isinstance(library.control_type, str) \
            else np.asarray(library.control_type, dtype=object)
        totals = np.asarray(matrix.sum(axis=0)).ravel()
        genes = np.diff(matrix.tocsc().indptr)
        mito_counts = np.asarray(matrix[np.flatnonzero(mito), :].sum(axis=0)).ravel()
        fraction = np.divide(mito_counts, totals, out=np.zeros(n_cells), where=totals > 0)
        passed = np.ones(n_cells, dtype=bool)
        for key, values, op in (("min_genes", genes, np.greater_equal), ("max_genes", genes, np.less_equal),
                                ("min_counts", totals, np.greater_equal), ("max_counts", totals, np.less_equal),
                                ("max_mito", fraction, np.less_equal)):
            if key in qc:
                passed &= op(values, qc[key])
        frames.append(pd.DataFrame({
            "dataset": dataset, "sample": library.sample, "batch": library.batch, "arm": arm, "control_type": control,
            "is_control": control != "", "guide": "", "n_guides": 0, "passed_qc": passed,
            "replicate_unit": library.unit or library.sample, "replicate_type": library.unit_type if library.unit else "technical",
            "pseudo_group": library.pseudo_group or library.sample,
        }, index=[f"{library.sample}:{barcode}" for barcode in library.barcodes])[keep])
        kept_matrices.append(matrix[:, np.flatnonzero(keep)])
        print(f"{library.sample}: {n_cells} cells in the deposit, {keep.sum()} in the listed arms, "
              f"{int(passed[keep].sum())} pass QC", flush=True)
    obs = pd.concat(frames)
    split = obs.replicate_unit.to_numpy() == obs.pseudo_group.to_numpy()
    split |= obs.replicate_unit.to_numpy() == obs["sample"].to_numpy()
    split &= (obs.replicate_type == "technical").to_numpy()
    if split.any():
        obs.loc[split, "replicate_unit"] = pseudo_replicate_units(
            obs.pseudo_group.to_numpy()[split], obs.arm.to_numpy()[split], obs.passed_qc.to_numpy()[split])
    counts = sparse.hstack(kept_matrices, format="csr").T.tocsr()
    out.parent.mkdir(parents=True, exist_ok=True)
    write_cell_store(out, counts, obs[list(OBS_COLUMNS)], features.key, gtf=GTF, source=f"{source}; cell QC: {qc_source}",
                     symbols=features.symbol, fractional=fractional)
    print(f"{dataset}: {obs.shape[0]} cells, {int(obs.passed_qc.sum())} pass QC, {features.shape[0]} features -> {out}", flush=True)


def read_dense_counts(path: pathlib.Path, sep: str = "\t") -> tuple[sparse.csc_matrix, pd.DataFrame, np.ndarray]:
    """Text table with the gene in the first column and one column per cell, in the form of read_10x_mtx."""
    table = pd.read_csv(path, sep=sep, index_col=0)
    counts = sparse.csc_matrix(table.to_numpy())
    return counts, pd.DataFrame({0: table.index.astype(str)}), table.columns.astype(str).to_numpy()


def store_path() -> pathlib.Path:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=pathlib.Path, required=True)
    return parser.parse_args().out
