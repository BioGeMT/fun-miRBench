"""Cell store: one h5ad per dataset, the input of funmirbench-single-cell-final, and the helpers adapters need.

X holds raw integer counts (cells x genes, sparse). var: index = feature as deposited, gene_id = unversioned Ensembl 115
gene ID where the feature maps uniquely (otherwise empty), symbol. obs holds the columns of OBS_COLUMNS for every cell.
uns["source"] records where the counts come from (authors' file or HKU reconstruction). A deposit with fractional counts
(multi-mapping reads shared between genes) keeps them, with uns["fractional_counts"] = True; they are rounded after
summing.
"""

from __future__ import annotations

import gzip
import pathlib
import re
import zlib

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

from funmirbench.single_cell_pipeline import (
    read_10x_mtx,  # noqa: F401 (adapters read mtx deposits from here, next to read_10x_h5)
    symbol_ensembl_ids,
    symbol_map_from_gtf,
    unversioned_ensembl_ids,
)

TEXT_COLUMNS = ("dataset", "sample", "batch", "replicate_unit", "replicate_type", "arm", "control_type", "guide")
OBS_COLUMNS = (*TEXT_COLUMNS, "is_control", "n_guides", "passed_qc")
NON_EMPTY_COLUMNS = ("dataset", "sample", "batch", "replicate_unit", "replicate_type", "arm")
REPLICATE_TYPES = ("biological", "lane_set", "technical")
ARM_SEPARATOR = ";"  # a cell carrying guides of several targets at high guide load: arm = "DICER1;DROSHA"
ENSEMBL_ID = r"^ENSG\d{11}(\.\d+)?$"


def gene_ids_from_gtf(gtf: pathlib.Path) -> set[str]:
    """Unversioned gene IDs of the gene records of an Ensembl GTF."""
    ids = set()
    with gzip.open(gtf, "rt") if gtf.suffix == ".gz" else gtf.open() as handle:
        for line in handle:
            fields = line.split("\t")
            if len(fields) > 8 and fields[2] == "gene":
                ids.add(re.search(r'gene_id "([^".]+)', fields[8]).group(1))
    return ids


def map_features_to_gene_ids(features: pd.Series, gtf: pathlib.Path) -> pd.Series:
    """Unversioned Ensembl gene ID per feature (Ensembl IDs or gene symbols); empty unless exactly one feature maps
    to exactly one gene of the GTF."""
    features = pd.Series(features, dtype=str).reset_index(drop=True)
    if features.str.match(ENSEMBL_ID).mean() >= 0.9:
        mapped, _ = unversioned_ensembl_ids(features)
        mapped = mapped[mapped.isin(gene_ids_from_gtf(gtf))]
    else:
        mapped, _ = symbol_ensembl_ids(features, symbol_map_from_gtf(gtf))
    mapped = mapped[~mapped.duplicated(keep=False)]
    return mapped.reindex(features.index, fill_value="")


def check_counts(counts: sparse.spmatrix, *, fractional: bool = False) -> None:
    data = counts.tocsr().data
    if data.size and data.min() < 0:
        raise ValueError(f"counts must not be negative (minimum {data.min()})")
    if not fractional and not np.issubdtype(data.dtype, np.integer) and not np.array_equal(data, np.round(data)):
        raise ValueError("counts must be raw integers; found fractional values")


def validate_cell_store(store: ad.AnnData, gtf_gene_ids: set[str] | None = None) -> None:
    """Raises ValueError for anything funmirbench-single-cell-final relies on."""
    if not sparse.issparse(store.X):
        raise ValueError("X must be sparse")
    check_counts(store.X, fractional=bool(store.uns.get("fractional_counts", False)))
    missing = [column for column in OBS_COLUMNS if column not in store.obs]
    if missing:
        raise ValueError(f"obs is missing columns: {missing}")
    obs = store.obs
    for column in OBS_COLUMNS:
        if obs[column].isna().any():
            raise ValueError(f"obs.{column} has missing values")
    for column in NON_EMPTY_COLUMNS:
        if (obs[column].astype(str) == "").any():
            raise ValueError(f"obs.{column} has empty values")
    if obs.dataset.nunique() != 1:
        raise ValueError(f"a cell store holds one dataset, found {sorted(obs.dataset.unique())}")
    unknown = set(obs.replicate_type.astype(str)) - set(REPLICATE_TYPES)
    if unknown:
        raise ValueError(f"replicate_type must be one of {REPLICATE_TYPES}, found {sorted(unknown)}")
    for column in ("is_control", "passed_qc"):
        if obs[column].dtype != bool:
            raise ValueError(f"obs.{column} must be boolean")
    if not np.issubdtype(obs.n_guides.dtype, np.integer) or (obs.n_guides < 0).any():
        raise ValueError("obs.n_guides must be non-negative integers")
    if ((obs.control_type.astype(str) != "") != obs.is_control.to_numpy()).any():
        raise ValueError("control_type must be set for control cells and empty for all others")
    if not {"gene_id", "symbol"} <= set(store.var):
        raise ValueError("var needs gene_id and symbol")
    gene_ids = store.var.gene_id.astype(str)
    mapped = gene_ids[gene_ids != ""]
    if not mapped.str.fullmatch(r"ENSG\d{11}").all():
        raise ValueError(f"var.gene_id must be unversioned Ensembl gene IDs, for example "
                         f"{mapped[~mapped.str.fullmatch(r'ENSG[0-9]{11}')].head(3).tolist()}")
    if mapped.duplicated().any():
        raise ValueError(f"var.gene_id must map uniquely, for example {mapped[mapped.duplicated()].head(3).tolist()}")
    if gtf_gene_ids is not None and not set(mapped) <= gtf_gene_ids:
        raise ValueError(f"var.gene_id not in the GTF, for example {sorted(set(mapped) - gtf_gene_ids)[:3]}")
    if not store.uns.get("source"):
        raise ValueError('uns["source"] must record where the counts come from')


def write_cell_store(path: pathlib.Path, counts, obs: pd.DataFrame, features: pd.Series, *, gtf: pathlib.Path,
                     source: str, symbols: pd.Series | None = None, gene_ids: pd.Series | None = None,
                     fractional: bool = False) -> ad.AnnData:
    """Validates and writes one cell store. counts: cells x features; features: as deposited (Ensembl IDs or gene
    symbols), mapped to gene_id through the GTF unless gene_ids gives them (deposits that are mostly non-Ensembl
    features); symbols defaults to the features; fractional keeps non-integer counts as deposited."""
    counts = sparse.csr_matrix(counts)
    check_counts(counts, fractional=fractional)
    features = pd.Series(features, dtype=str).reset_index(drop=True)
    symbols = features if symbols is None else pd.Series(symbols, dtype=str).reset_index(drop=True)
    mapped = map_features_to_gene_ids(features, gtf) if gene_ids is None else pd.Series(gene_ids, dtype=str).fillna("")
    var = pd.DataFrame({"gene_id": mapped.to_numpy(), "symbol": symbols.to_numpy()},
                       index=features.to_numpy())
    obs = obs.copy()
    for column in TEXT_COLUMNS:
        if column in obs:
            obs[column] = obs[column].astype("category")
    uns = {"source": source, "fractional_counts": True} if fractional else {"source": source}
    store = ad.AnnData(X=counts if fractional else counts.astype(np.int32), obs=obs, var=var, uns=uns)
    validate_cell_store(store)
    store.write_h5ad(path)
    return store


def read_cell_store(path: pathlib.Path) -> ad.AnnData:
    store = ad.read_h5ad(path)
    validate_cell_store(store)
    return store


def read_10x_h5(path: pathlib.Path):
    """Counts (genes x cells, CSC), the features table and the barcodes of one 10x HDF5 file, in the form of
    read_10x_mtx. Both layouts: a group "matrix" (CellRanger 3 and later) or one group per genome (CellRanger 2)."""
    def text(dataset) -> np.ndarray:
        return np.array([value.decode() if isinstance(value, bytes) else value for value in dataset[()]], dtype=object)

    with h5py.File(path, "r") as handle:
        layout_3 = "matrix" in handle
        group = handle["matrix"] if layout_3 else handle[next(iter(handle))]
        counts = sparse.csc_matrix((group["data"][()], group["indices"][()], group["indptr"][()]),
                                   shape=tuple(group["shape"][()]))
        barcodes = text(group["barcodes"])
        if layout_3:
            features = pd.DataFrame({0: text(group["features/id"]), 1: text(group["features/name"]),
                                     2: text(group["features/feature_type"])})
        else:
            features = pd.DataFrame({0: text(group["genes"]), 1: text(group["gene_names"])})
    if 2 in features:
        keep = (features[2] == "Gene Expression").to_numpy()
        counts, features = counts[keep, :], features[keep].reset_index(drop=True)
    return counts, features, barcodes


def raw_counts_from_lognorm(lognorm, library_size, *, scale: float = 1e4, tolerance: float = 0.1) -> sparse.csr_matrix:
    """HKU files hold log-normalised values (cells x genes): raw counts are round(expm1(X) * library_size / scale).
    Raises ValueError when the reconstruction is not integer within tolerance, which means the library sizes or the
    scale do not belong to these values."""
    matrix = sparse.csr_matrix(lognorm, dtype=np.float64)
    library_size = np.asarray(library_size, dtype=np.float64)
    if library_size.shape != (matrix.shape[0],) or not (library_size > 0).all():
        raise ValueError("library_size needs one positive value per cell")
    per_entry = np.repeat(library_size, np.diff(matrix.indptr))
    values = np.expm1(matrix.data) * per_entry / scale
    rounded = np.round(values)
    if values.size and np.abs(values - rounded).max() > tolerance:
        raise ValueError(f"reconstructed counts are not integer: largest distance {np.abs(values - rounded).max():.3f} "
                         f"to an integer (tolerance {tolerance})")
    out = sparse.csr_matrix((rounded.astype(np.int32), matrix.indices, matrix.indptr), shape=matrix.shape)
    out.eliminate_zeros()
    return out


def pseudo_replicate_units(samples, arms, passed, *, n_replicates: int = 3, seed: int = 0) -> np.ndarray:
    """replicate_unit for one sample per arm: the passing cells of each sample and arm split at random into equal
    parts, named <sample>_pr<i>. Cells that fail QC keep <sample> (they are never summed). A group's split depends
    only on the seed, its name and its number of cells. The adapter records replicate_type = technical."""
    samples, arms = np.asarray(samples, dtype=object), np.asarray(arms, dtype=object)
    passed = np.asarray(passed, dtype=bool)
    units = samples.astype(object).copy()
    for sample, arm in sorted(set(zip(samples[passed], arms[passed]))):
        cells = np.flatnonzero((samples == sample) & (arms == arm) & passed)
        rng = np.random.default_rng([seed, zlib.crc32(f"{sample}\t{arm}".encode())])
        parts = rng.permutation(len(cells)) % n_replicates
        units[cells] = [f"{sample}_pr{part + 1}" for part in parts]
    return units


def assign_lane_sets(lanes, *, lanes_per_set: int = 8) -> np.ndarray:
    """Lane set per cell from integer lane numbers: consecutive lanes in sets of lanes_per_set; a last incomplete set
    joins the one before."""
    lanes = np.asarray(lanes)
    if not np.issubdtype(lanes.dtype, np.integer):
        raise ValueError("lane numbers must be integers (text lane names do not sort consecutively)")
    observed, position = np.unique(lanes, return_inverse=True)
    n_sets = max(1, len(observed) // lanes_per_set)
    sets = np.minimum(position // lanes_per_set, n_sets - 1)
    return np.array([f"lane_set_{value + 1}" for value in sets], dtype=object)
