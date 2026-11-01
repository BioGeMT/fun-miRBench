"""Single-cell miRNA perturbation experiments to pseudobulk counts for the bulk DE pipeline.

Each sample (one culture, donor, day or hashed arm) is read from its 10x count files, optionally restricted to the
cells a label table assigns to it (hashtag calls of a pooled lane), cells are kept by the configured QC thresholds,
and the counts of the kept cells are summed per sample. The summed counts and a
count_matrix-mode config are written for funmirbench-experiments, which fits DESeq2 exactly as for bulk RNA-seq
(design ~ batch + condition when the samples carry a batch, for example the day of a paired design).
"""

from __future__ import annotations

import argparse
import gzip
import logging
import pathlib
import re
import sys

import numpy as np
import pandas as pd
import requests
import yaml
from scipy import io as scipy_io
from scipy import sparse

from funmirbench.experiments_pipeline import (
    ensure_clean_dir,
    load_yaml,
    normalize_space,
    repo_root,
    require_fields,
    resolve_path,
    run_ingestion_config,
    utc_now_stamp,
    validate_mirna_name,
    write_json,
)
from funmirbench.logger import parse_log_level, setup_logging

logger = logging.getLogger(__name__)

GROUPS = ("control", "treated")
FILE_KEYS = ("matrix", "features", "barcodes")
QC_KEYS = ("min_genes", "min_counts", "max_counts", "max_mito_fraction")


def cached_download(url: str, dest: pathlib.Path, *, force: bool) -> pathlib.Path:
    if dest.exists() and not force:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    response = requests.get(url, timeout=120)
    response.raise_for_status()
    dest.write_bytes(response.content)
    return dest


def read_10x_mtx(matrix: pathlib.Path, features: pathlib.Path, barcodes: pathlib.Path):
    """Counts (genes x cells, CSC), the features table and the barcodes of one 10x sample. Only rows of type
    "Gene Expression" are kept when the features file has a type column (antibody or hashtag rows are dropped)."""
    counts = sparse.csc_matrix(scipy_io.mmread(str(matrix)))
    feature_table = pd.read_csv(features, sep="\t", header=None, dtype=str, compression="infer")
    cell_barcodes = pd.read_csv(barcodes, sep="\t", header=None, dtype=str, compression="infer")[0].to_numpy()
    if counts.shape != (len(feature_table), len(cell_barcodes)):
        raise ValueError(f"{matrix}: matrix shape {counts.shape} does not match {len(feature_table)} features and "
                         f"{len(cell_barcodes)} barcodes")
    if feature_table.shape[1] >= 3:
        keep = (feature_table[2] == "Gene Expression").to_numpy()
        counts, feature_table = counts[keep, :], feature_table[keep].reset_index(drop=True)
    return counts, feature_table, cell_barcodes


def cell_qc(counts: sparse.csc_matrix, symbols: pd.Series, mito_prefix: str, qc: dict) -> pd.DataFrame:
    """Per-cell total counts, detected genes and mitochondrial fraction, and whether the cell passes every
    configured threshold."""
    totals = np.asarray(counts.sum(axis=0)).ravel()
    detected = np.diff(counts.indptr)  # nonzero rows per column of a CSC matrix
    mito = symbols.str.upper().str.startswith(mito_prefix.upper()).to_numpy()
    if not mito.any():
        raise ValueError(f"no gene symbol starts with {mito_prefix!r}; check the features file or qc.mito_prefix")
    mito_counts = np.asarray(counts[mito, :].sum(axis=0)).ravel()
    fraction = np.divide(mito_counts, totals, out=np.zeros_like(totals, dtype=float), where=totals > 0)
    table = pd.DataFrame({"total_counts": totals, "detected_genes": detected, "mito_fraction": fraction})
    passed = np.ones(len(table), dtype=bool)
    if "min_genes" in qc:
        passed &= detected >= qc["min_genes"]
    if "min_counts" in qc:
        passed &= totals >= qc["min_counts"]
    if "max_counts" in qc:
        passed &= totals <= qc["max_counts"]
    if "max_mito_fraction" in qc:
        passed &= fraction <= qc["max_mito_fraction"]
    table["passed"] = passed
    return table


def unversioned_ensembl_ids(ids: pd.Series) -> tuple[pd.Series, list[str]]:
    """Unversioned IDs and the features that are not human Ensembl genes (transgenes such as a reporter cassette),
    which are left out of the pseudobulk. Mostly non-Ensembl features mean the file uses another ID scheme."""
    ids = ids.str.replace(r"\.\d+$", "", regex=True)
    ensembl = ids.str.match(r"^ENSG\d{11}$")
    if ensembl.mean() < 0.9:
        raise ValueError(f"features are not human Ensembl gene IDs (for example {ids[~ensembl].head(3).tolist()}); "
                         "set feature_ids: symbol and gene_symbols_gtf for files with gene symbols")
    return ids[ensembl], ids[~ensembl].tolist()


def symbol_map_from_gtf(gtf: pathlib.Path) -> dict[str, str]:
    """Gene symbol -> Ensembl gene ID from the gene records of an Ensembl GTF, for symbols of exactly one gene."""
    rows = []
    with gzip.open(gtf, "rt") if gtf.suffix == ".gz" else gtf.open() as handle:
        for line in handle:
            fields = line.split("\t")
            if len(fields) > 8 and fields[2] == "gene":
                gene_id = re.search(r'gene_id "([^"]+)"', fields[8])
                name = re.search(r'gene_name "([^"]+)"', fields[8])
                if gene_id and name:
                    rows.append((name.group(1), gene_id.group(1)))
    table = pd.DataFrame(rows, columns=["symbol", "gene_id"]).drop_duplicates()
    table = table[~table.symbol.duplicated(keep=False)]
    return dict(zip(table.symbol, table.gene_id))


def symbol_ensembl_ids(symbols: pd.Series, symbol_map: dict[str, str]) -> tuple[pd.Series, list[str]]:
    """Ensembl IDs of the features whose symbol names exactly one gene; the others are left out."""
    ids = symbols.map(symbol_map)
    if ids.notna().mean() < 0.5:
        raise ValueError(f"fewer than half of the feature symbols match the GTF (for example "
                         f"{symbols[ids.isna()].head(3).tolist()})")
    return ids.dropna(), symbols[ids.isna()].tolist()


def labelled_barcodes(cell_labels: dict) -> np.ndarray:
    """Barcodes that a label table (for example hashtag calls) assigns to one sample. The label is the label
    column itself, or its first group under label_regex (for example the tag before an underscore)."""
    table = pd.read_csv(cell_labels["path"], sep=None, engine="python", dtype=str, compression="infer")
    labels = table[cell_labels["label_column"]]
    if cell_labels.get("label_regex"):
        labels = labels.str.extract(cell_labels["label_regex"], expand=False)
    return table.loc[labels == cell_labels["label"], cell_labels["barcode_column"]].to_numpy()


def load_samples(config: dict, *, config_path: pathlib.Path, repo: pathlib.Path, force: bool) -> list[dict]:
    """Validated sample entries with local file paths; files given as URLs are downloaded into data/single_cell/raw."""
    samples = config.get("samples") or []
    if not samples:
        raise ValueError("Config needs a samples list.")
    seen, out = set(), []
    for sample in samples:
        sample_id = normalize_space(sample.get("sample_id", ""))
        group = normalize_space(sample.get("group", ""))
        if not sample_id or sample_id in seen:
            raise ValueError(f"Missing or duplicated sample_id {sample_id!r}.")
        if group not in GROUPS:
            raise ValueError(f"{sample_id}: group must be one of {GROUPS}, not {group!r}.")
        seen.add(sample_id)
        entry = {"sample_id": sample_id, "group": group, "batch": normalize_space(sample.get("batch", ""))}
        if sample.get("cell_labels"):
            labels = dict(sample["cell_labels"])
            missing = {"path", "barcode_column", "label_column", "label"} - set(labels)
            if missing:
                raise ValueError(f"{sample_id}: cell_labels needs {sorted(missing)}.")
            labels["path"] = resolve_path(str(labels["path"]), root=config_path.parent, repo=repo)
            entry["cell_labels"] = labels
        for key in FILE_KEYS:
            value = normalize_space(sample.get(key, ""))
            if not value:
                raise ValueError(f"{sample_id}: {key} is required.")
            if value.startswith(("http://", "https://", "ftp://")):
                dest = repo / "data" / "single_cell" / "raw" / config["dataset_id"] / value.rsplit("/", 1)[-1]
                entry[key] = cached_download(value, dest, force=force)
            else:
                entry[key] = resolve_path(value, root=config_path.parent, repo=repo)
            if not entry[key].exists():
                raise FileNotFoundError(f"{sample_id}: {entry[key]} not found.")
        out.append(entry)
    for group in GROUPS:
        if not any(s["group"] == group for s in out):
            raise ValueError(f"No {group} samples.")
    batches = [s["batch"] for s in out]
    if any(batches) and not all(batches):
        raise ValueError("Either every sample has a batch or none does.")
    return out


def pseudobulk(samples: list[dict], qc: dict, mito_prefix: str, min_cells: int,
               symbol_map: dict[str, str] | None = None, pseudo_replicates: int = 1):
    """Summed counts of the passing cells (unversioned Ensembl gene IDs x columns), the per-cell QC table, the
    per-sample summary and the design (column, sample, group, batch). With symbol_map the features are gene
    symbols, mapped to Ensembl IDs. Only genes listed in every sample's features file are kept. With
    pseudo_replicates k > 1 the passing cells of each sample are split at random into k columns (technical
    replicates: their spread is cell sampling, not culture-to-culture variation)."""
    columns, cell_tables, summary, files, design = {}, [], [], {}, []
    rng = np.random.default_rng(0)
    for sample in samples:
        key = tuple(str(sample[k]) for k in FILE_KEYS)
        if key not in files:  # the samples of a hashed lane share one file
            files[key] = read_10x_mtx(sample["matrix"], sample["features"], sample["barcodes"])
        counts, features, barcodes = files[key]
        in_file = len(barcodes)
        if "cell_labels" in sample:
            selected = np.isin(barcodes, labelled_barcodes(sample["cell_labels"]))
            if not selected.any():
                raise ValueError(f"{sample['sample_id']}: no barcode of the label table is in the matrix.")
            counts, barcodes = counts[:, selected], barcodes[selected]
        cells = cell_qc(counts, features[1], mito_prefix, qc)
        cells.insert(0, "barcode", barcodes)
        cells.insert(0, "sample_id", sample["sample_id"])
        cell_tables.append(cells)
        kept = int(cells.passed.sum())
        if kept < min_cells:
            raise ValueError(f"{sample['sample_id']}: {kept} cells pass QC, fewer than {min_cells}.")
        passed = np.flatnonzero(cells.passed.to_numpy())
        summed = np.asarray(counts[:, passed].sum(axis=1)).ravel()
        if symbol_map is None:
            gene_ids, other = unversioned_ensembl_ids(features[0])
        else:
            gene_ids, other = symbol_ensembl_ids(features[0], symbol_map)
        part = rng.permutation(kept) % pseudo_replicates
        if pseudo_replicates > 1:
            cells.loc[cells.passed, "pseudo_replicate"] = part + 1
        for i in range(pseudo_replicates):
            column = sample["sample_id"] if pseudo_replicates == 1 else f"{sample['sample_id']}_pr{i + 1}"
            part_sum = np.asarray(counts[:, passed[part == i]].sum(axis=1)).ravel()
            columns[column] = pd.Series(part_sum[gene_ids.index], index=gene_ids.to_numpy()).groupby(level=0).sum()
            design.append({"column": column, "sample_id": sample["sample_id"], "group": sample["group"],
                           "batch": sample["batch"]})
        summary.append({"sample_id": sample["sample_id"], "group": sample["group"], "batch": sample["batch"],
                        "cells_in_file": in_file, "cells_selected": len(cells), "cells_passed": kept,
                        "features_left_out": len(other), "features_left_out_examples": ";".join(other[:5]),
                        "summed_counts": int(round(summed.sum())),
                        "fractional_counts": bool(not np.array_equal(counts.data, np.round(counts.data))),
                        "median_counts_per_cell": float(cells.total_counts[cells.passed].median())})
    # samples quantified against different annotations: a gene missing from a file is not measured there, not zero
    common = sorted(set.intersection(*(set(column.index) for column in columns.values())))
    first_column = {row["sample_id"]: row["column"] for row in reversed(design)}
    for row in summary:
        row["genes_not_in_every_file"] = len(columns[first_column[row["sample_id"]]]) - len(common)
    # sums rounded: some pipelines deposit fractional counts (multi-mapping reads shared between genes)
    matrix = pd.DataFrame({name: column.reindex(common) for name, column in columns.items()}).round().astype(np.int64)
    matrix.index.name = "gene_id"
    return matrix, pd.concat(cell_tables, ignore_index=True), pd.DataFrame(summary), pd.DataFrame(design)


def experiments_config(config: dict, design: pd.DataFrame, counts_path: pathlib.Path) -> dict:
    """count_matrix-mode config for funmirbench-experiments; the batch of each column when the design has one."""
    comparison = {
        "control_columns": design.column[design.group == "control"].tolist(),
        "treated_columns": design.column[design.group == "treated"].tolist(),
    }
    if design.batch.nunique() > 1:
        comparison["batch_by_sample"] = dict(zip(design.column, design.batch))
    out = {key: config[key] for key in ("dataset_id", "mirna_name", "experiment_type", "gse") if key in config}
    out["source"] = {"mode": "count_matrix", "count_matrix_path": str(counts_path), "gene_id_column": "gene_id"}
    out["comparison"] = comparison
    out["metadata"] = {"method": "scRNA-seq, pseudobulk", **(config.get("metadata") or {})}
    return out


def run_single_cell_config(config_path: pathlib.Path, repo: pathlib.Path | None = None, *, force: bool = False,
                           run_de: bool = False) -> dict:
    config_path = config_path.expanduser().resolve()
    repo = (repo or repo_root()).resolve()
    config = load_yaml(config_path)
    require_fields(config, ["dataset_id", "mirna_name", "experiment_type"])
    validate_mirna_name(config, repo=repo)
    qc = config.get("qc") or {}
    unknown = set(qc) - set(QC_KEYS) - {"mito_prefix"}
    if unknown:
        raise ValueError(f"Unknown qc keys: {sorted(unknown)}")
    min_cells = int(config.get("min_cells_per_sample", 20))
    pseudo_replicates = int(config.get("pseudo_replicates", 1))
    if pseudo_replicates < 1:
        raise ValueError("pseudo_replicates must be at least 1.")
    feature_ids = config.get("feature_ids", "ensembl")
    if feature_ids not in ("ensembl", "symbol"):
        raise ValueError("feature_ids must be ensembl or symbol.")
    symbol_map = None
    if feature_ids == "symbol":
        if not config.get("gene_symbols_gtf"):
            raise ValueError("feature_ids: symbol needs gene_symbols_gtf (the Ensembl GTF of the benchmark).")
        symbol_map = symbol_map_from_gtf(resolve_path(config["gene_symbols_gtf"], root=config_path.parent, repo=repo))

    run_dir = repo / "pipelines" / "single_cell" / "runs" / f"{utc_now_stamp()}_{config['dataset_id']}"
    ensure_clean_dir(run_dir, force=force)
    (run_dir / "config.yaml").write_text(config_path.read_text(encoding="utf-8"), encoding="utf-8")
    logger.info("Run dir: %s", run_dir)

    samples = load_samples(config, config_path=config_path, repo=repo, force=force)
    matrix, cells, summary, design = pseudobulk(samples, qc, qc.get("mito_prefix", "MT-"), min_cells, symbol_map,
                                                pseudo_replicates)
    counts_path = run_dir / "pseudobulk_counts.tsv.gz"
    matrix.to_csv(counts_path, sep="\t")
    with gzip.open(run_dir / "cell_qc.tsv.gz", "wt") as handle:
        cells.to_csv(handle, sep="\t", index=False)
    summary.to_csv(run_dir / "sample_summary.tsv", sep="\t", index=False)
    de_config_path = run_dir / "experiments_config.yaml"
    design.to_csv(run_dir / "design.tsv", sep="\t", index=False)
    de_config_path.write_text(yaml.safe_dump(experiments_config(config, design, counts_path), sort_keys=False),
                              encoding="utf-8")
    manifest = {"dataset_id": config["dataset_id"], "config_path": str(config_path), "qc": qc,
                "min_cells_per_sample": min_cells, "pseudo_replicates": pseudo_replicates,
                "genes": int(matrix.shape[0]),
                "samples": summary.to_dict(orient="records"),
                "inputs": {s["sample_id"]: {k: str(s[k]) for k in FILE_KEYS} for s in samples},
                "experiments_config": str(de_config_path)}
    result = {"run_dir": str(run_dir), "counts_path": str(counts_path), "experiments_config": str(de_config_path)}
    if run_de:
        de = run_ingestion_config(de_config_path, repo=repo, force=force)
        manifest["de_run_dir"], manifest["de_table_path"] = de["run_dir"], de["de_table_path"]
        result.update(de_run_dir=de["run_dir"], de_table_path=de["de_table_path"])
    write_json(manifest, run_dir / "run_manifest.json")
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Pseudobulk one single-cell experiment config.")
    parser.add_argument("--config", type=pathlib.Path, required=True)
    parser.add_argument("--run-de", action="store_true")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"])
    return parser.parse_args(sys.argv[1:])


def main() -> None:
    args = parse_args()
    setup_logging(parse_log_level(args.log_level))
    result = run_single_cell_config(args.config, force=args.force, run_de=args.run_de)
    for key, value in result.items():
        logger.info("%s: %s", key, value)


if __name__ == "__main__":
    main()
