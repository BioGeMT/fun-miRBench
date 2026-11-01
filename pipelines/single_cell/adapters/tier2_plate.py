"""Shared by the tier 2 plate adapters: reads one STARsolo SmartSeq matrix and writes the cell store of the cells listed in
a table. cells: indexed by the manifest's cell name, columns sample, batch, arm, control_type, passed_qc; optional
replicate_unit and replicate_type. Without a unit, the passing cells of each sample and arm become three pseudo-replicates."""
from __future__ import annotations

import pathlib

import pandas as pd

from funmirbench.single_cell_store import OBS_COLUMNS, pseudo_replicate_units, read_10x_mtx, write_cell_store

GTF = pathlib.Path("/SCRATCH/dtzim01/fun-experiment-data/references/ensembl_v115/Homo_sapiens.GRCh38.115.gtf.gz")
FINAL = pathlib.Path("/SCRATCH/dtzim01/fun-experiment-data/single-cell-final/20261031")
DATA = pathlib.Path("/SCRATCH/dtzim01/single_cell_mirna_data")


def write_plate_store(dataset: str, solo: pathlib.Path, cells: pd.DataFrame, out: pathlib.Path, source: str) -> None:
    counts, features, barcodes = read_10x_mtx(solo / "matrix.mtx", solo / "features.tsv", solo / "barcodes.tsv")
    position = pd.Index(barcodes).get_indexer(cells.index)
    if (position < 0).any():
        raise ValueError(f"{dataset}: {int((position < 0).sum())} listed cells are not in {solo}")
    obs = pd.DataFrame({"dataset": dataset, "sample": cells["sample"].astype(str), "batch": cells.batch.astype(str), "arm": cells.arm.astype(str),
                        "control_type": cells.control_type.astype(str), "is_control": cells.control_type.astype(str) != "", "guide": "", "n_guides": 0,
                        "passed_qc": cells.passed_qc.astype(bool)}, index=cells.index)
    if "replicate_unit" in cells:
        obs["replicate_unit"], obs["replicate_type"] = cells.replicate_unit.astype(str), cells.replicate_type.astype(str)
    else:
        obs["replicate_unit"] = pseudo_replicate_units(obs["sample"].to_numpy(), obs.arm.to_numpy(), obs.passed_qc.to_numpy())
        obs["replicate_type"] = "technical"
    out.parent.mkdir(parents=True, exist_ok=True)
    write_cell_store(out, counts.tocsc()[:, position].T.tocsr(), obs[list(OBS_COLUMNS)], features[0], gtf=GTF, source=source, symbols=features[1])
    summary = obs.groupby(["sample", "arm"], observed=True).passed_qc.agg(["size", "sum"])
    print(f"{dataset}: {len(obs)} cells, {int(obs.passed_qc.sum())} pass QC -> {out}\n{summary.to_string()}", flush=True)
