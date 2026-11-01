"""miR-126 sensor in AML PDX: GFP-low sorted cells (high miR-126 activity) against bulk human CD45+ cells, hashed into one
capture, per donor and treatment. Two captures of each pool are technical units, capture is the batch. Cells in the
authors' metadata passed their QC; the metadata itself was read from the deposit's rds to /SCRATCH/dtzim01/tmp."""

import numpy as np
import pandas as pd
import pathlib

from _sample_arms import DATA, Library, build_store, store_path
from funmirbench.single_cell_store import read_10x_mtx

RAW = DATA / "GSE185834" / "extracted"
cells = pd.read_csv(DATA / "GSE185834" / "gse185834_cells.tsv", sep="\t").set_index("cell")
libraries = []
for path in sorted(RAW.glob("GSM*_barcodes.tsv.gz")):
    stem = str(path).removesuffix("_barcodes.tsv.gz")
    capture = stem.rsplit("_", 1)[1]
    counts, features, barcodes = read_10x_mtx(pathlib.Path(f"{stem}_matrix.mtx.gz"), pathlib.Path(f"{stem}_features.tsv.gz"), pathlib.Path(f"{stem}_barcodes.tsv.gz"))
    names = np.array([f"{capture}_{barcode}" for barcode in barcodes])
    wanted = np.flatnonzero(cells.index.reindex(names)[1] >= 0)
    if not len(wanted):
        raise ValueError(f"{capture}: no barcode in the metadata")
    info = cells.loc[names[wanted]]
    arm = (info.donor + "_" + info.cond + "_" + info.gate).tolist()
    control_type = ["bulk_sorted" if gate == "Bulk" else "" for gate in info.gate]
    libraries.append(Library(capture, capture, arm, counts[:, wanted], features, barcodes[wanted], control_type, unit=capture, unit_type="technical"))
build_store("GSE185834", libraries, qc={}, qc_source="authors' metadata: only cells that passed their filter are listed", source="authors' matrices and metadata, GEO GSE185834", out=store_path())
