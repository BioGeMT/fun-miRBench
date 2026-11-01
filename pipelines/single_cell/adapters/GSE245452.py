"""del(5q) myelodysplastic syndrome CD34+ cells: del(5q) against normal cells of the same patient, as called by the authors
in their per-cell metadata (cell_5q, from RNA-based copy number). Patient is the batch and the real replicate. Run with
openpyxl available. QC: deposit's per-cell metadata (cells listed passed the authors' filter)."""

import re
import numpy as np
import pandas as pd
import pathlib

from _sample_arms import DATA, Library, build_store, store_path
from funmirbench.single_cell_store import read_10x_mtx

RAW = DATA / "GSE245452" / "extracted"
SAMPLES = [("GSM7842961", "SMD34459", "dx"), ("GSM7842962", "SMD35109", "dx"), ("GSM7842963", "SMD35303", "dx"), ("GSM7842964", "SMD37209", "dx"),
           ("GSM7842965", "FS-0406", "lenalidomide"), ("GSM7842966", "FS-0634", "lenalidomide"), ("GSM7842967", "SMD132114579", "lenalidomide")]
libraries = []
for gsm, patient, group in SAMPLES:
    stem = next(RAW.glob(f"{gsm}_*_barcodes.tsv.gz")).as_posix().removesuffix("_barcodes.tsv.gz")
    counts, features, barcodes = read_10x_mtx(pathlib.Path(f"{stem}_matrix.mtx.gz"), pathlib.Path(f"{stem}_features.tsv.gz"), pathlib.Path(f"{stem}_barcodes.tsv.gz"))
    meta = pd.read_excel(next(RAW.glob(f"{gsm}_Metadata_*.xlsx")))
    if "cell_passed_QC" in meta and not meta.cell_passed_QC.astype(bool).all():
        meta = meta[meta.cell_passed_QC.astype(bool)]
    called = dict(zip(meta.iloc[:, 0].astype(str).str.extract(r"([ACGT]{16})")[0], meta.cell_5q))
    calls = [called.get(re.search(r"[ACGT]{16}", barcode).group(0)) for barcode in barcodes]
    wanted = np.flatnonzero([call is not None for call in calls])
    if not len(wanted):
        raise ValueError(f"{gsm}: no barcode matched the metadata")
    arm = [f"{calls[i]}_{group}" for i in wanted]
    control_type = ["same_patient_normal_cells" if calls[i] == "normal" else "" for i in wanted]
    libraries.append(Library(patient, patient, arm, counts[:, wanted], features, barcodes[wanted], control_type, unit=patient))
build_store("GSE245452", libraries, qc={}, qc_source="authors' per-cell metadata", source="authors' matrices and per-cell metadata, GEO GSE245452", out=store_path())
