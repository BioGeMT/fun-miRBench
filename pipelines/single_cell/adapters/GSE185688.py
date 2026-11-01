"""miR-365 in iPSC cardiomyocytes, 25 nM, 48 h. Healthy line: mimic-365 against antimiR-365 in separate libraries, no
neutral control. LQT1 line: antimiR-365 and antimiR-Ctrl hashed into one lane. One arm each, three pseudo-replicates.
QC: authors (10.1038/s41467-021-27856-7)."""

import numpy as np
import pandas as pd
from _sample_arms import DATA, Library, build_store, store_path
from funmirbench.single_cell_store import read_10x_mtx

RAW = DATA / "GSE185688" / "matrices"
libraries = []
for stem, sample, arm, control_type in (("GSM5621619_mimic365", "healthy_mimic365", "mimic-365_healthy", ""),
                                        ("GSM5621620_antimiR", "healthy_antimiR365", "antimiR-365_healthy", "antimiR_opposite_direction")):
    counts, features, barcodes = read_10x_mtx(RAW / f"{stem}_matrix.mtx.gz", RAW / f"{stem}_genes.tsv.gz", RAW / f"{stem}_barcodes.tsv.gz")
    libraries.append(Library(sample, "healthy", arm, counts, features, barcodes, control_type))

counts, features, barcodes = read_10x_mtx(RAW / "GSM5621621_LQT1_matrix.mtx.gz", RAW / "GSM5621621_LQT1_features.tsv.gz", RAW / "GSM5621621_LQT1_barcodes.tsv.gz")
calls = pd.read_csv(RAW / "GSM5621622_LQT1_HTOtags.csv.gz", dtype=str)
calls = calls[calls.num_features == "1"].set_index("cell_barcode").feature_call.reindex(barcodes)
# the record's 1:2:3 pooling of untreated, control and antimiR-365 matches the 668, 1,791 and 5,056 calls of tags 2, 3, 4
TAGS = {"HumanHashtag2": ("untreated_LQT1", "untreated"), "HumanHashtag3": ("antimiR-Ctrl_LQT1", "control_antimir"), "HumanHashtag4": ("antimiR-365_LQT1", "")}
wanted = np.flatnonzero(calls.isin(TAGS).to_numpy())
arm = calls.iloc[wanted].map(lambda value: TAGS[value][0]).tolist()
control_type = calls.iloc[wanted].map(lambda value: TAGS[value][1]).tolist()
libraries.append(Library("LQT1_lane", "LQT1", arm, counts[:, wanted], features, barcodes[wanted], control_type))
build_store("GSE185688", libraries, qc={"min_genes": 300, "max_genes": 6500, "max_mito": 0.20}, qc_source="authors, 10.1038/s41467-021-27856-7",
            source="authors' Cell Ranger matrices and hashtag calls, GEO GSE185688", out=store_path())
