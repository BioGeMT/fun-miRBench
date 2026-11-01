"""miR-192 mimic in lipid nanoparticles against saline, PBMC, hashed into one BD Rhapsody lane; one arm each, three
pseudo-replicates. Barcodes in the authors' hashtag table are the deposit's cell filter; no thresholds are stated. The
matrix holds fractional STARsolo counts (multi-mapping reads shared between genes): the store keeps them, and they are
rounded after summing."""

import numpy as np
import pandas as pd
from _sample_arms import DATA, Library, build_store, store_path
from funmirbench.single_cell_store import read_10x_mtx

RAW = DATA / "GSE292981" / "matrices"
TAGS = {"A955": ("saline", "saline"), "A965": ("miR-192_LNP_mimic", "")}

counts, features, barcodes = read_10x_mtx(RAW / "GSM8872121_WTA1_matrix.mtx.gz", RAW / "GSM8872121_WTA1_features.tsv.gz", RAW / "GSM8872121_WTA1_barcodes.tsv.gz")
fractional_share = counts.data[counts.data != np.round(counts.data)].sum() / counts.data.sum()
print(f"share of counts in fractional entries: {fractional_share:.4f}", flush=True)
tags = pd.read_csv(RAW / "GSM8872122_WTA1_Tag1_barcode.csv.gz", dtype=str)
tag = pd.Series(tags.Nreads_CellBC.str.split("_").str[0].to_numpy(), index=tags.Barcode).reindex(barcodes)
wanted = np.flatnonzero(tag.isin(TAGS).to_numpy())
arm = tag.iloc[wanted].map(lambda value: TAGS[value][0]).to_numpy()
control_type = tag.iloc[wanted].map(lambda value: TAGS[value][1]).to_numpy()
features = features[[0, 1]].reset_index(drop=True)
library = Library("WTA1", "lane", list(arm), counts[:, wanted], features, barcodes[wanted], list(control_type))
build_store("GSE292981", [library], qc={}, qc_source="deposit's filter: barcodes with a hashtag call in the authors' table; no thresholds stated in the paper",
            source="authors' STARsolo matrix, GEO GSE292981", out=store_path(), fractional=True)
