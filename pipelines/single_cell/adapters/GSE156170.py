"""miR-124 suicide circuit in H1 teratomas grown in mice: ganciclovir kills cells without miR-124 activity. Human cells as
called by the authors' species table. The two injection routes are arms; each arm's cells are split into three
pseudo-replicates, since the control is one teratoma. QC: deposit's filter."""

import numpy as np
import pandas as pd

from _sample_arms import DATA, Library, build_store, store_path
from funmirbench.single_cell_store import read_10x_mtx

RAW = DATA / "GSE156170" / "extracted" / "miR-124"
SAMPLES = [("dm-ter-miRNA-neg", "no_GCV", "no_ganciclovir"), ("dm-ter-miRNA-pos-it1", "GCV_IT", ""), ("dm-ter-miRNA-pos-it2", "GCV_IT", ""),
           ("dm-ter-miRNA-pos-ipit1", "GCV_IPIT", ""), ("dm-ter-miRNA-pos-ipit2", "GCV_IPIT", "")]
libraries = []
for name, arm, control_type in SAMPLES:
    folder = RAW / name
    counts, features, barcodes = read_10x_mtx(folder / "matrix.mtx.gz", folder / "features.tsv.gz", folder / "barcodes.tsv.gz")
    human_genes = features[0].str.startswith("hg19_").to_numpy()
    counts, features = counts[np.flatnonzero(human_genes), :], features[human_genes].reset_index(drop=True)
    features = pd.DataFrame({0: features[0].str.removeprefix("hg19_"), 1: features[1].str.removeprefix("hg19_")})
    calls = pd.read_csv(RAW / f"{name}_species_class.csv").set_index("barcode").call.reindex(barcodes)
    human = np.flatnonzero((calls == "hg19").to_numpy())
    libraries.append(Library(name, "series", arm, counts[:, human], features, barcodes[human], control_type, pseudo_group=arm))
build_store("GSE156170", libraries, qc={}, qc_source="deposit's filter and authors' species call", source="authors' matrices, GEO GSE156170", out=store_path())
