"""Doxycycline-inducible miR-199a precursor with GFP on one episome, HEK293 lines i199 and i199-KTN1, nine doses pooled per 10x
library. GFP is in every cell (leaky vector), so arms are defined after the fact: tertiles of eGFP UMIs per 10,000 gene UMIs
within each library, high against low; the middle tertile is stored but not compared. Three pseudo-replicates per arm.
Counts: our STARsolo alignment (human_benchmark/source/align_srp150046.sh), filtered cells; no thresholds applied."""

import pandas as pd
from _sample_arms import DATA, Library, build_store, store_path
from funmirbench.single_cell_store import read_10x_mtx

ROOT = DATA / "SRP150046"
libraries = []
for library, line in (("i199", "i199"), ("i199_ktn1", "i199-KTN1")):
    folder = ROOT / "starsolo" / f"{library}_genes" / "Solo.out" / "Gene" / "filtered"
    counts, features, barcodes = read_10x_mtx(folder / "matrix.mtx", folder / "features.tsv", folder / "barcodes.tsv")
    gfp = pd.read_csv(ROOT / "arms" / f"{library}_gfp_arms.csv").set_index("barcode").gfp_per_10k.reindex(barcodes)
    tertile = pd.qcut(gfp.rank(method="first"), 3, labels=["low", "mid", "high"]).astype(str)
    arm = [f"GFP_{value}_{line}" for value in tertile]
    control_type = ["low_GFP_same_library" if value == "low" else "" for value in tertile]
    libraries.append(Library(line, line, arm, counts, features, barcodes, control_type))
build_store("SRP150046", libraries, qc={}, qc_source="STARsolo Cell Ranger-like cell filter, our alignment; no further thresholds", source="our STARsolo 2.7.11b alignment against Ensembl 115, reads SRA SRP150046", out=store_path())
