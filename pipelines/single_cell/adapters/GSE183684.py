"""CRISPR deletion of the enhancer of the miR-9-2 primary transcript against isogenic control, human retinal organoids,
weeks 5, 12, 20, 28. Replicates are independent clones per the paper. QC: deposit's filter (authors' barcode lists)."""

import re
import pathlib

from _sample_arms import DATA, Library, build_store, store_path
from funmirbench.single_cell_store import read_10x_mtx

RAW = DATA / "GSE183684" / "raw"
libraries = []
for path in sorted(RAW.glob("GSM*_barcodes.tsv.gz")):
    stem = str(path).removesuffix("_barcodes.tsv.gz")
    genotype, week, clone = re.search(r"_(wt|ko)(\d+)-(\d+)$", stem).groups()
    counts, features, barcodes = read_10x_mtx(pathlib.Path(f"{stem}_matrix.mtx.gz"), pathlib.Path(f"{stem}_features.tsv.gz"), pathlib.Path(f"{stem}_barcodes.tsv.gz"))
    arm = f"{genotype}_wk{week}"
    libraries.append(Library(f"{genotype}{week}-{clone}", f"wk{week}", arm, counts, features, barcodes,
                             "isogenic_control" if genotype == "wt" else "", unit=f"{genotype}{week}-{clone}"))
build_store("GSE183684", libraries, qc={}, qc_source="deposit's filter: barcode lists of the authors' matrices; thresholds not read (paper full text not available)",
            source="authors' snRNA-seq matrices, GEO GSE183684", out=store_path())
