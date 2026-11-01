"""PSCSR-seq V2 (A549, HeLa, K562, HEK293T mix; two chips). Deposit: log1p per 10,000 for miRNA and for mRNA."""
import pathlib, sys
sys.path.insert(0, str(pathlib.Path(__file__).parent))
import numpy as np, pandas as pd
from scipy import sparse
from cooccurrence_common import DATA, write_study

root = DATA / "GSE226714"
mirna_parts, mrna_parts, cell_parts, features = [], [], [], {}
tables = {}
for chip, mirna_file, mrna_file in (("chip1", "GSM7082512_cells_chip1_miRNA.txt.gz", "GSM7082513_cells_chip1_mRNA.txt.gz"), ("chip2", "GSM7082514_cells_chip2_miRNA.txt.gz", "GSM7082515_cells_chip2_mRNA.txt.gz")):
    mirna = pd.read_csv(root / "extracted" / mirna_file, sep="\t", index_col=0)
    mrna = pd.read_csv(root / "extracted" / mrna_file, sep="\t", index_col=0)
    barcodes = pd.read_csv(root / "deposit" / f"GSE226714_cells_{chip}_barcodes.txt.gz", sep="\t", header=None, names=["barcode", "label"])
    common = barcodes[barcodes.barcode.isin(mirna.columns) & barcodes.barcode.isin(mrna.columns)]
    y = np.expm1(mirna[common.barcode].to_numpy())
    smallest = np.where(y > 0, y, np.inf).min(axis=0)  # one count per cell: every value is an integer multiple of it
    counts = np.round(np.divide(y, smallest, where=np.isfinite(smallest), out=np.zeros_like(y)))
    if np.abs(counts * np.where(np.isfinite(smallest), smallest, 0) - y).max() > 1e-6:
        raise ValueError("miRNA values are not integer multiples of the smallest value of their cell")
    tables[chip] = (pd.DataFrame(counts, index=mirna.index, columns=[f"{chip}_{b}" for b in common.barcode]),
                    pd.DataFrame(np.expm1(mrna[common.barcode].to_numpy()), index=mrna.index, columns=[f"{chip}_{b}" for b in common.barcode]))
    line = common.label.str.replace("-MS", "", regex=False)
    cell_parts.append(pd.DataFrame({"cell": [f"{chip}_{b}" for b in common.barcode], "group": line + "_" + chip, "cell_line": line, "chip": chip}))
mirna = pd.concat([tables[c][0] for c in tables], axis=1).fillna(0)
mrna = pd.concat([tables[c][1] for c in tables], axis=1).fillna(0)
cells = pd.concat(cell_parts, ignore_index=True)
write_study("GSE226714_lines", mirna=sparse.csr_matrix(mirna[cells.cell].to_numpy().T), mirna_features=mirna.index, mrna=sparse.csr_matrix(mrna[cells.cell].to_numpy().T), mrna_features=mrna.index, cells=cells, library_totals=None, library_total_source="none: the deposit holds miRNA rows only, with no small RNA library total; the mRNA is a separate library", info={
    "feature_level": "mature_arm", "same_cell": "parallel small RNA and mRNA libraries from each nanowell (GSM7082512-15, barcodes shared)",
    "values": "miRNA: integer counts recovered exactly from the log1p per 10,000 values (each cell's values are integer multiples of its smallest; assumes that smallest is one count). mRNA: per 10,000 from the same transform, fractional (ratios such as 1/7 and 1/5, so multi-mapper shares), no integer counts available",
    "grouping": "cell line within chip (the deposit's barcode files)", "qc_source": "deposit (authors' barcode files)", "source": "authors, GSE226714 supplementary files",
    "caveats": "Mixed culture; mouse lung chips left out. Features of a chip missing on the other are zero there, which has no effect because groups are within chip."})
