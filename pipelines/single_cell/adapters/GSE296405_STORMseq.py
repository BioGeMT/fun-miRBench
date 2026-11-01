"""STORM-seq and VASA-seq plates (K-562, HEK293T, RMG-2): one total RNA library per cell, kb counts on Ensembl 101."""
import pathlib, sys
sys.path.insert(0, str(pathlib.Path(__file__).parent))
import pandas as pd
from scipy import sparse
from cooccurrence_common import DATA, mirna_ensembl_ids, write_study

root = DATA / "GSE296405"
storm = pd.read_csv(root / "deposit/GSE296405_storm_hek293t_k562_rmg2_raw_umi_counts_with_multimappers_kb_ens101.txt.gz", sep="\t", index_col=0)
parts = {"STORM": storm}
for plate, path in (("VASAplate1", next((root / "extracted").glob("GSM8969810*"))), ("VASAplate2", next((root / "extracted").glob("GSM8969811*")))):
    parts[plate] = pd.read_csv(path, sep="\t", index_col=0)
cell_rows = []
for assay, table in parts.items():
    for column in table.columns:
        cell_rows.append((f"{assay}_{column}", column.rsplit("_", 1)[1], assay))
    table.columns = [f"{assay}_{c}" for c in table.columns]
table = pd.concat(parts.values(), axis=1).fillna(0)
cells = pd.DataFrame(cell_rows, columns=["cell", "cell_line", "assay"])
print(cells.groupby(["assay", "cell_line"]).size().to_dict(), flush=True)
cells = cells[cells.cell_line.isin(["K562", "HEK293T", "RMG2"])].reset_index(drop=True)  # the rest are the plate's controls
cells["group"] = cells.cell_line + "_" + cells.assay
table.index = table.index.str.replace(r"\.\d+$", "", regex=True)
table = table.groupby(level=0).sum()[cells.cell]
is_mirna = table.index.isin(mirna_ensembl_ids())
write_study("GSE296405_STORMseq", mirna=sparse.csr_matrix(table[is_mirna].to_numpy().T), mirna_features=table.index[is_mirna], mrna=sparse.csr_matrix(table[~is_mirna].to_numpy().T), mrna_features=table.index[~is_mirna], cells=cells, library_totals=table.sum(axis=0), library_total_source="sum of every row of the deposited table of the cell's library", info={
    "feature_level": "precursor_locus", "same_cell": "one total RNA library per cell (GEO design)", "values": "kb UMI counts with multimappers, fractional",
    "grouping": "cell line and assay", "qc_source": "deposit (cells in the authors' matrices)", "source": "authors, GSE296405 supplementary files",
    "caveats": "Random-primed total RNA with ribosomal reduction, so reads at miRNA loci are mostly precursor or primary transcript; the authors do not claim mature miRNA quantification. miRNA features come from the Ensembl 115 biotype of the Ensembl 101 IDs. Smart-seq3xpress arm left out: poly(A) priming, no small RNA."})
