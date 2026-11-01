"""Tabula Sapiens total RNA-seq, donor TSP33: whole cells and nuclei from 22 tissues (figshare 31827436)."""
import pathlib, sys
sys.path.insert(0, str(pathlib.Path(__file__).parent))
import h5py, numpy as np, pandas as pd
from anndata.io import read_elem
from cooccurrence_common import DATA, write_study

with h5py.File(DATA / "tabula_sapiens_total/deposit/TSP33_total_RNA.h5ad", "r") as handle:
    counts = read_elem(handle["layers"]["raw_counts"])
    obs, var = read_elem(handle["obs"]), read_elem(handle["var"])
print(counts.shape, counts.nnz, var.biotype.value_counts().head(8).to_dict(), flush=True)
is_mirna = (var.biotype.astype(str) == "miRNA").to_numpy()
gene_ids = var.gene_ids.astype(str)
is_gene = ~is_mirna & gene_ids.str.startswith("ENSG").to_numpy()
print("miRNA", is_mirna.sum(), "ENSG genes", is_gene.sum(), flush=True)
cells = pd.DataFrame({"cell": obs.index.astype(str), "group": (obs.method.astype(str) + "_" + obs.tissue.astype(str) + "_" + obs.cell_type.astype(str)).to_numpy(), "method": obs.method.astype(str).to_numpy(), "tissue": obs.tissue.astype(str).to_numpy(), "cell_type": obs.cell_type.astype(str).to_numpy(), "batch": obs.batch.astype(str).to_numpy()})
write_study("TabulaSapiens_totalRNA", mirna=counts[:, np.flatnonzero(is_mirna)], mirna_features=var.index.astype(str)[is_mirna], mrna=counts[:, np.flatnonzero(is_gene)], mrna_features=gene_ids[is_gene], cells=cells, library_totals=pd.Series(np.asarray(counts.sum(axis=1)).ravel(), index=cells.cell.to_numpy()), library_total_source="sum of every feature in the authors' raw_counts layer for the cell", info={
    "feature_level": "precursor", "same_cell": "one total RNA library per cell or nucleus (figshare record)", "values": "the authors' raw_counts layer, partly fractional (multi-mappers)",
    "grouping": "method (cells or nuclei), tissue and the authors' cell type", "qc_source": "deposit (cells in the authors' h5ad)", "source": "authors, figshare 31827436",
    "caveats": "One donor. Features are miRBase precursors (hsa-mir-*). Nuclei hold nascent rather than mature cytoplasmic miRNA. Raw reads need a data transfer agreement."})
