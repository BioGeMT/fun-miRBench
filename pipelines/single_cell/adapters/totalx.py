"""TotalX droplet total RNA-seq: figshare 29803628 and the GEO HEK293T arms (GSE315939). Argument: HEK293T, Huh7, PBMC or brain."""
import pathlib, sys, tarfile
sys.path.insert(0, str(pathlib.Path(__file__).parent))
import anndata as ad, h5py, numpy as np, pandas as pd
from anndata.io import read_elem
from cooccurrence_common import DATA, SCRATCH, write_study

deposit = DATA / "totalx/deposit"
SOURCE = "authors, figshare 29803628 (GEO GSE315939 for the HEK293T arms)"
CAVEATS = "Features are genes of biotype miRNA in the authors' annotation, so reads are assigned at the locus, not to an arm. Droplet data: a cell without miRNA counts is left out by the step."


def load(path, layer=None):
    """Only obs, var and the count matrix: the whole object can hold several copies of the matrix."""
    with h5py.File(path, "r") as handle:
        return ad.AnnData(X=read_elem(handle["layers"][layer] if layer else handle["X"]), obs=read_elem(handle["obs"]), var=read_elem(handle["var"]))


arm = sys.argv[1]
if arm == "HEK293T":
    extracted = SCRATCH / "totalx_hek"
    extracted.mkdir(exist_ok=True)
    with tarfile.open(deposit / "GSE315939_RAW.tar") as archive:
        archive.extractall(extracted, filter="data")
    tables = []
    for name, path in (("TotalX", extracted / "GSM9440658_SSTX_counts_adata.h5ad"), ("miRNAplus_noDASH", deposit / "TotalX_HEK293T.h5ad"), ("miRNAplus_2xDASH", extracted / "GSM9440660_SSTX_counts_adata_2XDASH.h5ad")):
        table = load(path)
        print(name, table.shape, table.var.columns.tolist(), flush=True)
        table.var_names = table.var.gene_id.astype(str).to_numpy()
        table.obs = pd.DataFrame({"group": f"HEK293T_{name}"}, index=[f"{name}_{b}" for b in table.obs_names])
        table.var = table.var[["gene_id", "biotype"]].astype(str)
        table.obs["library_total"] = np.asarray(table.X.sum(axis=1)).ravel()
        tables.append(table)
    table = ad.concat(tables, join="inner", merge="same")
    grouping = "GEO arm"
elif arm == "Huh7":
    table = load(deposit / "TotalX_Huh7_DENV2_raw_counts_scrublet.h5ad")
    table.obs["group"] = table.obs.source.astype(str)
    grouping = "infection state (source)"
elif arm == "PBMC":
    table = load(deposit / "TotalX_PBMC_final_3conditions.h5ad", "counts")
    table.obs["group"] = table.obs.source.astype(str) + "_" + table.obs.cell_type_lvl2.astype(str)
    grouping = "protocol arm (source) and the authors' cell type (cell_type_lvl2)"
    CAVEATS += " The three source values are protocol arms: SSTX, SSTX-pA+ and 10X-dUTSO; not all capture non-polyadenylated RNA."
elif arm == "brain":
    table = load(deposit / "TotalX_Brain_log1p_scvi_full_annotated.h5ad", "counts")
    table.obs["group"] = table.obs.brain.astype(str) + "_" + table.obs.cell_type.astype(str)
    grouping = "donor and the authors' cell type"
else:
    raise SystemExit(f"unknown arm {arm}")
library_totals = table.obs.library_total if arm == "HEK293T" else pd.Series(np.asarray(table.X.sum(axis=1)).ravel(), index=table.obs_names)
gene_ids = table.var.gene_id.astype(str).to_numpy()
is_mirna = (table.var.biotype.astype(str) == "miRNA").to_numpy()
is_gene = ~is_mirna & pd.Series(gene_ids).str.startswith("ENSG").to_numpy()
cells = pd.DataFrame({"cell": table.obs_names, "group": table.obs.group.astype(str).to_numpy()})
write_study(f"TotalX_{arm}", mirna=table.X[:, np.flatnonzero(is_mirna)], mirna_features=gene_ids[is_mirna], mrna=table.X[:, np.flatnonzero(is_gene)], mrna_features=gene_ids[is_gene], cells=cells, library_totals=library_totals, library_total_source="sum of every feature in the authors' count matrix for the cell", info={
    "feature_level": "precursor_locus", "same_cell": "one droplet library per cell (TotalX)", "values": "the authors' raw counts", "grouping": grouping,
    "qc_source": "deposit (cells in the authors' h5ad, after their doublet filter)", "source": SOURCE, "caveats": CAVEATS})
