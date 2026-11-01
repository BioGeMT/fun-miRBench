"""Holo-seq (HCC cells, one patient). Aligned raw counts from the source data of Velut et al. 2025, not the authors."""
import pathlib, sys, zipfile
sys.path.insert(0, str(pathlib.Path(__file__).parent))
import pandas as pd
from scipy import sparse
from cooccurrence_common import DATA, write_study

archive = zipfile.ZipFile(DATA / "CRA001133/deposit/velut_source_data.zip")
mirna = pd.read_csv(archive.open("Xiao_miRNA_raw_count.csv"), sep=r"\s+")
mrna = pd.read_csv(archive.open("Xiao_mRNA_raw_count.csv"), sep=r"\s+")
print(mirna.shape, mrna.shape, mirna.columns[:3].tolist(), mrna.columns[:3].tolist(), mirna.index[:3].tolist(), flush=True)
common = [c for c in mirna.columns if c in set(mrna.columns)]
cells = pd.DataFrame({"cell": common, "group": "HCC", "cell_type": "hepatocellular carcinoma"})
write_study("CRA001133_HCC", mirna=sparse.csr_matrix(mirna[common].to_numpy().T), mirna_features=mirna.index, mrna=sparse.csr_matrix(mrna[common].to_numpy().T), mrna_features=mrna.index, cells=cells, library_totals=None, library_total_source="none: the source data hold miRNA rows only, with no small RNA library total; the mRNA is a separate library", info={
    "feature_level": "mature_arm", "same_cell": "Holo-seq: small RNA and mRNA libraries from one cell (Xiao et al. 2018)",
    "values": "aligned raw counts, not the scImpute table", "grouping": "one tumour", "qc_source": "Velut et al. 2025 source data",
    "source": "third party: Velut et al. 2025 source data (Xiao_*_raw_count.csv); raw reads on GSA CRA001133, 14.2 GB, not aligned here (tier 2)",
    "caveats": "One patient. The original paper's supplement lists detected miRNAs only, with no counts."})
