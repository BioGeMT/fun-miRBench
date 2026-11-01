"""DMF-DM-seq (Chen et al. 2024, Anal Chem): mRNA and miRNA libraries of the same single cell, MCF10A, MCF7, MDA-MB-231, K562.
Deposit: reads only (SRP500050). `prepare` writes the manifests; tier2_SRP500050_DMFDMseq.sh aligns; `build` writes the study.

The ten mRNA samples named MDA-MB-21-mRNA-N are MDA-MB-231: each BioSample says "mRNA data of the Nth single MDA-MB-231
cell" and sample_type "MDA-MB-231 cell culture-N", the same as its miRNA partner. Cells pair by cell line and number."""
import pathlib
import re
import sys

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import numpy as np
import pandas as pd
from cooccurrence_common import DATA, write_study
from funmirbench.single_cell_store import read_10x_mtx

ROOT = DATA / "SRP500050"
ALIGNED = ROOT / "aligned"
NAME = re.compile(r"^(MCF10A|MCF7|K562|MDA-MB-21|MDA-MB-231)-(mRNA|miRNA)-(\d+)$")


def runs() -> pd.DataFrame:
    report = pd.read_csv(ROOT / "meta/SRP500050.tsv", sep="\t")
    parts = report.library_name.str.extract(NAME)
    if parts.isna().any(axis=None):
        raise ValueError("library names outside the known pattern")
    report["line"] = parts[0].replace({"MDA-MB-21": "MDA-MB-231"})
    report["kind"] = parts[1]
    report["cell"] = report.line + "_" + parts[2]
    return report


def prepare() -> None:
    ALIGNED.mkdir(exist_ok=True)
    report = runs()
    fastq = ROOT / "fastq"
    mrna, mirna = report[report.kind == "mRNA"], report[report.kind == "miRNA"]
    (ALIGNED / "mrna.manifest").write_text("".join(f"{fastq}/{r.run_accession}_1.fastq.gz\t{fastq}/{r.run_accession}_2.fastq.gz\t{r.cell}\n" for r in mrna.itertuples()))
    (ALIGNED / "mirna.manifest").write_text("".join(f"{fastq}/{r.run_accession}_1.fastq.gz\t{r.cell}\n" for r in mirna.itertuples()))


def build() -> None:
    report = runs()
    both = report.pivot(index="cell", columns="kind", values="run_accession")
    if both.isna().any(axis=None):
        raise ValueError(f"cells without both libraries: {both[both.isna().any(axis=1)].index.tolist()}")
    lines = report.drop_duplicates("cell").set_index("cell").line
    cells = pd.DataFrame({"cell": both.index, "group": lines[both.index].to_numpy(), "cell_line": lines[both.index].to_numpy(),
                          "run_mrna": both.mRNA.to_numpy(), "run_mirna": both.miRNA.to_numpy()})
    solo = ALIGNED / "mrna/Solo.out/Gene/raw"
    mrna, features, barcodes = read_10x_mtx(solo / "matrix.mtx", solo / "features.tsv", solo / "barcodes.tsv")
    order = pd.Index(barcodes).get_indexer(cells.cell)
    if (order < 0).any():
        raise ValueError("cells missing from the mRNA alignment")
    mirna = pd.read_csv(ALIGNED / "mirna/mirna_counts.tsv.gz", sep="\t", index_col=0).loc[cells.cell]
    stats = pd.read_csv(ALIGNED / "mirna/library_stats.tsv", sep="\t", index_col=0).loc[cells.cell]
    write_study("SRP500050_DMFDMseq", mirna=mirna.to_numpy(np.int64), mirna_features=mirna.columns, mrna=mrna.T.tocsr()[order, :], mrna_features=features[0],
                cells=cells, library_totals=stats.inserts, library_total_source="reads of the small RNA library with the adapter and an insert of 18-30 nt (our count from the reads; the deposit has no table)", info={
        "feature_level": "mature_arm", "same_cell": "yes: the two libraries are separated from one lysate on the chip (BioProject description; paper paywalled, methods not read)",
        "values": "integers from our alignment. mRNA: STARsolo SmartSeq, Ensembl 115, reads in the exons of one gene, unique mappers, no deduplication. miRNA: read 1, TGGAATTCTCGGGTGCCAAGG trimmed, inserts of 18-30 nt, bowtie 1 to miRBase 22.1 hairpins, reads on one arm only",
        "grouping": "cell line", "qc_source": "none: paper paywalled, the deposit has no filter; all 35 pairs kept", "source": "our alignment of SRP500050 reads",
        "caveats": "mRNA samples are named MDA-MB-21 but their BioSamples say MDA-MB-231, so all 35 cells pair. K562 has 5 cells, below the step's 10. Few cells per line, not replicated across cultures as far as the record shows."})


if __name__ == "__main__":
    {"prepare": prepare, "build": build}[sys.argv[1]]()
