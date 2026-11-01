"""Human breast milk cells from seven samples of gestational ages 27 to 39 weeks, 10x 3' v3 (BGISEQ-500 and NovaSeq 6000), and bulk exosomal miRNA-seq of the same milk samples. Reads only, paper not located.
No perturbation, so no comparison. Cell store from our STARsolo alignment with the EmptyDrops_CR cell filter (no thresholds from a paper or deposit); miRNA table per sample in tables/.
Usage: SRP690994_milk.py prepare | build <final dir>"""
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import pandas as pd
from _sample_arms import Library, build_store
from funmirbench.single_cell_store import read_10x_mtx
from tier2_plate import DATA

ROOT = DATA / "SRP690994"
ALIGNED = ROOT / "aligned"


def report() -> pd.DataFrame:
    return pd.read_csv(ROOT / "meta/SRP690994.tsv", sep="\t")


def prepare() -> None:
    ALIGNED.mkdir(exist_ok=True)
    table = report()
    cells = table[~table.library_name.str.startswith("miRNA")]
    (ALIGNED / "samples.tsv").write_text("".join(f"{r.library_name}\t{ROOT}/fastq/{r.run_accession}_1.fastq.gz\t{ROOT}/fastq/{r.run_accession}_2.fastq.gz\n" for r in cells.itertuples()))
    small = table[table.library_name.str.startswith("miRNA")]
    trimmed = small[small.library_layout == "SINGLE"]
    with_adapter = small[small.library_layout == "PAIRED"]
    (ALIGNED / "mirna_trimmed.manifest").write_text("".join(f"{ROOT}/fastq/{r.run_accession}.fastq.gz\t{r.library_name}\n" for r in trimmed.itertuples()))
    (ALIGNED / "mirna_adapter.manifest").write_text("".join(f"{ROOT}/fastq/{r.run_accession}_1.fastq.gz\t{r.library_name}\n" for r in with_adapter.itertuples()))


def build(final: pathlib.Path) -> None:
    (final / "store").mkdir(parents=True)
    (final / "tables").mkdir()
    table = report()
    libraries = []
    for row in table[~table.library_name.str.startswith("miRNA")].itertuples():
        folder = ALIGNED / row.library_name / "Solo.out/Gene/filtered"
        counts, features, barcodes = read_10x_mtx(folder / "matrix.mtx", folder / "features.tsv", folder / "barcodes.tsv")
        age = row.library_name.split("_")[0]
        libraries.append(Library(row.library_name, row.instrument_model.replace(" ", "_"), f"milk_{age}", counts, features, barcodes, unit=row.library_name, unit_type="biological"))
    build_store("SRP690994_milk", libraries, qc={}, qc_source="EmptyDrops_CR in STARsolo; no thresholds from a paper or deposit", source="our STARsolo 2.7.11b alignment against Ensembl 115, unique mappers, reads SRA SRP690994", out=final / "store/SRP690994_milk.h5ad")
    frames = [pd.read_csv(ALIGNED / name / "mirna_counts.tsv.gz", sep="\t", index_col=0) for name in ("mirna_trimmed", "mirna_adapter")]
    pd.concat(frames).fillna(0).astype("int64").T.rename_axis("mirna").to_csv(final / "tables/sample_mirna_counts.tsv.gz", sep="\t")
    pd.concat([pd.read_csv(ALIGNED / name / "library_stats.tsv", sep="\t") for name in ("mirna_trimmed", "mirna_adapter")]).to_csv(final / "tables/sample_mirna_library_stats.tsv", sep="\t", index=False)


if __name__ == "__main__":
    prepare() if sys.argv[1] == "prepare" else build(pathlib.Path(sys.argv[2]))
