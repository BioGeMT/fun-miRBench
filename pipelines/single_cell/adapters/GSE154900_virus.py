"""BC-1 primary effusion lymphoma, two 10x 3' v2 replicates (Landis et al. 2020), KSHV and EBV wild type. The GEO tables are integer but
SCTransform-style corrected (every cell of a replicate holds the same total), so they are not raw counts; raw counts come from our
STARsolo alignment of SRP273071 against Ensembl 115, which has no viral genes. Barcodes in the authors' tables are the deposit's cell filter.
No comparison: no sample-level arm exists; latent and lytic cells share one culture and arms would have to come from viral expression."""
import gzip
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from _sample_arms import DATA, Library, build_store, store_path
from funmirbench.single_cell_store import read_10x_mtx

ROOT = DATA / "GSE154900"
libraries = []
for replicate, deposit in (("rep1", "GSM4682311_bc1_replicate1"), ("rep2", "GSM4682312_bc1_replicate2")):
    with gzip.open(ROOT / "deposit" / f"{deposit}_normalized_counts.csv.gz", "rt") as handle:
        kept = {name.strip('"').rsplit(".", 1)[0] for name in handle.readline().rstrip("\n").split(",")[1:]}
    folder = ROOT / "aligned" / replicate / "Solo.out/Gene/raw"
    counts, features, barcodes = read_10x_mtx(folder / "matrix.mtx", folder / "features.tsv", folder / "barcodes.tsv")
    arm = ["BC-1" if barcode in kept else None for barcode in barcodes]
    found = sum(value is not None for value in arm)
    if found != len(kept):
        raise ValueError(f"{replicate}: {len(kept)} barcodes in the authors' table, {found} in the alignment")
    libraries.append(Library(replicate, "BC-1", arm, counts, features, barcodes, unit=replicate, unit_type="biological"))
build_store("GSE154900_virus", libraries, qc={}, qc_source="deposit's filter: barcodes in the authors' GEO tables; no thresholds applied here",
            source="our STARsolo 2.7.11b alignment against Ensembl 115, unique mappers, reads SRA SRP273071", out=store_path())
