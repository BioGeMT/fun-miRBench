"""agoTRIBE, Smart-seq3 (Sethi et al. 2024): HEK-293T transfected with the Argonaute2-ADAR2 fusion (540 cells) and control (26), K562 undifferentiated (41) and differentiated 24 h (36),
all transfected. Only the expression side is built here: editing events are not called. ENA holds one FASTQ per cell (read 2 by its file name), so reads are single-end; only cells the authors kept are in the deposit (paper: 566 HEK-293T cells), so
no further QC. Reads are counted, not UMIs: Smart-seq3 UMIs sit in 5' reads of read 1, which the deposit does not hold. Bulk libraries of the same project are left out.
Usage: PRJNA994505.py prepare | build HEK293T|K562 <store path>"""
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import pandas as pd
from tier2_plate import DATA, write_plate_store

ROOT = DATA / "PRJNA994505"
ALIGNED = ROOT / "aligned"


def cells_of_report() -> pd.DataFrame:
    table = pd.read_csv(ROOT / "meta/PRJNA994505.tsv", sep="\t")
    table = table[table.library_source.str.contains("SINGLE")].copy()
    table["kind"] = table.library_name.str.replace(r"^[ACGT+-]+_", "", regex=True)
    table["pool"] = table.run_alias.str.extract(r"^(SS3_21_\d+)")[0].fillna("K562_pool")
    table["cell"] = table.run_accession
    if len(table) != 643 or table.fastq_ftp.str.contains(";").any():
        raise ValueError("the run report changed")
    return table.set_index("cell")


def prepare() -> None:
    ALIGNED.mkdir(exist_ok=True)
    table = cells_of_report()
    runs = ROOT / "meta/single_cell_runs.txt"
    runs.write_text("\n".join(table.run_accession) + "\n")
    (ALIGNED / "cells.manifest").write_text("".join(f"{ROOT}/fastq/{r.run_accession}.fastq.gz\t-\t{cell}\n" for cell, r in table.iterrows()))


def build(line: str, out: pathlib.Path) -> None:
    table = cells_of_report()
    if line == "HEK293T":
        table = table[table.kind.isin(["agoTRIBE", "control"])]
        arm = table.kind.map({"agoTRIBE": "agoTRIBE", "control": "control"})
        control = table.kind.map({"agoTRIBE": "", "control": "untransfected"})
    else:
        table = table[table.kind.str.startswith("K-562")]
        arm = table.kind.map({"K-562_undifferentiated": "K562_undifferentiated", "K-562_differentiated_24hr": "K562_differentiated_24h"})
        control = table.kind.map({"K-562_undifferentiated": "undifferentiated", "K-562_differentiated_24hr": ""})
    cells = pd.DataFrame({"sample": table.pool + "_" + arm, "batch": table.pool, "arm": arm, "control_type": control, "passed_qc": True}, index=table.index)
    write_plate_store(f"PRJNA994505_{line}", ALIGNED / "cells/Solo.out/Gene/raw", cells, out,
                      "our STARsolo SmartSeq alignment of SRA reads against Ensembl 115, unique mappers, reads not UMIs; cell QC: the authors' filter, only kept cells are deposited")


if __name__ == "__main__":
    prepare() if sys.argv[1] == "prepare" else build(sys.argv[2], pathlib.Path(sys.argv[3]))
