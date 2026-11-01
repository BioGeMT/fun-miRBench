"""GSE225807 K562 CRISPRi RNA-binding-protein CROP-seq, one sample. Authors' mtx and barcode-to-guide map; label negative
is the control. One sample per arm, so three technical pseudo-replicates. Usage: python gse225807.py <extracted RAW folder>"""

import sys
import pathlib

import numpy as np
import pandas as pd

from funmirbench.single_cell_store import pseudo_replicate_units, read_10x_mtx
from screen_common import Targets, comparison_rows, plan_cells, write_screen

DATASET = "GSE225807_K562"
CONTROL = "negative"


def main(root: pathlib.Path) -> None:
    targets = Targets()
    counts, features, barcodes = read_10x_mtx(root / "GSM7056649_matrix.mtx.gz", root / "GSM7056649_features.tsv.gz",
                                              root / "GSM7056649_barcodes.tsv.gz")
    calls = pd.read_csv(root / "GSM7056650_bc_to_sgrna_mapping_3_06.csv.gz", sep="\t", index_col=0).set_index("barcode")
    bare = pd.Series(barcodes).str.replace(r"-\d+$", "", regex=True)
    table = calls.reindex(bare)
    rows = np.flatnonzero(table.RBP.notna())
    table = table.iloc[rows]
    control = (table.RBP == CONTROL).to_numpy()
    obs = pd.DataFrame({"sample": "GSM7056649", "batch": "GSM7056649", "replicate_unit": "GSM7056649",
                        "replicate_type": "technical", "arm": table.RBP.astype(str).to_numpy(), "is_control": control,
                        "control_type": np.where(control, "non_targeting", ""), "guide": table.sgRNA.astype(str).to_numpy(),
                        "n_guides": 1, "passed_qc": True})
    kept, positions = plan_cells(obs, targets)
    kept["replicate_unit"] = pseudo_replicate_units(kept["sample"], kept.arm, kept.passed_qc)
    caveats = ("CRISPRi knockdown; one sample, so replicate units are three technical pseudo-replicates: labels mean "
               "changed in this culture; few controls; barcodes as the authors' guide map keeps them")
    comparisons = comparison_rows(DATASET, kept.arm[~kept.is_control], targets, control_arm=CONTROL, caveats=caveats)
    write_screen(DATASET, counts.T.tocsr()[rows[positions]], kept, features[0].astype(str), comparisons=comparisons,
                 symbols=features[1].astype(str), source="authors' raw counts, GSE225807 GSM7056649 mtx")


if __name__ == "__main__":
    main(pathlib.Path(sys.argv[1]))
