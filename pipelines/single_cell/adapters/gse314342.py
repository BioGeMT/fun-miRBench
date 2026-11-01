"""GSE314342 genome-scale CRISPRi in primary CD4 T cells (authors' per donor-condition h5ad, S3). Target and non-targeting cells were
extracted per file; authors' low_quality flag is the cell QC; multi-guide cells are out. Donors are the replicates and the batches.
Usage: python gse314342.py <Rest|Stim8hr|Stim48hr>"""

import sys

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from screen_common import DATA, Targets, comparison_rows, plan_cells, write_screen

CONTROL = "non-targeting"
DONORS = ("D1", "D2", "D3", "D4")


def main(condition: str) -> None:
    dataset = f"GSE314342_CD4_{condition}"
    targets = Targets()
    parts, tables, var = [], [], None
    for donor in DONORS:
        path = DATA / "GSE314342" / "extracted" / f"{donor}_{condition}.h5ad"
        piece = ad.read_h5ad(path)
        var = piece.var
        meta = piece.obs
        control = (meta.is_nt == "True").to_numpy()
        tables.append(pd.DataFrame({"sample": meta.lane_id.astype(str).to_numpy(), "batch": donor, "replicate_unit": donor,
                                    "replicate_type": "biological", "arm": np.where(control, CONTROL, meta.perturbed_gene_name.astype(str)),
                                    "is_control": control, "control_type": np.where(control, "non_targeting", ""),
                                    "guide": meta.guide_id.astype(str).to_numpy(), "n_guides": 1,
                                    "passed_qc": (meta.low_quality != "True").to_numpy()}))
        parts.append(piece.X.tocsr())
    obs = pd.concat(tables, ignore_index=True)
    kept, rows = plan_cells(obs, targets)
    caveats = (f"CRISPRi knockdown in primary CD4 T cells, {condition}; replicate units are the donors; authors' low-quality flag; "
               "multi-guide cells left out; targets with few cells stay listed")
    comparisons = comparison_rows(dataset, kept.arm[(~kept.is_control) & kept.passed_qc], targets, control_arm=CONTROL, caveats=caveats)
    write_screen(dataset, sparse.vstack(parts).tocsr()[rows].astype(np.int32), kept, var.gene_ids.astype(str).reset_index(drop=True), comparisons=comparisons,
                 symbols=var.gene_name.astype(str).reset_index(drop=True), source="authors' h5ad per donor and condition, s3://genome-scale-tcell-perturb-seq/marson2025_data")


if __name__ == "__main__":
    main(sys.argv[1])
