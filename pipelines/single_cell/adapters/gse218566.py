"""GSE218566 HEK293 PerturbSci-Kinetics, one pooled culture. Authors' whole-transcriptome counts and on-target cell metadata.
Replicate culture labels are not deposited, so the library date is the unit (technical sub-libraries). Usage: python gse218566.py"""

import numpy as np
import pandas as pd

from screen_common import DATA, Targets, comparison_rows, plan_cells, read_mtx_rows, write_screen

DATASET = "GSE218566_HEK293"
CONTROL = "NO-TARGET"


def main() -> None:
    targets = Targets()
    root = DATA / "GSE218566"
    barcodes = pd.read_csv(root / "GSM6752591_on_target_whole_tx.Barcodes.tsv.gz", header=None)[0]
    genes = pd.read_csv(root / "GSM6752591_on_target_whole_tx.Genes.tsv.gz", header=None)[0].astype(str)
    meta = pd.read_csv(root / "GSM6752591_on_target_cell_metadata.csv.gz").set_index("cell_names").reindex(barcodes)
    if meta.target_genes.isna().any():
        raise ValueError("cells without metadata")
    date = barcodes.str[:6]
    control = (meta.target_genes == CONTROL).to_numpy()
    obs = pd.DataFrame({"sample": barcodes.str.split(".").str[0].to_numpy(), "batch": date.to_numpy(),
                        "replicate_unit": date.to_numpy(), "replicate_type": "lane_set", "arm": meta.target_genes.to_numpy(),
                        "is_control": control, "control_type": np.where(control, "non_targeting", ""),
                        "guide": meta.target.astype(str).to_numpy(), "n_guides": 1, "passed_qc": True})
    kept, rows = plan_cells(obs, targets)
    counts = read_mtx_rows(root / "GSM6752591_on_target_whole_tx_count_matrix.mtx.gz", rows, cells_on_rows=False)
    caveats = ("CRISPRi knockdown; culture replicates are not labelled per cell, so replicate units are library dates, "
               "technical sub-libraries of one pooled culture; cells as the authors filtered them")
    comparisons = comparison_rows(DATASET, kept.arm[~kept.is_control], targets, control_arm=CONTROL, caveats=caveats)
    write_screen(DATASET, counts, kept, genes, comparisons=comparisons,
                 source="authors' whole-transcriptome counts, GSM6752591, GEO GSE218566")


if __name__ == "__main__":
    main()
