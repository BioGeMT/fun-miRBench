"""GSE345058 A549 inducible Cas9 CROP-seq, 1,000 genes. Authors' raw mtx (cells x genes) and cell annotations. The four
authors' batches (eight lanes each) are the lane sets. Usage: python gse345058.py"""

import numpy as np
import pandas as pd

from screen_common import DATA, Targets, comparison_rows, plan_cells, read_mtx_rows, write_screen

DATASET = "GSE345058_A549"
CONTROL = "Nontargeting"


def main() -> None:
    targets = Targets()
    root = DATA / "GSE345058"
    cells = pd.read_csv(root / "GSE345058_cell_annotations.csv.gz")
    genes = pd.read_csv(root / "GSE345058_gene_annotations.csv.gz")
    control = (cells.Perturbation == CONTROL).to_numpy()
    obs = pd.DataFrame({"sample": cells.Lane.astype(str), "batch": cells.Batch.astype(str),
                        "replicate_unit": cells.Batch.astype(str), "replicate_type": "lane_set",
                        "arm": cells.Perturbation.astype(str), "is_control": control,
                        "control_type": np.where(control, "non_targeting", ""), "guide": cells.sgRNA.astype(str),
                        "n_guides": 1, "passed_qc": True})
    kept, rows = plan_cells(obs, targets)
    counts = read_mtx_rows(root / "GSE345058_raw_counts.mtx.gz", rows)
    caveats = "Cas9 knockout; replicate units are the authors' four batches of eight lanes; cells as the authors filtered them"
    comparisons = comparison_rows(DATASET, kept.arm[~kept.is_control], targets, control_arm=CONTROL, caveats=caveats)
    write_screen(DATASET, counts, kept, genes.ensembl_id.astype(str), comparisons=comparisons, symbols=genes.gene_symbol.astype(str),
                 source="authors' raw counts, GSE345058_raw_counts.mtx.gz, GEO GSE345058")


if __name__ == "__main__":
    main()
