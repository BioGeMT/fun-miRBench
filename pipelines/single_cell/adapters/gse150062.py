"""GSE150062 CRISPRi iPSC, day 8 of neural induction. Raw RNA counts and authors' singlet guide calls from the deposited Seurat
object (extracted with seurat_counts.R; the deposited mtx is SCT-corrected, not raw). The six 10x samples are the units.
lncRNA library IDs map to hosts as in experiments.tsv. Usage: python gse150062.py"""

import numpy as np
import pandas as pd
from scipy import io, sparse

from funmirbench.single_cell_pipeline import symbol_ensembl_ids, symbol_map_from_gtf
from screen_common import DATA, GTF, Targets, comparison_rows, plan_cells, write_screen

DATASET = "GSE150062_iPSC"
CONTROL = "Non-Targeting"
LIBRARY_IDS = {"LH12331": "MIR302CHG", "LH14698": "ENSG00000270933", "LH14970": "LINC-PINT", "LH15346": "PVT1", "LH15347": "PVT1",
               "LH07082": "PLEKHM1P1", "LH08564": "MIR4754"}


def main() -> None:
    targets = Targets()
    root = DATA / "GSE150062" / "x_rna"
    meta = pd.read_csv(root / "meta.tsv", sep="\t", dtype=str, keep_default_na=False)
    features = (root / "features.tsv").read_text().splitlines()
    arm = meta.guide_target.replace(LIBRARY_IDS)
    arm = arm.map(lambda value: targets.symbol_of(value) or value)
    control = (arm == CONTROL).to_numpy()
    obs = pd.DataFrame({"sample": meta["sample"], "batch": meta["sample"], "replicate_unit": meta["sample"], "replicate_type": "lane_set",
                        "arm": arm, "is_control": control, "control_type": np.where(control, "non_targeting", ""),
                        "guide": meta.guide_identity, "n_guides": 1, "passed_qc": True})
    kept, rows = plan_cells(obs, targets)
    counts = sparse.csr_matrix(io.mmread(root / "counts.mtx").T)[rows].astype(np.int32)
    caveats = ("CRISPRi knockdown in day 8 neural induction, one differentiating culture; replicate units are the six 10x samples; "
               "high MOI, singlet calls as the authors made them; only day 8 is in the deposited Seurat object; PVT1 pools LH15346 and LH15347")
    comparisons = comparison_rows(DATASET, kept.arm[~kept.is_control], targets, control_arm=CONTROL, caveats=caveats)
    # most of the 60,497 features are GENCODE 19 clones and lncRNA library IDs: only symbols the GTF knows are mapped
    names = pd.Series(features, dtype=str)
    symbol_map = symbol_map_from_gtf(GTF)
    known = names[names.isin(symbol_map.index if hasattr(symbol_map, "index") else list(symbol_map))]
    mapped = symbol_ensembl_ids(known, symbol_map)[0].reindex(names.index, fill_value="")
    write_screen(DATASET, counts, kept, names, comparisons=comparisons, gene_ids=mapped,
                 source="authors' raw RNA counts and singlet calls, GSE150062_seurat.rds.gz, GEO GSE150062")


if __name__ == "__main__":
    main()
