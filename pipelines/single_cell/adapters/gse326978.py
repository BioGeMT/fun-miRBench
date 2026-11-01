"""GSE326978 Cas13 array screen in MDA-MB-231 and HAP1 (Cell Genomics 2026). Authors' per-lane mtx and per-cell guide UMI
tables (g1_g2_umi); a cell belongs to the target of its one array; NT_NT arrays are the controls. Each 10x lane is a batch and a
replicate unit. Usage: python gse326978.py <MDAMB231|HAP1>"""

import sys

import numpy as np
import pandas as pd
from scipy import sparse

from funmirbench.single_cell_store import read_10x_mtx
from screen_common import DATA, Targets, comparison_rows, plan_cells, write_screen

CONTROL = "NT_NT"
TABLES = {"MDAMB231": "GSE326978_mda.g1_g2_umi.csv.gz", "HAP1": "GSE326978_hap.g1_g2_umi.csv.gz"}


def main(line: str) -> None:
    dataset = f"GSE326978_{line}"
    targets = Targets()
    root = DATA / "GSE326978"
    calls = pd.read_csv(root / TABLES[line], index_col=0)
    reference = pd.read_csv(root / "GSE326978_feature_reference.csv.gz").set_index("id")
    values = calls.to_numpy()
    n_called = (values > 0).sum(axis=1)
    best = calls.columns.to_numpy()[values.argmax(axis=1)]
    name = reference.name.reindex(best).astype(str).to_numpy()
    # target_gene_name in the reference is a placeholder (CD55 for all); XLOC names map to Ensembl through Table S1A of the paper
    xloc = pd.read_csv(root / "tables" / "table_S1A_targets.tsv", sep="\t", dtype=str).set_index("lncRNA").ENSEMBL_ID
    gene = pd.Series(name).map(lambda value: xloc.get(value, value)).astype(str).str.split(".").str[0].to_numpy()
    control = name == CONTROL
    arm = np.where(control, CONTROL, gene)
    print(pd.Series(arm).map(targets.symbol_of).value_counts().to_dict(), n_called.max(), (n_called == 1).mean(), flush=True)
    cells = pd.Series(calls.index, dtype=str)
    lane, bare = cells.str.split("_", n=1).str[0].to_numpy(), cells.str.split("_", n=1).str[1].to_numpy()
    obs = pd.DataFrame({"sample": lane, "batch": lane, "replicate_unit": lane, "replicate_type": "lane_set", "arm": arm,
                        "is_control": control, "control_type": np.where(control, "non_targeting", ""), "guide": best,
                        "n_guides": n_called, "passed_qc": n_called == 1})
    kept, rows = plan_cells(obs, targets)
    chosen_lane, chosen_bare = lane[rows], bare[rows]
    parts, order, features = [], [], None
    for code in pd.unique(chosen_lane):
        prefix = next((root / "x").glob(f"*_{code}_GEX_matrix.mtx.gz")).name.replace("_matrix.mtx.gz", "")
        counts, features, barcodes = read_10x_mtx(root / "x" / f"{prefix}_matrix.mtx.gz", root / "x" / f"{prefix}_features.tsv.gz",
                                                  root / "x" / f"{prefix}_barcodes.tsv.gz")
        index = pd.Series(np.arange(len(barcodes)), index=barcodes)
        where = np.flatnonzero(chosen_lane == code)
        parts.append(counts.T.tocsr()[index.loc[chosen_bare[where]].to_numpy()])
        order.append(where)
    matrix = sparse.vstack(parts).tocsr()[np.argsort(np.concatenate(order))]
    caveats = ("Cas13 RNA knockdown; a related paper from this line of work, Liang et al., Cell 2024, was retracted, while this "
               "series belongs to the later Cell Genomics 2026 paper; Cas13 cuts spliced and nuclear RNA, so the pri-miRNA may or "
               "may not fall; replicate units are the 10x lanes; few cells per target")
    comparisons = comparison_rows(dataset, kept.arm[(~kept.is_control) & kept.passed_qc], targets, control_arm=CONTROL, caveats=caveats)
    write_screen(dataset, matrix, kept, features[0].astype(str), comparisons=comparisons, symbols=features[1].astype(str),
                 source=f"authors' raw counts and {TABLES[line]}, GEO GSE326978")


if __name__ == "__main__":
    main(sys.argv[1])
