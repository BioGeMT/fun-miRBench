"""DepMap Perturb-seq pilot in 16 cancer lines (figshare 33273600): CRISPR knockout, 100 genes, 2 guides each. Authors' 10x
files and their cell table (processed/all_cells_table.csv: assigned guide, pass_qc, control_guide). Controls are the
authors' control guides, which cut olfactory-receptor genes. One sample per line, so three technical pseudo-replicates.

Usage: python cancer16.py <cell line>
"""

import sys

import pandas as pd

from funmirbench.single_cell_store import pseudo_replicate_units, read_10x_h5
from screen_common import DATA, Targets, comparison_rows, plan_cells, write_screen

ROOT = DATA / "cancer16" / "x"
CONTROL = "control_guide"


def main(line: str) -> None:
    dataset = f"Cancer16_{line}"
    targets = Targets()
    manifest = pd.read_csv(ROOT / "single_cell_data" / "raw_data_file_manifest.csv").set_index("CellLine")
    counts, features, barcodes = read_10x_h5(ROOT / "single_cell_data" / manifest.loc[line, "GeneExpression"])
    cells = pd.read_csv(ROOT / "processed" / "all_cells_table.csv", low_memory=False)
    cells = cells[cells.cell_line == line].set_index("barcode")
    barcodes = pd.Series(barcodes)
    if not barcodes.isin(cells.index).any():  # the authors' cell table drops the -1 suffix for most lines
        barcodes = barcodes.str.replace(r"-\d+$", "", regex=True)
    table = cells.reindex(barcodes.to_numpy())
    known = table.assigned_ko.notna().to_numpy() & (table.assigned_ko != "unperturbed").to_numpy()
    rows = known.nonzero()[0]
    table = table.iloc[rows]
    control = table.control_guide.astype(bool).to_numpy()
    arm = table.assigned_ko.astype(str).where(~control, CONTROL).to_numpy()
    obs = pd.DataFrame({"sample": line, "batch": line, "replicate_unit": line, "replicate_type": "technical", "arm": arm,
                        "is_control": control, "control_type": ["cutting_control" if flag else "" for flag in control],
                        "guide": table.assigned_grna.astype(str).to_numpy(), "n_guides": 1,
                        "passed_qc": table.pass_qc.astype(bool).to_numpy()})
    kept, positions = plan_cells(obs, targets)
    kept["replicate_unit"] = pseudo_replicate_units(kept["sample"], kept.arm, kept.passed_qc)
    matrix = counts.T.tocsr()[rows[positions]]
    caveats = ("Cas9 knockout; controls cut olfactory-receptor genes; one library per line, so replicate units are three "
               "technical pseudo-replicates: labels mean changed in this culture; cell QC and guide calls from the "
               "authors' cell table")
    comparisons = comparison_rows(dataset, kept.arm[~kept.is_control], targets, control_arm=CONTROL, caveats=caveats)
    write_screen(dataset, matrix, kept, pd.Series(features[0], dtype=str), comparisons=comparisons,
                 symbols=pd.Series(features[1], dtype=str),
                 source=f"authors' 10x counts {manifest.loc[line, 'GeneExpression']}, figshare 33273600")


if __name__ == "__main__":
    main(sys.argv[1])
