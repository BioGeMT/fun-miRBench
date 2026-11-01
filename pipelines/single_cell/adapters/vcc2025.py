"""Virtual Cell Challenge 2025 H1 CRISPRi (Arc), authors' raw counts in three files. The same 38,176 control cells are in all
three files, so controls are read from Training only. Unit = lane set of eight Flex batches. Usage: python vcc2025.py"""

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

from screen_common import DATA, Targets, comparison_rows, h5_frame, h5_rows, lane_sets, natural_order, plan_cells, write_screen

DATASET = "VCC2025_H1"
CONTROL = "non-targeting"
FILES = ("adata_Training.h5ad", "adata_Validation.h5ad", "adata_Test.h5ad")


def main() -> None:
    targets = Targets()
    tables, parts, var_index = [], [], None
    seen_controls = False
    all_batches = set()
    for name in FILES:
        with h5py.File(DATA / "vcc2025" / name, "r") as handle:
            all_batches |= set(decode_batches(handle))
    numbers = dict(zip(sorted(all_batches), natural_order(sorted(all_batches))))
    for name in FILES:
        with h5py.File(DATA / "vcc2025" / name, "r") as handle:
            deposit = h5_frame(handle["obs"])
            var = h5_frame(handle["var"])
            if var_index is None:
                var_index = var.index
            elif not var.index.equals(var_index):
                raise ValueError(f"{name}: genes differ between files")
            keep = deposit.target_gene.astype(str).map(targets.symbol_of).notna().to_numpy()
            if not seen_controls:
                keep |= deposit.target_gene.to_numpy() == CONTROL
                seen_controls = True
            rows = np.flatnonzero(keep)
            parts.append(h5_rows(handle, rows))
            tables.append(deposit.iloc[rows])
    table = pd.concat(tables)
    if table.index.duplicated().any():
        raise ValueError("a cell appears twice")
    units, rule = lane_sets(table.batch.map(numbers).to_numpy())
    control = table.target_gene.to_numpy() == CONTROL
    obs = pd.DataFrame({"sample": table.batch.astype(str).to_numpy(), "batch": units, "replicate_unit": units,
                        "replicate_type": "lane_set", "arm": table.target_gene.astype(str).to_numpy(), "is_control": control,
                        "control_type": np.where(control, "non_targeting", ""), "guide": table.guide_id.astype(str).to_numpy(),
                        "n_guides": 1, "passed_qc": True})
    kept, positions = plan_cells(obs, targets)
    caveats = (f"CRISPRi knockdown; {rule}; controls counted once though in all three files; perturbations were chosen by "
               "the authors for effect strength; cells as the authors filtered them")
    comparisons = comparison_rows(DATASET, kept.arm[~kept.is_control], targets, control_arm=CONTROL, caveats=caveats)
    write_screen(DATASET, sparse.vstack(parts).tocsr()[positions], kept, pd.Series(var_index, dtype=str), comparisons=comparisons,
                 source="authors' raw counts, adata_Training/Validation/Test.h5ad, Arc virtual-cell-challenge/2025")


def decode_batches(handle) -> np.ndarray:
    return pd.unique(h5_frame(handle["obs"]).batch.astype(str))


if __name__ == "__main__":
    main()
