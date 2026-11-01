"""Feng et al. 2026 CRISPRi map in HipSci iPSC lines (figshare 27989294), authors' dense RNA UMI CSV (genes x cells) and
cell metadata. Controls are only the cells whose guide call is NonTarget; unassigned cells are not controls. Replicate
unit is the donor line; batch is the pool the lines were cultured in.

Usage: python feng2026.py <gwsf|gwsnf|ts>
"""

import gzip
import re
import sys

import numpy as np
import pandas as pd
from scipy import sparse

from screen_common import DATA, Targets, comparison_rows, plan_cells, write_screen

PREFIX = {"gwsf": "GenomeWideScreen_FitnessGenes", "gwsnf": "GenomeWideScreen_NonFitnessGenes", "ts": "TargetedScreen"}
CONTROL = "NonTarget"


def column_name(cell_id: str) -> str:
    """Cell metadata IDs carry a screen prefix the count files do not: MP-ST-P3-D3_I1_AAA-1 is ST-P3-D3_I1_AAA-1, and in
    the targeted screen PC-P4-D3_I73_AAA-1 is P4_I73_AAA-1."""
    if cell_id.startswith("MP-"):
        return cell_id[3:]
    return re.sub(r"^PC-(P\d+)-D\d+_", r"\1_", cell_id)


def read_cells(path, columns: np.ndarray) -> tuple[sparse.csr_matrix, list[str], list[str]]:
    """Dense genes x cells CSV to cells x genes for the given 0-based cell columns; the first field of a row is the gene."""
    features, parts = [], []
    with gzip.open(path, "rt") as handle:
        next(handle)
        for line in handle:
            name, rest = line.split(",", 1)
            values = np.fromstring(rest, sep=",", dtype=np.int32)
            features.append(name.strip('"'))
            parts.append(sparse.csc_matrix(values[columns].reshape(-1, 1)))
    return sparse.hstack(parts, format="csr", dtype=np.int32), features


def main(key: str) -> None:
    dataset = {"gwsf": "Feng2026_fitness", "gwsnf": "Feng2026_nonfitness", "ts": "Feng2026_targeted"}[key]
    targets = Targets()
    root = DATA / "feng2026"
    meta = pd.read_csv(root / f"{PREFIX[key]}_Cell-Metadata.tsv.gz", sep="\t", dtype=str, keep_default_na=False)
    target = meta.Guide_Call.str.rsplit("_", n=1).str[0]
    control = target == CONTROL
    pools = meta.Batch.str.extract(r"(Pool\d+)", expand=False).fillna(meta.Batch)
    obs = pd.DataFrame({"sample": meta.Batch, "batch": pools, "replicate_unit": meta.Cell_Line, "replicate_type": "biological",
                        "arm": target, "is_control": control.to_numpy(), "control_type": np.where(control, "non_targeting", ""),
                        "guide": meta.Guide_Call, "n_guides": 1, "passed_qc": True})
    obs.index = meta.Cell_ID
    kept, positions = plan_cells(obs, targets)
    kept_ids = meta.Cell_ID.to_numpy()[positions]
    path = root / f"{PREFIX[key]}_RNA-UMI-Counts.csv.gz"
    with gzip.open(path, "rt") as handle:
        header = [value.strip('"') for value in next(handle).rstrip("\n").split(",")][1:]
    column_of = {name: position for position, name in enumerate(header)}
    columns = np.array([column_of[column_name(cell)] for cell in kept_ids])
    counts, features = read_cells(path, columns)
    ids = pd.Series(features).str.split(":")
    caveats = (f"CRISPRi knockdown; controls are the {int(kept.is_control.sum())} NonTarget cells kept, unassigned cells "
               "are not used; replicate units are donor lines within pools; cells as the authors filtered them")
    comparisons = comparison_rows(dataset, kept.arm[~kept.is_control], targets, control_arm=CONTROL, caveats=caveats)
    write_screen(dataset, counts, kept, ids.str[0], comparisons=comparisons, symbols=ids.str[1],
                 source=f"authors' raw counts {path.name}, figshare 27989294")


if __name__ == "__main__":
    main(sys.argv[1])
