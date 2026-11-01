"""GSE210681 TeloHAEC CRISPRi 2,620 targets, 20 lanes from three runs of one culture. Authors' per-lane mtx and guide calls
(noguidemut_calls, Good_singlet = 1). Safe-targeting and negative-control guides are the controls. Each run is one batch
and one replicate unit. Targets labelled <gene>-TSS<n> are pooled into the gene. Usage: python gse210681.py <extracted RAW folder>"""

import pathlib
import re
import sys

import numpy as np
import pandas as pd
from scipy import sparse

from funmirbench.single_cell_store import read_10x_mtx
from screen_common import Targets, comparison_rows, plan_cells, write_screen

DATASET = "GSE210681_TeloHAEC"
CONTROL = "control"
CONTROLS = {"safe-targeting": "safe_targeting", "negative-control": "non_targeting"}


def main(root: pathlib.Path) -> None:
    targets = Targets()
    tables, lanes = [], []
    for calls_path in sorted(root.glob("*noguidemut_calls.tsv.gz")):
        lane = calls_path.name.replace("_noguidemut_calls.tsv.gz", "")
        run = re.search(r"2kG_([0-9A-Za-z]+Dox\d?)_", lane).group(1)
        calls = pd.read_csv(calls_path, sep="\t", index_col=0)
        calls = calls[calls.Good_singlet == 1]
        calls = calls.assign(lane=lane, run=run, gene=calls.Gene.astype(str).str.replace(r"-TSS\d+$", "", regex=True))
        tables.append(calls), lanes.append(lane)
    table = pd.concat(tables, ignore_index=True)
    control = table.gene.isin(CONTROLS).to_numpy()
    obs = pd.DataFrame({"sample": table.lane, "batch": table.run, "replicate_unit": table.run, "replicate_type": "lane_set",
                        "arm": np.where(control, CONTROL, table.gene), "is_control": control,
                        "control_type": table.gene.map(CONTROLS).fillna("").to_numpy(), "guide": table.Guide.astype(str),
                        "n_guides": 1, "passed_qc": True})
    kept, positions = plan_cells(obs, targets)
    chosen = table.iloc[positions]
    parts, order = [], []
    features = None
    for lane in lanes:
        where = np.flatnonzero(chosen.lane == lane)
        counts, features, barcodes = read_10x_mtx(root / f"{lane}_matrix.mtx.gz", root / f"{lane}_features.tsv.gz",
                                                  root / f"{lane}_barcodes.tsv.gz")
        index = pd.Series(np.arange(len(barcodes)), index=barcodes)
        parts.append(counts.T.tocsr()[index.loc[chosen.CBC_10x.to_numpy()[where] + "-1"].to_numpy()])
        order.append(where)
    matrix = sparse.vstack(parts).tocsr()[np.argsort(np.concatenate(order))]
    caveats = ("CRISPRi knockdown; one culture, lanes are technical; replicate units are the three runs of the 20 lanes; "
               "singlets as the authors called them; TSS-specific targets pooled into the gene")
    comparisons = comparison_rows(DATASET, kept.arm[~kept.is_control], targets, control_arm=CONTROL, caveats=caveats)
    write_screen(DATASET, matrix, kept, features[0].astype(str), comparisons=comparisons, symbols=features[1].astype(str),
                 source="authors' raw counts, GSE210681 per-lane mtx and noguidemut_calls")


if __name__ == "__main__":
    main(pathlib.Path(sys.argv[1]))
