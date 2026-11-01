"""GSE249595 Jurkat genome-wide CRISPRi, 16 channels, high MOI, cell hashing for untreated and activated cells . Authors' per-channel mtx
of raw UMIs over all barcodes.
Channels form lane sets of eight. Usage: python gse249595.py <untreated|activated> <RAW folder>"""

import gzip
import pathlib
import re
import sys

import numpy as np
import pandas as pd
import scipy.io as scipy_io
from scipy import sparse

from screen_common import Targets, call_guides, comparison_rows, lane_sets, natural_order, pathway_arms, target_arms, write_screen

CONTROL = "CTRL"
EXPECTED_CELLS_PER_CHANNEL = 586_000 // 16  # the paper's cell count over its 16 channels; the deposit does not give the expected number


def mtx(path: pathlib.Path) -> sparse.csr_matrix:
    with gzip.open(path) as handle:
        return sparse.csr_matrix(scipy_io.mmread(handle)).T.tocsr()


def main(condition: str, root: pathlib.Path) -> None:
    dataset = f"GSE249595_Jurkat_{condition}"
    targets = Targets()
    prefixes = sorted({path.name.split("_transcriptome")[0] for path in root.glob("*_transcriptome_matrix.mtx.gz")},
                      key=lambda name: int(re.search(r"channel(\d+)", name).group(1)))
    parts, tables, features = [], [], None
    for prefix in prefixes:
        channel = re.search(r"channel\d+", prefix).group(0)
        labels = mtx(root / f"{prefix}_labels_matrix.mtx.gz")
        label_names = pd.read_csv(root / f"{prefix}_labels_features.tsv.gz", sep="\t", header=None, dtype=str)[0].to_numpy()
        whole = mtx(root / f"{prefix}_transcriptome_matrix.mtx.gz")
        totals = np.asarray(whole.sum(axis=1)).ravel()
        # the authors' cell filter (GEO): UMIs at least 10% of the 99th percentile of the top expected cells
        top = np.sort(totals)[-EXPECTED_CELLS_PER_CHANNEL:]
        rows = np.flatnonzero(totals >= 0.1 * np.percentile(top, 99))
        hashes = labels[rows].toarray()
        # the matrices hold raw hash and guide UMIs; a hash is called with the top hash at least half of the cell's hash UMIs
        one_hash = hashes.max(axis=1) >= 0.5 * np.maximum(hashes.sum(axis=1), 1)
        is_wanted = one_hash & (hashes.sum(axis=1) > 0) & np.array([label_names[value].startswith(condition) for value in hashes.argmax(axis=1)])
        rows = rows[is_wanted]
        counts = whole[rows]
        guides = mtx(root / f"{prefix}_guides_matrix.mtx.gz")[rows].tocsr()
        guide_names = pd.read_csv(root / f"{prefix}_guides_features.tsv.gz", sep="\t", header=None, dtype=str)[0].to_numpy()
        guide_targets = np.array([name.rsplit("_", 1)[0] for name in guide_names], dtype=object)
        arms, numbers, _ = call_guides(guides, guide_targets)
        arms = [arm or "no_guide" for arm in arms]
        features = pd.read_csv(root / f"{prefix}_transcriptome_features.tsv.gz", sep="\t", header=None, dtype=str)
        parts.append(counts.astype(np.int32))
        tables.append(pd.DataFrame({"channel": channel, "arm": arms, "n_guides": numbers}))
        print(channel, counts.shape, np.median(numbers) if len(numbers) else 0, flush=True)
    table = pd.concat(tables, ignore_index=True)
    units, rule = lane_sets(natural_order(table.channel.to_numpy()))
    table["arm"] = [";".join(targets.symbol_of(part) or part for part in value.split(";")) for value in table.arm]
    control = (table.arm == CONTROL).to_numpy()
    obs = pd.DataFrame({"sample": table.channel, "batch": units, "replicate_unit": units, "replicate_type": "lane_set", "arm": table.arm,
                        "is_control": control, "control_type": np.where(control, "non_targeting", ""), "guide": "", "n_guides": table.n_guides,
                        "passed_qc": True})
    found = target_arms(obs.arm, targets)
    caveats = (f"CRISPRi knockdown, Jurkat {condition}; high MOI, median {obs.n_guides.median():g} guides per cell; {rule}; "
               "the deposit holds raw hash and guide UMIs over all barcodes, so the authors' UMI filter was reapplied with about 36,600 expected "
               "cells per channel, hashes called with the top hash at least half of the cell's hash UMIs, and guides with the authors' binomial test "
               "at BH 0.05 per channel with at least 3 UMIs; their mitochondrial outlier filter gives no threshold and is not reapplied")
    comparisons = comparison_rows(dataset, found, targets, control_arm=CONTROL, caveats=caveats,
                                  other_pathway_arms=pathway_arms(obs.arm, targets))
    write_screen(dataset, sparse.vstack(parts).tocsr(), obs, features[0].astype(str), comparisons=comparisons, symbols=features[1].astype(str),
                 source="authors' STARsolo counts and calls, GSE249595 per-channel mtx")


if __name__ == "__main__":
    main(sys.argv[1], pathlib.Path(sys.argv[2]))
