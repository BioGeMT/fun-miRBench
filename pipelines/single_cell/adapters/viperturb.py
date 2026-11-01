"""VIPerturb-seq K562 (Zenodo 18460279), raw counts of the three full bins, extracted with seurat_counts.R (the filtered
object is response-filtered and is not used). Each bin carries its own control cells, kept once by cell name. Unit:
10x lane (two lanes in one culture). Usage: python viperturb.py"""

import numpy as np
import pandas as pd
from scipy import io, sparse

from screen_common import DATA, Targets, comparison_rows, plan_cells, write_screen

DATASET = "VIPerturb_K562"
CONTROL = "NO-TARGET"


def main() -> None:
    targets = Targets()
    root = DATA / "viperturb"
    parts, metas, features = [], [], None
    seen = set()
    for bin_name in "ABC":
        folder = root / f"x_{bin_name}"
        genes = (folder / "features.tsv").read_text().splitlines()
        if features is None:
            features = genes
        elif genes != features:
            raise ValueError(f"{folder}: genes differ between bins")
        meta = pd.read_csv(folder / "meta.tsv", sep="\t", dtype=str, keep_default_na=False)
        counts = sparse.csr_matrix(io.mmread(folder / "counts.mtx").T)
        new = np.array([cell not in seen for cell in meta.cell])
        seen |= set(meta.cell)
        parts.append(counts[np.flatnonzero(new)])
        metas.append(meta[new])
    meta = pd.concat(metas, ignore_index=True)
    lane = meta["sample"].str.extract(r"_(L\d+)$", expand=False).fillna(meta["sample"]).to_numpy()
    control = (meta.gene == CONTROL).to_numpy()
    obs = pd.DataFrame({"sample": meta["sample"], "batch": lane, "replicate_unit": lane, "replicate_type": "lane_set",
                        "arm": meta.gene, "is_control": control, "control_type": np.where(control, "non_targeting", ""),
                        "guide": meta.guide, "n_guides": 1, "passed_qc": True})
    kept, rows = plan_cells(obs, targets)
    caveats = ("CRISPRi knockdown in K562 with KRAB-dCas9-MeCP2; one culture, replicate units are the 10x lanes; shallow "
               "per-gene depth; full bins used, not the response-filtered object; cells as the authors filtered them")
    comparisons = comparison_rows(DATASET, kept.arm[~kept.is_control], targets, control_arm=CONTROL, caveats=caveats)
    write_screen(DATASET, sparse.vstack(parts).tocsr()[rows].astype(np.int32), kept, pd.Series(features, dtype=str),
                 comparisons=comparisons, source="authors' raw counts, genome_wide_binA/B/C.RDS, Zenodo 18460279")


if __name__ == "__main__":
    main()
