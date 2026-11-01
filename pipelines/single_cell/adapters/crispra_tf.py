"""CRISPRa transcription-factor screens, RPE1 (Zenodo 15213619) and Hs27 (Zenodo 15200179), authors' final-population h5ad of
raw counts. Cells may carry several guides; thresholded_features lists them, and the arm holds their targets. Each GEM
group is a batch and a replicate unit. Gain of function: orientation OE. Usage: python crispra_tf.py <rpe1|hs27>"""

import sys

import h5py
import numpy as np
import pandas as pd

from screen_common import DATA, Targets, comparison_rows, h5_frame, h5_rows, pathway_arms, plan_cells, target_arms, write_screen

SCREENS = {"rpe1": ("CRISPRaTF_RPE1", "crispra_rpe1/RPE1_CRISPRa_final_population.h5ad", "Zenodo 15213619"),
           "hs27": ("CRISPRaTF_Hs27", "crispra_hs27/fibroblast_CRISPRa_final_pop.h5ad", "Zenodo 15200179")}
CONTROL = "non_targeting"


def target_of(guide: str) -> str:
    return CONTROL if guide.startswith(("non_targeting", "non-targeting", "off_target", "off-target")) else guide.rsplit("_", 1)[0]


def main(key: str) -> None:
    dataset, relative, deposit_name = SCREENS[key]
    targets = Targets()
    with h5py.File(DATA / relative, "r") as handle:
        deposit = h5_frame(handle["obs"])
        var = h5_frame(handle["var"])
        guides = deposit.thresholded_features.astype(str).str.split("|")
        arms = [";".join(sorted({target_of(guide) for guide in value})) for value in guides]
        control = np.array([arm == CONTROL for arm in arms])
        batch = "gem_" + deposit.gem_group.astype(str).to_numpy()
        obs = pd.DataFrame({"sample": batch, "batch": batch, "replicate_unit": batch, "replicate_type": "lane_set", "arm": arms,
                            "is_control": control, "control_type": np.where(control, "non_targeting_or_off_target", ""),
                            "guide": deposit.thresholded_features.astype(str).to_numpy(),
                            "n_guides": [len(value) for value in guides], "passed_qc": True})
        median_guides = obs.n_guides.median()
        kept, rows = plan_cells(obs, targets, all_cells=median_guides > 1.5)
        counts = h5_rows(handle, rows)
    found = target_arms(kept.arm, targets)
    caveats = ("CRISPRa gain of function; replicate units are single GEM groups; RAX in this library is the RAX gene, not "
               "PRKRA; cells as the authors filtered them")
    comparisons = comparison_rows(dataset, found, targets, control_arm=CONTROL, caveats=caveats, raised=True,
                                  other_pathway_arms=pathway_arms(kept.arm, targets))
    features = var.gene_id.astype(str) if "gene_id" in var else pd.Series(var.index, dtype=str)
    write_screen(dataset, counts, kept, features.reset_index(drop=True), comparisons=comparisons,
                 symbols=pd.Series(var.gene_name.astype(str).to_numpy()), source=f"authors' raw counts, {relative.split('/')[-1]}, {deposit_name}")


if __name__ == "__main__":
    main(sys.argv[1])
