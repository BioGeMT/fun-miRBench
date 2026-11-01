"""Screens read from scPerturb copies (pertpy mirror), which hold the authors' counts with per-cell guide labels: Tian 2021
CRISPRi neurons (GSE152988, whose GEO files have no guide features), Frangieh 2021 melanoma (SCP1064). Only cells with
one guide or none are used, as the deposits label them. Usage: python scperturb_screens.py <GSE152988_neurons|SCP1064_melanoma>
"""

import sys

import anndata as ad
import numpy as np
import pandas as pd

from screen_common import DATA, Targets, comparison_rows, plan_cells, write_screen

CONTROL = "control"
SCREENS = {
    "GSE152988_neurons": dict(path="GSE152988/tian_2021_crispri.h5ad", batch="batch", symbol="_index",
                              caveats="CRISPRi knockdown in iPSC-derived neurons; replicate units are the four libraries; "
                                      "guide labels come from scPerturb, as the GEO files hold no guide features; few controls",
                              source="scPerturb TianKampmann2021_CRISPRi via the pertpy mirror; counts as in GSE152988"),
    "SCP1064_melanoma": dict(path="SCP1064/frangieh_2021_rna.h5ad", batch="perturbation_2", symbol="gene_symbol",
                             caveats="Cas9 knockout in patient-derived melanoma cells; replicate units are the three culture "
                                     "conditions (control, IFN-gamma, co-culture), each with its own controls; a Cas9 cut in a "
                                     "lncRNA gene may not remove the transcript",
                             source="scPerturb FrangiehIzar2021_RNA via the pertpy mirror; counts as in SCP1064"),
}


def main(dataset: str) -> None:
    screen = SCREENS[dataset]
    targets = Targets()
    deposit = ad.read_h5ad(DATA / screen["path"])
    meta = deposit.obs
    labels = meta.perturbation.astype(str)
    print(labels[meta.nperts == 0].value_counts().head(3).to_dict(), flush=True)
    control = (meta.nperts == 0).to_numpy()
    batch = meta[screen["batch"]].astype(str).to_numpy()
    obs = pd.DataFrame({"sample": batch, "batch": batch, "replicate_unit": batch, "replicate_type": "lane_set",
                        "arm": np.where(control, CONTROL, labels), "is_control": control,
                        "control_type": np.where(control, "non_targeting", ""), "guide": meta.guide_id.astype(str).to_numpy(),
                        "n_guides": meta.nperts.astype(int).to_numpy(), "passed_qc": True})
    kept, rows = plan_cells(obs, targets)
    counts = deposit.X.tocsr()[rows]
    if counts.data.size and not np.array_equal(counts.data, np.round(counts.data)):
        raise ValueError("X is not integer")
    symbols = deposit.var[screen["symbol"]].astype(str) if screen["symbol"] in deposit.var else deposit.var_names.astype(str)
    comparisons = comparison_rows(dataset, kept.arm[~kept.is_control], targets, control_arm=CONTROL, caveats=screen["caveats"])
    ensembl = pd.Series(deposit.var.ensembl_id.astype(str).to_numpy())
    # the Frangieh file has Ensembl IDs for under 90% of its features, so it maps by symbol
    features = ensembl if ensembl.str.startswith("ENSG").mean() >= 0.9 else pd.Series(np.asarray(symbols), dtype=str)
    write_screen(dataset, counts, kept, features, comparisons=comparisons, symbols=pd.Series(np.asarray(symbols)), source=screen["source"])


if __name__ == "__main__":
    main(sys.argv[1])
