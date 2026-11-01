"""IGVF CC-Perturb-seq deposits (Engreitz lab): WTC-11 endothelial differentiation D0 and D1, TeloHAEC genome-wide by condition.
Each h5ad holds layers/raw_counts and the pipeline's binary guide assignment in obsm; high guide load. The deposit's replicate is
a batch and a unit. SAFE and non-targeting guides are controls. TeloHAEC has no raw layer: its X holds fractional estimates, kept in the store
and rounded after summing. Usage: python igvf_screens.py <dataset> [store folder under FINAL]"""

import sys

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

from screen_common import DATA, Targets, comparison_rows, decode, h5_frame, h5_rows, pathway_arms, target_arms, write_screen

CONTROL = "control"
SCREENS = {
    "IGVF_WTC11_endo_D0": (["igvf_endo/IGVFFI7756RKHT.h5ad"], "biological_sample", "IGVFDS5257JGQO", "WTC-11 iPSC, day 0 of endothelial differentiation"),
    "IGVF_WTC11_endo_D1": (["igvf_endo/IGVFFI4783YGLZ.h5ad"], "biological_sample", "IGVFDS3530HHEU", "WTC-11 iPSC, day 1 of endothelial differentiation"),
    "IGVF_TeloHAEC_resting": (["igvf_telohaec/IGVFFI2522NCFR.h5ad", "igvf_telohaec/IGVFFI0721PXFB.h5ad"], None, "IGVFDS2068OODO, IGVFDS6017INUS", "TeloHAEC resting"),
    "IGVF_TeloHAEC_il1b": (["igvf_telohaec/IGVFFI1663PLGP.h5ad", "igvf_telohaec/IGVFFI5515HCMW.h5ad"], None, "IGVFDS5371SLCS, IGVFDS2589XTPD", "TeloHAEC with IL-1 beta"),
    "IGVF_TeloHAEC_oxldl": (["igvf_telohaec/IGVFFI2608RWDT.h5ad", "igvf_telohaec/IGVFFI2649OSOD.h5ad"], None, "IGVFDS9554DAPM, IGVFDS4026PWBE", "TeloHAEC with oxidised LDL"),
}


def target_of(guide: str) -> str:
    """Endothelial library: GENE_strand_position; TeloHAEC library: SEQUENCE_GENE_source_guide, with non-targeting and safe-targeting names."""
    parts = guide.split("_")
    gene = parts[1] if len(parts) > 2 and set(parts[0]) <= set("ACGT") else parts[0]
    return CONTROL if gene.lower().startswith(("safe", "non-targeting")) else gene


def read_assignment(handle: h5py.File) -> sparse.csr_matrix:
    node = handle["obsm/guide_assignment"]
    if isinstance(node, h5py.Dataset):
        return sparse.csr_matrix(node[()] > 0)
    return sparse.csr_matrix((node["data"][()] > 0, node["indices"][()], node["indptr"][()]), shape=tuple(node.attrs["shape"]))


def main(dataset: str, folder: str | None = None) -> None:
    files, unit_column, accession, context = SCREENS[dataset]
    targets = Targets()
    parts, tables, features = [], [], None
    fractional = False
    for number, relative in enumerate(files, start=1):
        with h5py.File(DATA / relative, "r") as handle:
            meta = h5_frame(handle["obs"])
            var = h5_frame(handle["var"])
            names = decode(handle["uns/guide_names"][()])
            assignment = read_assignment(handle)
            if "raw_counts" in handle["layers"]:
                counts = h5_rows(handle, np.arange(meta.shape[0]), key="layers/raw_counts")
            else:  # TeloHAEC: X holds fractional estimates (a quarter of the values, reads shared between genes; cell totals are whole numbers)
                counts = h5_rows(handle, np.arange(meta.shape[0]), integer=False).astype(np.float32)
                fractional = True
        guide_targets = np.array([target_of(name) for name in names], dtype=object)
        arms = [";".join(sorted(set(guide_targets[assignment.indices[assignment.indptr[cell]:assignment.indptr[cell + 1]]]))) or "no_guide"
                for cell in range(assignment.shape[0])]
        unit = meta[unit_column].astype(str).to_numpy() if unit_column else np.full(meta.shape[0], f"rep{number}", dtype=object)
        if features is None:
            features = var
        elif not var.index.equals(features.index):
            raise ValueError(f"{relative}: genes differ")
        parts.append(counts)
        tables.append(pd.DataFrame({"sample": meta.sample_id.astype(str).to_numpy() if "sample_id" in meta else unit, "batch": unit,
                                    "replicate_unit": unit, "replicate_type": "biological", "arm": arms, "guide": "",
                                    "n_guides": np.diff(assignment.indptr), "passed_qc": True}))
    obs = pd.concat(tables, ignore_index=True)
    obs["arm"] = [";".join(targets.symbol_of(part) or part for part in value.split(";")) for value in obs.arm]
    obs["is_control"] = (obs.arm == CONTROL).to_numpy()
    obs["control_type"] = np.where(obs.is_control, "non_targeting_or_safe_targeting", "")
    found = target_arms(obs.arm, targets)
    caveats = (f"CRISPRi knockdown, {context}; high guide load, median {obs.n_guides.median():g} guides per cell; guide assignment as the "
               "IGVF pipeline deposited; replicate units are the deposit's replicates")
    if fractional:
        caveats += "; the deposit holds fractional expression estimates per cell (reads shared between genes), rounded after summing"
    comparisons = comparison_rows(dataset, found, targets, control_arm=CONTROL, caveats=caveats, other_pathway_arms=pathway_arms(obs.arm, targets))
    symbol_column = next((column for column in ("gene_symbol", "gene_name") if column in var), None)
    symbols = var[symbol_column].astype(str) if symbol_column else pd.Series(var.index, dtype=str)
    id_column = next((column for column in ("ensembl_id", "gene_id") if column in var), None)
    feature_ids = var[id_column].astype(str) if id_column else pd.Series(var.index, dtype=str)
    write_screen(dataset, sparse.vstack(parts).tocsr(), obs, feature_ids.reset_index(drop=True), comparisons=comparisons,
                 symbols=symbols.reset_index(drop=True), folder=folder, fractional=fractional,
                 source=f"authors' {'X, fractional estimates, no raw layer' if fractional else 'raw_counts layer'}, IGVF {accession}")


if __name__ == "__main__":
    main(*sys.argv[1:])
