"""Screens deposited as 10x h5 files with guide counts and no per-cell guide calls: GSE327644 (macrophages, Ctrl and LPS),
GSE197452 (A375), GSE337988 (DLD-1, low, medium and high MOI). Guides are called with screen_common.call_guides. Each
10x lane is a batch and a replicate unit; with one lane, three technical pseudo-replicates.

Usage: python guide_h5_screens.py <dataset> <folder with the extracted RAW tar>
"""

import gzip
import pathlib
import re
import sys

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

from funmirbench.single_cell_store import pseudo_replicate_units
from screen_common import (DATA, Targets, call_guides, comparison_rows, decode, pathway_arms, plan_cells,
                           read_h5_with_guides, target_arms, write_screen)

CONTROL = "control"


def macrophage_lanes(root: pathlib.Path, condition: str) -> list[tuple[str, pathlib.Path, pathlib.Path | None]]:
    return [(re.search(r"Mac_(L\d+)_", path.name).group(1), path, None)
            for path in sorted(root.glob(f"*_Mac_{condition}_sample_filtered_feature_bc_matrix.h5"))]


def dld1_lanes(root: pathlib.Path, moi: str) -> list[tuple[str, pathlib.Path, pathlib.Path | None]]:
    lanes = []
    with gzip.open(DATA / "GSE337988" / "GSE337988_series_matrix.txt.gz", "rt", errors="replace") as handle:
        for line in handle:
            if line.startswith("!Sample_title"):
                titles = [value.strip('"') for value in line.rstrip("\n").split("\t")[1:]]
            if line.startswith("!Sample_geo_accession"):
                accessions = [value.strip('"') for value in line.rstrip("\n").split("\t")[1:]]
    for accession, title in zip(accessions, titles):
        if f"MOI {moi}" in title and title.endswith("GEX") and "sublib" in title:
            gex = sorted(root.glob(f"{accession}_*.filtered_feature_bc_matrix.h5"))
            crispr = sorted(root.glob(f"{accession}_*.crispr.umi_correct.h5"))
            if len(gex) == 1 and len(crispr) == 1:
                lanes.append((accession, gex[0], crispr[0]))
    return lanes


def read_separate_guides(path: pathlib.Path, barcodes: np.ndarray) -> tuple[sparse.csr_matrix, np.ndarray]:
    """Guide counts in their own h5, aligned to the given barcodes (zero for barcodes it lacks)."""
    with h5py.File(path, "r") as handle:
        group = handle["matrix"]
        whole = sparse.csc_matrix((group["data"][()], group["indices"][()], group["indptr"][()]),
                                  shape=tuple(group["shape"][()])).T.tocsr()
        own = pd.Series(np.arange(whole.shape[0]), index=[value.split("-")[0] for value in decode(group["barcodes"][()])])
        names = decode(group["features/name"][()])
    bare = [value.split("-")[0] for value in barcodes]
    where = own.reindex(bare).to_numpy()
    found = ~np.isnan(where)
    picked = whole[where[found].astype(int)]
    index = sparse.csr_matrix((np.ones(found.sum(), dtype=np.int32), (np.flatnonzero(found), np.arange(found.sum()))),
                              shape=(len(bare), found.sum()))
    return (index @ picked).tocsr().astype(np.int32), names


def build(dataset: str, lanes, guide_target, *, caveats: str, source: str, raised: bool = False) -> None:
    targets = Targets()
    parts, tables, features = [], [], None
    for lane, gex_path, crispr_path in lanes:
        counts, genes, barcodes, guides, names = read_h5_with_guides(gex_path)
        if crispr_path is not None:
            guides, names = read_separate_guides(crispr_path, barcodes)
        if features is None:
            features = genes
        elif not genes.id.equals(features.id):
            raise ValueError(f"{gex_path}: genes differ between lanes")
        arms, numbers, called = call_guides(guides, np.array([guide_target(name) for name in names], dtype=object))
        parts.append(counts)
        tables.append(pd.DataFrame({"sample": lane, "batch": lane, "replicate_unit": lane, "replicate_type": "lane_set",
                                    "arm": arms, "guide": called, "n_guides": numbers}))
        print(f"{lane}: {counts.shape[0]} cells, median guides {np.median(numbers):g}", flush=True)
    obs = pd.concat(tables, ignore_index=True)
    obs["arm"] = obs.arm.where(obs.arm != "", "no_guide")
    obs["is_control"] = (obs.arm == CONTROL).to_numpy()
    obs["control_type"] = np.where(obs.is_control, "non_targeting", "")
    obs["passed_qc"] = True
    median_guides = obs.n_guides.median()
    kept, rows = plan_cells(obs, targets, all_cells=median_guides > 1.5)
    if len(lanes) == 1:
        kept["replicate_unit"] = pseudo_replicate_units(kept["sample"], kept.arm, kept.passed_qc)
        kept["replicate_type"] = "technical"
    found = target_arms(kept.arm, targets)
    comparisons = comparison_rows(dataset, found, targets, control_arm=CONTROL, caveats=f"{caveats}; guides called from guide counts with "
                                  "a binomial test against the lane's background, BH 0.05, at least 3 UMIs; the authors' calls are not deposited",
                                  raised=raised, other_pathway_arms=pathway_arms(kept.arm, targets) if median_guides > 1.5 else "")
    write_screen(dataset, sparse.vstack(parts).tocsr()[rows], kept, features.id, comparisons=comparisons, symbols=features.name, source=source)


def main(dataset: str, root: pathlib.Path) -> None:
    if dataset.startswith("GSE327644_macrophage"):
        template = pd.read_csv(DATA / "GSE327644" / "GSE327644_PRJNA1450051_marson_mac_perturbseq_guide_template.txt.gz", sep="\t")
        mapping = dict(zip(template.guide_id, np.where(template.guide_role == "Control", CONTROL, template.guide_target_gene_name)))
        condition = dataset.rsplit("_", 1)[1]
        # h5 guide names outside the template (for example JUN-1) are <gene>-<number>
        build(dataset, macrophage_lanes(root, condition), lambda name: mapping.get(name, name.rsplit("-", 1)[0]),
              caveats=f"Cas9 knockout delivered with virus-like particles, macrophages {condition}; one donor on the sample pages read; "
                      "replicate units are single lanes",
              source="authors' 10x counts, GSE327644 per-lane filtered h5")
    elif dataset == "GSE197452_A375":
        control = re.compile(r"^(NO_SITE|ONE_NON-GENE_SITE)_")
        build(dataset, [("pert_illumina", next(root.glob("*filtered_feature_bc_matrix.pert.ill.h5")), None)],
              lambda name: CONTROL if control.match(name) else name.rsplit("_", 1)[0],
              caveats="Cas9 knockout, transcription factor and chromatin regulator guide library (Simmons 2022 methods), sequencing method comparison; one library, so "
                      "replicate units are three technical pseudo-replicates: labels mean changed in this culture; Illumina reads used",
              source="authors' 10x counts, filtered_feature_bc_matrix.pert.ill.h5, GEO GSE197452")
    elif dataset.startswith("GSE337988_DLD1"):
        moi = dataset.rsplit("_", 1)[1].capitalize()
        index_map = pd.read_csv(DATA / "GSE337988" / "GSE337988_NGS6475_crispr_demux_index_map_prod.csv.gz")
        mapping = dict(zip(index_map["index"], index_map["sample"]))

        def target(name: str) -> str:
            gene = mapping.get(name, name.split("_")[0])
            return CONTROL if gene == "NTC" else gene

        build(dataset, dld1_lanes(root, moi), target,
              caveats=f"CRISPRi knockdown, 5,000-gene sub-library, MOI {moi}; replicate units are single 10x channels",
              source=f"authors' 10x counts and crispr h5, GEO GSE337988, MOI {moi} channels")
    else:
        raise ValueError(dataset)


if __name__ == "__main__":
    main(sys.argv[1], pathlib.Path(sys.argv[2]))
