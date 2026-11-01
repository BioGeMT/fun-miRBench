"""GSE171452 STING-seq v2, K562 with KRAB-dCas9-MeCP2, four 10x lanes, high MOI. Authors' cDNA and GDO counts and their UMI
thresholds per guide. The guides at miRNA loci were found by matching protospacers to GRCh38: SNP-202 is 1.3 kb from MIR142
inside MIR142HG; SNP-199 and SNP-200 flank MIR144/MIR451 (the labels SNP-142, -144 and -451 are list numbers, not miRNAs).
Each lane is a batch and a unit. Usage: python gse171452.py <RAW folder>"""

import gzip
import pathlib
import sys

import numpy as np
import pandas as pd
import scipy.io as scipy_io
from scipy import sparse

from funmirbench.single_cell_store import write_cell_store
from screen_common import ADAPTERS, FINAL, GTF, MIN_CELLS_PER_UNIT

DATASET = "GSE171452_STINGseq"
CRE = {"SNP-202": "MIR142_CRE", "SNP-199": "MIR144_MIR451_CRE_upstream", "SNP-200": "MIR144_MIR451_CRE_downstream"}
TARGETS = {"MIR142_CRE": ("MIR142", "ENSG00000284353;ENSG00000234383", "1.3 kb from MIR142, inside MIR142HG"),
           "MIR144_MIR451_CRE_upstream": ("MIR144;MIR451A", "ENSG00000283819;ENSG00000284565", "about 3.9 kb from MIR451A"),
           "MIR144_MIR451_CRE_downstream": ("MIR144;MIR451A", "ENSG00000283819;ENSG00000284565", "about 1.4 kb from MIR144")}


def mtx(path: pathlib.Path) -> sparse.csr_matrix:
    with gzip.open(path) as handle:
        return sparse.csr_matrix(scipy_io.mmread(handle)).T.tocsr()


def main(root: pathlib.Path) -> None:
    parts, tables, features = [], [], None
    for lane in "ABCD":
        prefix = next(root.glob(f"*STINGseq-v2_cDNA-{lane}_barcodes.mtx.gz")).name.replace("_barcodes.mtx.gz", "")
        accession = prefix.split("_")[0]
        matrix_path = root / f"{prefix}_matrix.mtx.gz"
        if not matrix_path.exists():
            matrix_path = root / f"{accession}_matrix.mtx.gz"  # GEO's corrected matrix of lane B (June 2025)
        counts = mtx(matrix_path)
        barcodes = pd.read_csv(root / f"{prefix}_barcodes.mtx.gz", header=None, dtype=str)[0]
        if counts.shape[0] != len(barcodes):
            print(f"lane {lane} left out: matrix has {counts.shape[0]} cells, barcodes file {len(barcodes)}", flush=True)
            continue
        genes = pd.read_csv(root / f"{prefix}_features.mtx.gz", sep="\t", header=None, dtype=str)
        features = genes
        guide_prefix = next(root.glob(f"*STINGseq-v2_GDO-{lane}_barcodes.mtx.gz")).name.replace("_barcodes.mtx.gz", "")
        guide_all = mtx(root / f"{guide_prefix}_matrix.mtx.gz")
        guide_barcodes = pd.read_csv(root / f"{guide_prefix}_barcodes.mtx.gz", header=None, dtype=str)[0]
        guide_features = pd.read_csv(root / f"{guide_prefix}_features.mtx.gz", sep="\t", header=None, dtype=str)
        thresholds = pd.read_csv(next(root.glob(f"*STINGseq-v2_GDO-{lane}.umi_thresholds.txt.gz")))
        limit = dict(zip(thresholds["Protospacer"], thresholds["UMI threshold"]))
        columns = np.flatnonzero(guide_features[0].isin(limit).to_numpy())
        names = guide_features[0].to_numpy()[columns]
        guides = guide_all[:, columns].tocoo()
        called = guides.data >= np.array([limit[names[j]] for j in guides.col])
        where = pd.Series(np.arange(len(guide_barcodes)), index=guide_barcodes).reindex(barcodes).to_numpy()
        print(lane, counts.shape, guide_all.shape, len(columns), np.isnan(where).sum(), flush=True)
        per_cell = [[] for _ in range(len(guide_barcodes))]
        for cell, column in zip(guides.row[called], guides.col[called]):
            per_cell[cell].append(names[column])
        arms, numbers = [], []
        for position in where:
            mine = per_cell[int(position)] if position == position else []
            labels = {CRE.get(name.rsplit("-", 1)[0], name.rsplit("-", 1)[0]) for name in mine}
            arms.append(";".join(sorted(labels)) or "no_guide")
            numbers.append(len(mine))
        tables.append(pd.DataFrame({"sample": f"lane_{lane}", "batch": f"lane_{lane}", "replicate_unit": f"lane_{lane}",
                                    "replicate_type": "lane_set", "arm": arms, "guide": "", "n_guides": numbers}))
        parts.append(counts)
    obs = pd.concat(tables, ignore_index=True).assign(dataset=DATASET, passed_qc=True)
    obs["is_control"] = obs.arm.str.fullmatch(r"(NTC[^;]*)(;NTC[^;]*)*").to_numpy()
    obs["control_type"] = np.where(obs.is_control, "non_targeting", "")
    print("median guides", obs.n_guides.median(), flush=True)
    rows = [{"comparison_id": f"{DATASET}__{arm}_KO", "dataset": DATASET, "category": "host_gene", "mirna_or_gene": mirna, "orientation": "KO",
             "treated_arm": arm, "control_arm": "NTC", "min_cells": MIN_CELLS_PER_UNIT,
             "caveats": f"CRISPRi of a blood-trait cis-regulatory element {distance} (sequence match of the guides to GRCh38); high MOI: "
                        "control is the cells of the lane without these guides; the mature miRNA was not measured; replicate units are lanes A, C and D: GEO's corrected lane B matrix "
                        "has 19,639 cells against 18,954 barcodes, so its cells cannot be labelled and lane B is left out",
             "target_gene_id": gene_ids, "other_pathway_arms": ""} for arm, (mirna, gene_ids, distance) in TARGETS.items()]
    store_path = FINAL / DATASET / "store" / f"{DATASET}.h5ad"
    store_path.parent.mkdir(parents=True, exist_ok=True)
    write_cell_store(store_path, sparse.vstack(parts).tocsr(), obs, features[0].astype(str), gtf=GTF, symbols=features[1].astype(str),
                     source="authors' cDNA and GDO counts with their UMI thresholds, STING-seq v2, GEO GSE171452")
    pd.DataFrame(rows).to_csv(ADAPTERS.parent / "comparisons" / f"{DATASET}.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main(pathlib.Path(sys.argv[1]))
