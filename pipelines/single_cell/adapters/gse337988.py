"""GSE337988 DLD-1 CRISPRi 5,000-gene sub-library at three MOI levels, authors' processed objects: raw counts in assay001 (cells after
their QC and hash demultiplexing) and their guide calls (perturbationIndicator; extracted by gse337988_extract.R). Each 10x channel is a
batch and a unit. Counts and gene names are checked against the Cell Ranger h5 on one cell. Usage: python gse337988.py <High|Med|Low>"""

import sys

import h5py
import numpy as np
import pandas as pd
from scipy import io, sparse

from screen_common import DATA, Targets, comparison_rows, h5_rows, pathway_arms, plan_cells, read_h5_with_guides, target_arms, write_screen

CONTROL = "control"


def main(moi: str) -> None:
    dataset = f"GSE337988_DLD1_{moi}"
    root = DATA / "GSE337988"
    targets = Targets()
    folder = root / f"x_{moi}"
    cells = pd.read_csv(folder / "cells.tsv", sep="\t", dtype=str, keep_default_na=False)
    guide_names = (folder / "guides.tsv").read_text().splitlines()
    guide_cells = (folder / "guide_cells.tsv").read_text().splitlines()
    if guide_cells != cells.cell.tolist():
        raise ValueError("guide calls and cells differ in order")
    index_map = pd.read_csv(root / "GSE337988_NGS6475_crispr_demux_index_map_prod.csv.gz")
    mapping = dict(zip(index_map["index"], index_map["sample"]))
    target = np.array([CONTROL if mapping.get(name, name.split("_")[0]) == "NTC" else mapping.get(name, name.split("_")[0]) for name in guide_names], dtype=object)
    indicator = sparse.csr_matrix(io.mmread(folder / "guides.mtx").T)
    arms = [";".join(sorted(set(target[indicator.indices[indicator.indptr[cell]:indicator.indptr[cell + 1]]]))) or "no_guide"
            for cell in range(indicator.shape[0])]
    channel = cells.cell.str.split("-").str[0].to_numpy()
    control = np.array([arm == CONTROL for arm in arms])
    obs = pd.DataFrame({"sample": channel, "batch": channel, "replicate_unit": channel, "replicate_type": "lane_set", "arm": arms,
                        "is_control": control, "control_type": np.where(control, "non_targeting", ""), "guide": "",
                        "n_guides": np.diff(indicator.indptr), "passed_qc": True})
    median_guides = obs.n_guides.median()
    kept, rows = plan_cells(obs, targets, all_cells=median_guides > 1.5)
    with h5py.File(root / f"GSE337988_sublib2_processed_objects_{moi}_assays.h5", "r") as handle:
        counts = h5_rows(handle, rows, key="assay001", block=2000)
    gene_ids = pd.Series((folder / "genes.tsv").read_text().splitlines(), dtype=str)
    sample, barcode = cells.cell.iloc[rows[0]].split("-")
    library = cells.NGS_ID.iloc[0]
    library_of = pd.read_csv(root / "GSE337988_sublib2_samplesheet.csv.gz")
    library = library_of[(library_of.Sample == sample) & (library_of.DataType == "rna")].Library.iloc[0]
    whole, table, barcodes, _, _ = read_h5_with_guides(next((root / "x").glob(f"*_{library}.filtered_feature_bc_matrix.h5")))
    position = np.flatnonzero(pd.Series(barcodes).str.replace("-1", "", regex=False) == barcode)[0]
    order = pd.Series(np.arange(len(gene_ids)), index=gene_ids).reindex(table.id).to_numpy()
    if not np.array_equal(whole[position].toarray().ravel(), counts[0].toarray().ravel()[order]):
        raise ValueError("assay001 and the Cell Ranger h5 differ on the check cell")
    print("gene names checked on one cell", flush=True)
    names = table.set_index("id").name.reindex(gene_ids).astype(str).reset_index(drop=True)
    found = target_arms(kept.arm, targets)
    caveats = (f"CRISPRi knockdown, 5,000-gene sub-library, MOI {moi}; median {median_guides:g} guides per cell; replicate units are 10x channels; "
               "cells and guide calls as the authors deposited them")
    comparisons = comparison_rows(dataset, found, targets, control_arm=CONTROL, caveats=caveats,
                                  other_pathway_arms=pathway_arms(kept.arm, targets) if median_guides > 1.5 else "")
    write_screen(dataset, counts, kept, gene_ids, comparisons=comparisons, symbols=names,
                 source=f"authors' processed object, GSE337988 sublib2 {moi}")


if __name__ == "__main__":
    main(sys.argv[1])
