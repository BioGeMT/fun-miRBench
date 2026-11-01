"""Write the co-occurrence outcomes and studies tables from each study's study.json and groups.tsv, plus the rows that
are not built. Argument: the run's units.tsv."""
import json
import pathlib
import subprocess
import sys

import pandas as pd

FINAL = pathlib.Path("/SCRATCH/dtzim01/fun-experiment-data/single-cell-final/20261031/cooccurrence_library_total")
OUTCOMES = pathlib.Path("/SCRATCH/dtzim01/fun-miRBind-work/single_cell_mirna/final_build/outcomes/cooccurrence.tsv")
DATA = pathlib.Path("/SCRATCH/dtzim01/single_cell_mirna_data")
DEPOSIT_DIR = {"GSE114071_K562": "GSE114071", "GSE226714_lines": "GSE226714", "GSE151334_SmartSeqTotal": "GSE151334",
               "CRA001133_HCC": "CRA001133", "GSE249713_embryo": "GSE249713", "GSE158880_scDASH": "GSE158880",
               "GSE296405_STORMseq": "GSE296405", "GSE176588_HEK293T": "GSE176588", "GSE119984_neuroblastoma": "GSE119984",
               "GSE174293_brain": "GSE174293", "GSE261371_lines": "GSE261371", "GSE261371_PBMC": "GSE261371",
               "TotalX_HEK293T": "totalx", "TotalX_Huh7": "totalx", "TotalX_PBMC": "totalx", "TotalX_brain": "totalx",
               "TabulaSapiens_totalRNA": "tabula_sapiens_total"}
NOT_BUILT = [
    ("SRP500050_DMFDMseq", "blocked", "tier 2: raw reads only, not aligned (70 runs, 25.1 GB, ENA). 35 miRNA and 35 mRNA runs; the mRNA runs titled MDA-MB-21 pair with the miRNA runs titled MDA-MB-231 if that is a typo, so 35 cells, or 25 without them", "same cell per the paper's title and run names; paper not open access", "unknown until aligned"),
    ("GSE81287_cooccurrence", "population_level_only", "different cells for the two assays; not a same-cell study", "no: Small-seq and Smart-seq2 from different cells", "mature arms and precursors, separate tables"),
    ("E-MTAB-10098", "population_level_only", "same embryos, different cells: Small-seq was done on cells collected from a subset of embryos, separately from scM&T-seq (Genome Res 2022 methods). Labels such as C2_s1 and C2_1 look alike but are not the same cell", "no: same embryos only", "mature arms (mature.counts.tsv.gz, not downloaded)"),
    ("GSE102734", "rejected", "microfluidic RT-qPCR of selected assays; miRNA and mRNA panels in separate experiments, not sequencing", "no", "not applicable"),
    ("GSE134004", "rejected", "PSCSR-seq version 1: all 8 samples miRNA-Seq, no mRNA library", "no mRNA", "not applicable"),
    ("GSE162514", "rejected", "all 245 samples miRNA-Seq, no mRNA library", "no mRNA", "not applicable"),
    ("SRP279094", "rejected", "the SRA study of GSE162514; miRNA-seq only", "no mRNA", "not applicable"),
    ("GSE264138", "rejected", "all 33 samples miRNA-Seq, no mRNA library", "no mRNA", "not applicable"),
    ("GSE95218", "rejected", "miRNA-Seq and ncRNA-Seq only (CAS-seq), no mRNA library", "no mRNA", "not applicable"),
    ("SRP453419", "rejected", "ENA: miRNA-seq runs only, human lines and mouse", "no mRNA", "not applicable"),
    ("PRJNA417304", "rejected", "ENA and SRA return no runs for this accession, so nothing could be read; the search's record says small RNA only", "not verifiable", "not applicable"),
]


def size_gb(path: pathlib.Path) -> float:
    return round(int(subprocess.run(["du", "-sb", str(path)], capture_output=True, text=True, check=True).stdout.split()[0]) / 1e9, 3)


units = pd.read_csv(sys.argv[1], sep="\t", keep_default_na=False).set_index("unit")
outcomes, studies = [], []
for dataset in sorted(DEPOSIT_DIR):
    folder = FINAL / dataset
    row = {"dataset": dataset, "category": "cooccurrence", "source_files": str(DATA / DEPOSIT_DIR[dataset] / "deposit"),
           "size_gb": size_gb(DATA / DEPOSIT_DIR[dataset]), "final_dir": str(folder)}
    if dataset not in units.index or units.status[dataset] != "ok":
        reason = units.note[dataset] if dataset in units.index else "not run"
        outcomes.append({**row, "cells_in_store": 0, "comparisons": 0, "fitted": 0, "not_fitted": 0, "failed": 1, "status": "failed", "reason": reason})
        studies.append({"dataset": dataset, "status": "failed", "reason": reason})
        continue
    info = json.loads((folder / "study.json").read_text())
    groups = pd.read_csv(folder / "correlation/groups.tsv", sep="\t", keep_default_na=False)
    correlated = groups[groups.status == "correlated"]
    outcomes.append({**row, "cells_in_store": info["cells_listed"], "comparisons": len(groups), "fitted": len(correlated),
                     "not_fitted": len(groups) - len(correlated), "failed": 0, "status": "done",
                     "reason": f"feature level {info['feature_level']}; miRNA normalisation {info['mirna_normalisation']}; allow_fractional {info['allow_fractional']}; {info['values']}; groups, not comparisons: fitted means correlated"})
    studies.append({"dataset": dataset, "status": "done", "same_cell": info["same_cell"], "feature_level": info["feature_level"],
                    "values": info["values"], "allow_fractional": info["allow_fractional"], "mirna_normalisation": info["mirna_normalisation"],
                    "library_total_source": info["library_total_source"], "cells_with_zero_mirna": int(groups.n_cells_zero_mirna.sum()), "cells_left_out": int(groups.n_cells_without_counts.sum()),
                    "mirna_integer": info["mirna_values_integer"], "mrna_integer": info["mrna_values_integer"],
                    "grouping": info["grouping"], "qc_source": info["qc_source"], "source": info["source"], "cells_listed": info["cells_listed"],
                    "cells_with_counts": int(groups.n_cells.sum()), "groups": len(groups), "groups_correlated": len(correlated),
                    "median_expressed_mirnas_correlated": correlated.n_expressed_mirnas.median() if len(correlated) else 0,
                    "max_expressed_mirnas": int(groups.n_expressed_mirnas.max()), "caveats": info["caveats"]})
for dataset, status, reason, same_cell, level in NOT_BUILT:
    outcome = "blocked" if status == "blocked" else "done"
    outcomes.append({"dataset": dataset, "category": "cooccurrence", "source_files": "", "size_gb": 0, "cells_in_store": 0, "comparisons": 0,
                     "fitted": 0, "not_fitted": 0, "failed": 0, "final_dir": "", "status": outcome, "reason": f"{status}: {reason}"})
    studies.append({"dataset": dataset, "status": status, "same_cell": same_cell, "feature_level": level, "reason": reason})
columns = ["dataset", "category", "source_files", "size_gb", "cells_in_store", "comparisons", "fitted", "not_fitted", "failed", "final_dir", "status", "reason"]
pd.DataFrame(outcomes)[columns].to_csv(OUTCOMES, sep="\t", index=False)
pd.DataFrame(studies).to_csv(FINAL / "studies.tsv", sep="\t", index=False)
print(pd.DataFrame(outcomes)[["dataset", "cells_in_store", "comparisons", "fitted", "not_fitted", "failed", "status"]].to_string(index=False))
