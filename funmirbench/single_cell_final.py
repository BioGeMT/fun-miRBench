"""Single-cell final data: cell store and comparison list to per-comparison tables.

For each comparison the passing cells are summed per replicate unit and arm within each batch both arms share, the
sums are fitted with the bulk DESeq2 script and labelled with the bulk label_states.R on the benchmark gene list. A
comparison that cannot be fitted gets a status, a reason and the log2 ratio of the sums. All sums of arms are computed
once per store and shared between comparisons, so the control cells of a screen are read and summed once.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import os
import pathlib
import shutil
import subprocess
import sys
import traceback
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from scipy import sparse

from funmirbench.experiments_pipeline import repo_root, run_deseq2_from_counts
from funmirbench.logger import parse_log_level, setup_logging
from funmirbench.single_cell_store import ARM_SEPARATOR, read_cell_store

logger = logging.getLogger(__name__)

LABEL_SCRIPT = pathlib.Path("pipelines") / "processing" / "src" / "label_states.R"
LIST_COLUMNS = ("comparison_id", "dataset", "category", "mirna_or_gene", "orientation", "treated_arm", "control_arm",
                "min_cells", "caveats")
OPTIONAL_LIST_COLUMNS = ("other_pathway_arms", "target_gene_id")
CATEGORIES = ("raise", "lower", "pooled_mirna", "pathway", "host_gene", "reporter", "viral_mirna", "ko", "other")
ORIENTATIONS = ("OE", "KO")
HIGH_GUIDE_LOAD = 1.5
MIN_UNITS_PER_SIDE = 2
STATES = ("responsive", "no_response_above_margin", "unresolved")
# fit_failed is a result (failures stay listed), so run.sh must not mark the run FAILED for it; error is infrastructure
UNIT_STATUS = {"fitted": "ok", "not_fitted": "skipped", "fit_failed": "skipped", "error": "failed"}


def read_comparison_list(path: pathlib.Path, dataset: str) -> pd.DataFrame:
    table = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    missing = [column for column in LIST_COLUMNS if column not in table]
    if missing:
        raise ValueError(f"{path}: missing columns {missing}")
    for column in OPTIONAL_LIST_COLUMNS:
        if column not in table:
            table[column] = ""
    if table.empty or table.comparison_id.duplicated().any():
        raise ValueError(f"{path}: needs at least one row and unique comparison_id values")
    for column in ("comparison_id", "mirna_or_gene", "treated_arm", "control_arm"):
        if (table[column] == "").any():
            raise ValueError(f"{path}: empty {column}")
    if not table.category.isin(CATEGORIES).all():
        raise ValueError(f"{path}: category must be one of {CATEGORIES} (co-occurrence datasets use "
                         f"funmirbench-single-cell-cooccurrence), found {sorted(set(table.category) - set(CATEGORIES))}")
    if not table.orientation.isin(ORIENTATIONS).all():
        raise ValueError(f"{path}: orientation must be one of {ORIENTATIONS}")
    if not (table.dataset == dataset).all():
        raise ValueError(f"{path}: every row must have dataset {dataset} (the cell store's)")
    if (table.treated_arm == table.control_arm).any():
        raise ValueError(f"{path}: treated_arm and control_arm must differ")
    if not table.min_cells.str.fullmatch(r"[1-9]\d*").all():
        raise ValueError(f"{path}: min_cells must be a positive integer")
    return table.assign(min_cells=table.min_cells.astype(int))


def group_sums(counts: sparse.csr_matrix, obs: pd.DataFrame, keys: list[str]) -> tuple[pd.DataFrame, sparse.csr_matrix]:
    """Per group of keys: a table with n_cells and replicate_type, and the summed counts (groups x genes)."""
    grouped = obs.groupby(keys, sort=True)
    table = grouped.agg(n_cells=("arm", "size"), replicate_type=("replicate_type", lambda s: ";".join(sorted(set(s)))),
                        control_type=("control_type", lambda s: ";".join(sorted(set(s) - {""})))).reset_index()
    codes = grouped.ngroup().to_numpy()
    # int64 membership so the product of integer counts is int64: a gene's sum over a whole screen arm can pass 2^31
    membership = sparse.csr_matrix((np.ones(len(codes), dtype=np.int64), (codes, np.arange(len(codes)))),
                                   shape=(len(table), len(codes)))
    return table, (membership @ counts).tocsr()


@dataclass
class Units:
    """One side of a comparison: one row per replicate unit (batch, replicate_unit, n_cells, ...) and its sums."""
    table: pd.DataFrame
    sums: sparse.csr_matrix

    def where(self, mask) -> "Units":
        mask = np.asarray(mask, dtype=bool)
        return Units(self.table[mask].reset_index(drop=True), self.sums[np.flatnonzero(mask)])


class CellSums:
    """The passing cells of a cell store and their sums, shared read-only by all comparisons."""

    def __init__(self, store):
        passing = store[store.obs.passed_qc.to_numpy()]
        has_id = passing.var.gene_id.astype(str).to_numpy() != ""
        self.gene_ids = passing.var.gene_id.astype(str).to_numpy()[has_id]
        self.counts = sparse.csr_matrix(passing.X)[:, np.flatnonzero(has_id)]
        self.obs = pd.DataFrame({column: passing.obs[column].astype(str).to_numpy() for column in
                                 ("batch", "replicate_unit", "replicate_type", "arm", "control_type")})
        self.source = str(store.uns["source"])
        self.fractional = bool(store.uns.get("fractional_counts", False))
        self.median_guides = float(passing.obs.n_guides.median()) if passing.n_obs else 0.0
        self.high_guide_load = self.median_guides > HIGH_GUIDE_LOAD
        self._arm_table, self._arm_sums = group_sums(self.counts, self.obs, ["batch", "replicate_unit", "arm"])
        if self.high_guide_load:
            self._unit_table, self._unit_sums = group_sums(self.counts, self.obs, ["batch", "replicate_unit"])
            targets = self.obs.arm.str.split(ARM_SEPARATOR).explode()
            self._carriers = {target: cells.index.to_numpy() for target, cells in targets.groupby(targets)}

    def _carrying(self, targets: list[str]) -> np.ndarray:
        found = [self._carriers[target] for target in targets if target in self._carriers]
        return np.unique(np.concatenate(found)) if found else np.array([], dtype=np.int64)

    def _units_of_cells(self, cells: np.ndarray) -> Units:
        if not len(cells):
            return Units(self._unit_table.iloc[:0], self._unit_sums[:0])
        return Units(*group_sums(self.counts[cells], self.obs.iloc[cells], ["batch", "replicate_unit"]))

    def _units_without_cells(self, cells: np.ndarray) -> Units:
        """Every replicate unit minus the given cells: the unit sums are computed once, only the few cells carrying
        the comparison's guides are summed and subtracted."""
        removed = self._units_of_cells(cells)
        rows = self._unit_table.reset_index().merge(removed.table[["batch", "replicate_unit"]].reset_index(),
                                                    on=["batch", "replicate_unit"], suffixes=("", "_removed"))
        n_cells = self._unit_table.n_cells.to_numpy().copy()
        n_cells[rows["index"].to_numpy()] -= removed.table.n_cells.to_numpy()[rows.index_removed.to_numpy()]
        position = sparse.csr_matrix((np.ones(len(rows), dtype=np.int64),
                                      (rows["index"].to_numpy(), rows.index_removed.to_numpy())),
                                     shape=(len(self._unit_table), len(removed.table)))
        return Units(self._unit_table.assign(n_cells=n_cells), (self._unit_sums - position @ removed.sums).tocsr())

    def _whole_counts(self, units: Units) -> Units:
        """Fractional counts are rounded here, after summing: the fit takes integers."""
        if not self.fractional:
            return units
        sums = units.sums.copy()
        sums.data = np.round(sums.data)
        sums.eliminate_zeros()
        return Units(units.table, sums.astype(np.int64))

    def sides(self, comparison) -> tuple[Units, Units]:
        """Treated and control units of one comparison, before any min_cells or batch rule."""
        if not self.high_guide_load:
            arm_units = Units(self._arm_table, self._arm_sums)
            return (self._whole_counts(arm_units.where(self._arm_table.arm == comparison.treated_arm)),
                    self._whole_counts(arm_units.where(self._arm_table.arm == comparison.control_arm)))
        others = [arm for arm in comparison.other_pathway_arms.split(ARM_SEPARATOR) if arm]
        return (self._whole_counts(self._units_of_cells(self._carrying([comparison.treated_arm]))),
                self._whole_counts(self._units_without_cells(self._carrying([comparison.treated_arm, *others]))))


def select_units(treated: Units, control: Units, min_cells: int) -> tuple[Units, Units, str]:
    """Units with at least min_cells cells in the batches both arms share, and why they cannot be fitted ("" if they
    can)."""
    if treated.table.n_cells.sum() == 0 or control.table.n_cells.sum() == 0:
        side = "treated" if treated.table.n_cells.sum() == 0 else "control"
        return treated.where(np.zeros(len(treated.table))), control.where(np.zeros(len(control.table))), \
            f"the {side} arm has no passing cells"
    treated_enough = treated.table.n_cells >= min_cells
    control_enough = control.table.n_cells >= min_cells
    shared = set(treated.table.batch[treated_enough]) & set(control.table.batch[control_enough])
    kept_treated = treated.where(treated_enough & treated.table.batch.isin(shared))
    kept_control = control.where(control_enough & control.table.batch.isin(shared))
    if not shared:
        return kept_treated, kept_control, f"no batch holds both arms in units of at least {min_cells} cells"
    n_treated, n_control = len(kept_treated.table), len(kept_control.table)
    if min(n_treated, n_control) < MIN_UNITS_PER_SIDE:
        return kept_treated, kept_control, (f"{n_treated} treated and {n_control} control units of at least "
                                            f"{min_cells} cells in shared batches; {MIN_UNITS_PER_SIDE} per side needed")
    return kept_treated, kept_control, ""


def write_pseudobulk(directory: pathlib.Path, gene_ids: np.ndarray, treated: Units, control: Units) -> pd.DataFrame:
    """pseudobulk_counts.tsv.gz and design.tsv of one comparison. Columns are named <condition>_<i> because the DESeq2
    script takes them as a comma-separated list; design.tsv gives the batch and replicate unit of each."""
    design, columns = [], {}
    for condition, units in (("control", control), ("treated", treated)):
        dense = units.sums.toarray()
        for i, row in enumerate(units.table.itertuples(index=False)):
            column = f"{condition}_{i + 1}"
            columns[column] = dense[i]
            design.append({"column": column, "condition": condition, "batch": row.batch,
                           "replicate_unit": row.replicate_unit, "replicate_type": row.replicate_type,
                           "n_cells": row.n_cells})
    pd.DataFrame(columns, index=pd.Index(gene_ids, name="gene_id")).to_csv(directory / "pseudobulk_counts.tsv.gz",
                                                                            sep="\t")
    design = pd.DataFrame(design, columns=["column", "condition", "batch", "replicate_unit", "replicate_type",
                                           "n_cells"])
    design.to_csv(directory / "design.tsv", sep="\t", index=False)
    return design


def ratio_table(genes: pd.DataFrame, gene_ids: np.ndarray, treated: Units, control: Units) -> pd.DataFrame:
    """Per listed gene, the summed counts and the library-normalised log2 ratio of the sums (0.5 added to each sum).
    No ratio for a gene with no count on either side or not in the store."""
    treated_sum = np.asarray(treated.sums.sum(axis=0)).ravel()
    control_sum = np.asarray(control.sums.sum(axis=0)).ravel()
    if treated_sum.sum() > 0 and control_sum.sum() > 0:
        ratio = np.log2((treated_sum + 0.5) / (treated_sum.sum() + 1)) \
            - np.log2((control_sum + 0.5) / (control_sum.sum() + 1))
        ratio[(treated_sum + control_sum) == 0] = np.nan
    else:
        ratio = np.full(len(gene_ids), np.nan)
    by_gene = pd.DataFrame({"treated_sum": treated_sum, "control_sum": control_sum, "log2_ratio_of_sums": ratio},
                           index=gene_ids).reindex(genes.gene_id)
    return pd.DataFrame({"gene_id": genes.gene_id.to_numpy(), "gene_name": genes.gene_name.to_numpy(),
                         "treated_sum": by_gene.treated_sum.astype("Int64").to_numpy(),
                         "control_sum": by_gene.control_sum.astype("Int64").to_numpy(),
                         "log2_ratio_of_sums": by_gene.log2_ratio_of_sums.to_numpy(),
                         "analysis_status": "not_fitted"})


def last_lines(path: pathlib.Path, n: int = 2, width: int = 400) -> str:
    lines = [line.strip() for line in path.read_text(errors="replace").splitlines() if line.strip()]
    return " | ".join(lines[-n:])[-width:]


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def fit_and_label(directory: pathlib.Path, design: pd.DataFrame, comparison, *, gene_list: pathlib.Path,
                  n_genes: int, out: pathlib.Path, repo: pathlib.Path) -> tuple[str, dict]:
    """DESeq2 and label_states.R on one comparison's pseudobulk. Returns the reason the fit failed ("" if it
    worked) and the final table's summary. Rscript killed by a signal raises: that is the machine, not the data."""
    covariates = None
    if design.batch.nunique() > 1:
        covariates = directory / "sample_covariates.tsv"
        design.rename(columns={"column": "sample_id"})[["sample_id", "batch"]].to_csv(covariates, sep="\t", index=False)
    try:
        run_deseq2_from_counts(repo=repo, run_dir=directory, counts_path=directory / "pseudobulk_counts.tsv.gz",
                               gene_id_column="gene_id", control_cols=design.column[design.condition == "control"].tolist(),
                               treated_cols=design.column[design.condition == "treated"].tolist(),
                               output_path=directory / "de_results.tsv", sample_covariates_path=covariates)
    except RuntimeError as error:
        if "exit code -" in str(error):
            raise
        return f"DESeq2: {last_lines(directory / 'deseq2.stderr.txt')}", {}

    prefix = out / "genes" / comparison.comparison_id
    log = directory / "label_states.log"
    with open(log, "w") as handle:
        done = subprocess.run(["Rscript", str(repo / LABEL_SCRIPT), str(directory / "complete_analysis.rds"),
                               str(gene_list), str(prefix), comparison.orientation],
                              stdout=handle, stderr=subprocess.STDOUT, check=False)
    if done.returncode < 0:
        raise RuntimeError(f"label_states.R killed by signal {-done.returncode}")
    if done.returncode > 0:
        return f"label_states.R: {last_lines(log)}", {}
    for suffix in ("_events.tsv", "_metadata.txt"):
        pathlib.Path(f"{prefix}{suffix}").replace(directory / f"label_states{suffix}")
    table_path = pathlib.Path(f"{prefix}.tsv.gz")
    table = pd.read_csv(table_path, sep="\t", usecols=["eligible", "analysis_status", "state_d05"])
    if len(table) != n_genes:
        raise RuntimeError(f"{table_path}: {len(table)} rows, expected {n_genes}")
    summary = {"n_genes_fitted": int((table.analysis_status != "not_fitted_all_zero").sum()),
               "n_genes_eligible": int(table.eligible.sum()), "table_sha256": sha256(table_path)}
    states = table.state_d05.value_counts()
    summary.update({f"d05_{state}": int(states.get(state, 0)) for state in STATES})
    return "", summary


def target_change(comparison, gene_ids: np.ndarray, treated: Units, control: Units, directory: pathlib.Path,
                  table_path: pathlib.Path) -> dict:
    """The target's own change, as information only: its summed counts, the log2 ratio of the sums, and, when fitted,
    its DESeq2 log2FC and shrunken log2FC. target_gene_id may list several IDs of one symbol; the first in the store
    counts."""
    wanted = [gene for gene in comparison.target_gene_id.split(ARM_SEPARATOR) if gene]
    where = {gene: position for position, gene in enumerate(gene_ids)}
    found = next((gene for gene in wanted if gene in where), "")
    if not found:
        return {"target_in_store": False}
    treated_sum = int(treated.sums[:, where[found]].sum())
    control_sum = int(control.sums[:, where[found]].sum())
    treated_total, control_total = int(treated.sums.sum()), int(control.sums.sum())
    change = {"target_in_store": True, "target_used_gene_id": found, "target_treated_sum": treated_sum,
              "target_control_sum": control_sum}
    if treated_total > 0 and control_total > 0 and treated_sum + control_sum > 0:
        change["target_log2_ratio_of_sums"] = float(np.log2((treated_sum + 0.5) / (treated_total + 1))
                                                    - np.log2((control_sum + 0.5) / (control_total + 1)))
    results = directory / "de_results.tsv"
    if results.exists():
        fit = pd.read_csv(results, sep="\t", usecols=["gene_id", "logFC"]).set_index("gene_id").logFC
        change["target_log2fc"] = fit.get(found, np.nan)
    if table_path.exists():
        table = pd.read_csv(table_path, sep="\t")
        if "logFC_shrunk" in table:
            shrunk = table.set_index("gene_id").logFC_shrunk
            change["target_log2fc_shrunk"] = shrunk.get(found, np.nan)
    return change


def run_comparison(comparison, cell_sums: CellSums, *, genes: pd.DataFrame, gene_list: pathlib.Path,
                   out: pathlib.Path, repo: pathlib.Path) -> dict:
    record = {column: getattr(comparison, column) for column in (*LIST_COLUMNS, *OPTIONAL_LIST_COLUMNS)}
    record["source"] = cell_sums.source
    record["median_guides_per_cell"] = cell_sums.median_guides
    record["control_selection"] = "arm"
    if cell_sums.high_guide_load:
        record["control_selection"] = "cells_of_the_batch_without_these_guides"
        note = (f"high guide load: median {cell_sums.median_guides:g} guides per cell; control is the cells of the "
                f"batch with no guide of {comparison.treated_arm} or of the arms in other_pathway_arms")
        record["caveats"] = "; ".join(part for part in (comparison.caveats, note) if part)
    directory = out / "comparisons" / comparison.comparison_id
    directory.mkdir(parents=True)

    all_treated, all_control = cell_sums.sides(comparison)
    treated, control, reason = select_units(all_treated, all_control, comparison.min_cells)
    # not fitted with nothing kept on a side: the ratio still uses every passing cell of the two arms
    if len(treated.table) == 0 or len(control.table) == 0:
        treated, control = all_treated.where(all_treated.table.n_cells > 0), all_control.where(all_control.table.n_cells > 0)
    design = write_pseudobulk(directory, cell_sums.gene_ids, treated, control)
    record.update(
        n_treated_cells=int(treated.table.n_cells.sum()), n_control_cells=int(control.table.n_cells.sum()),
        n_treated_units=len(treated.table), n_control_units=len(control.table),
        n_treated_units_below_min_cells=int((all_treated.table.n_cells < comparison.min_cells).sum()),
        n_control_units_below_min_cells=int((all_control.table.n_cells < comparison.min_cells).sum()),
        n_batches=design.batch.nunique(),
        replicate_type=";".join(sorted({part for value in design.replicate_type for part in value.split(";")})),
        control_type=";".join(sorted({part for value in control.table.control_type for part in value.split(";") if part})),
    )

    status, summary = ("not_fitted", {}) if reason else ("fitted", {})
    if not reason:
        record["design"] = "~ batch + condition" if design.batch.nunique() > 1 else "~ condition"
        reason, summary = fit_and_label(directory, design, comparison, gene_list=gene_list, n_genes=len(genes),
                                        out=out, repo=repo)
        status = "fit_failed" if reason else "fitted"
    if status != "fitted":
        table_path = out / "genes" / f"{comparison.comparison_id}.tsv.gz"
        ratio_table(genes, cell_sums.gene_ids, treated, control).to_csv(table_path, sep="\t", index=False, na_rep="NA")
        summary = {"table_sha256": sha256(table_path)}
    record.update(status=status, reason=reason, table=f"genes/{comparison.comparison_id}.tsv.gz", **summary)
    if comparison.target_gene_id:
        record.update(target_change(comparison, cell_sums.gene_ids, treated, control, directory,
                                    out / "genes" / f"{comparison.comparison_id}.tsv.gz"))
    return record


def preflight(repo: pathlib.Path, genes: pd.DataFrame) -> None:
    if not {"gene_id", "gene_name"} <= set(genes):
        raise ValueError("--genes needs gene_id and gene_name")
    if not (repo / LABEL_SCRIPT).exists():
        raise FileNotFoundError(repo / LABEL_SCRIPT)
    if shutil.which("Rscript") is None:
        raise FileNotFoundError("Rscript is not on PATH; activate the environment with DESeq2 and apeglm")
    check = subprocess.run(["Rscript", "-e", "suppressPackageStartupMessages({library(DESeq2); library(apeglm)})"],
                           capture_output=True, text=True, check=False)
    if check.returncode != 0:
        raise RuntimeError(f"Rscript cannot load DESeq2 and apeglm: {check.stderr.strip()[-400:]}")


def run_final(store_path: pathlib.Path, comparisons_path: pathlib.Path, gene_list: pathlib.Path, out: pathlib.Path,
              *, workers: int, units_path: pathlib.Path | None = None, repo: pathlib.Path | None = None) -> pd.DataFrame:
    repo = (repo or repo_root()).resolve()
    genes = pd.read_csv(gene_list, sep="\t", dtype=str)
    preflight(repo, genes)
    store = read_cell_store(store_path)
    dataset = str(store.obs.dataset.iloc[0])
    comparisons = read_comparison_list(comparisons_path, dataset)
    (out / "genes").mkdir(parents=True, exist_ok=False)

    cell_sums = CellSums(store)
    del store
    print(f"{dataset}: {cell_sums.counts.shape[0]} passing cells, {len(cell_sums.gene_ids)} genes with a gene_id, "
          f"median guides per cell {cell_sums.median_guides:g}", flush=True)

    def run_one(comparison) -> dict:
        try:
            return run_comparison(comparison, cell_sums, genes=genes, gene_list=gene_list, out=out, repo=repo)
        except Exception as error:  # one comparison's machine failure must not stop the others; main exits non-zero
            logger.error("%s\n%s", comparison.comparison_id, traceback.format_exc())
            return {**{column: getattr(comparison, column) for column in (*LIST_COLUMNS, *OPTIONAL_LIST_COLUMNS)},
                    "source": cell_sums.source, "status": "error", "reason": f"{type(error).__name__}: {error}"[:400]}

    rows = {}
    with ThreadPoolExecutor(workers) as pool:
        futures = [pool.submit(run_one, comparison) for comparison in comparisons.itertuples(index=False)]
        for done, future in enumerate(as_completed(futures), start=1):
            row = future.result()
            rows[row["comparison_id"]] = row
            print(f"PROGRESS {done}/{len(futures)} {row['comparison_id']} {row['status']}", flush=True)
    table = pd.DataFrame([rows[comparison_id] for comparison_id in comparisons.comparison_id])
    table.to_csv(out / "comparisons.tsv", sep="\t", index=False)
    units = table.assign(unit_status=table.status.map(UNIT_STATUS))[["comparison_id", "unit_status", "reason"]]
    units.to_csv(units_path or out / "units.tsv", sep="\t", index=False, header=False)
    manifest = {"built": datetime.now(timezone.utc).isoformat(timespec="seconds"), "dataset": dataset,
                "cell_store": str(store_path.resolve()), "cell_store_bytes": store_path.stat().st_size,
                "comparison_list": str(comparisons_path.resolve()), "comparison_list_sha256": sha256(comparisons_path),
                "gene_list": str(gene_list.resolve()), "gene_list_sha256": sha256(gene_list),
                "label_states_R_sha256": sha256(repo / LABEL_SCRIPT), "workers": workers,
                "status_counts": {key: int(value) for key, value in table.status.value_counts().items()}}
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return table


def combine_comparisons(root: pathlib.Path) -> pd.DataFrame:
    """One table of the comparisons.tsv of every dataset folder under root; root/comparisons.tsv is replaced. A folder
    with a SUPERSEDED file (a replaced build, naming its replacement) stays out."""
    frames = []
    for path in sorted(root.glob("*/comparisons.tsv")):
        if (path.parent / "SUPERSEDED").exists():
            continue
        frames.append(pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False).assign(final_dir=path.parent.name))
    if not frames:
        raise ValueError(f"{root}: no <dataset>/comparisons.tsv")
    table = pd.concat(frames, ignore_index=True)
    if table.comparison_id.duplicated().any():
        raise ValueError(f"{root}: duplicated comparison_id, for example "
                         f"{table.comparison_id[table.comparison_id.duplicated()].head(3).tolist()}")
    temporary = root / "comparisons.tsv.tmp"
    table.to_csv(temporary, sep="\t", index=False)
    temporary.replace(root / "comparisons.tsv")
    return table


def combine_main() -> None:
    parser = argparse.ArgumentParser(description=combine_comparisons.__doc__.splitlines()[0])
    parser.add_argument("--root", type=pathlib.Path, required=True)
    table = combine_comparisons(parser.parse_args().root)
    print(table.groupby(["category", "status"]).size().to_string())


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--store", type=pathlib.Path, required=True)
    parser.add_argument("--comparisons", type=pathlib.Path, required=True)
    parser.add_argument("--genes", type=pathlib.Path, required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--units", type=pathlib.Path)
    parser.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"])
    args = parser.parse_args()
    setup_logging(parse_log_level(args.log_level))
    for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[variable] = "1"  # one R process per worker
    table = run_final(args.store, args.comparisons, args.genes, args.out, workers=args.workers, units_path=args.units)
    print(table.status.value_counts().to_string())
    failed = table[table.status != "fitted"]
    if len(failed):
        print("Comparisons not fitted:")
        print(failed[["comparison_id", "status", "reason"]].to_string(index=False, max_colwidth=120))
    if (table.status == "error").any():
        sys.exit(1)


if __name__ == "__main__":
    main()
