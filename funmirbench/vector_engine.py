"""Batched rank-vector calculations and compact caches for exploratory analysis."""
from __future__ import annotations

import json
import logging
import pathlib
import shutil
import tempfile
import time
import zipfile
from collections import OrderedDict

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from funmirbench.response_vectors import build_vectors, specificity, unit_rank_vector
from funmirbench.vector_sources import atomic_json, digest, fingerprint, prediction_stream

logger = logging.getLogger(__name__)
ENGINE_VERSION = "rank-vectors-1"


def load_matrix(source, direction, mirnas, genes, cache):
    """Read only three columns, filter each chunk and aggregate directly into a small matrix."""
    identity = {"version": ENGINE_VERSION, "source": fingerprint(source["path"]),
                "member": source.get("member"), "direction": direction, "mirnas": list(mirnas), "genes": list(genes)}
    path = cache / "matrices" / (digest(identity) + ".npz")
    if path.is_file():
        try:
            with np.load(path, allow_pickle=False) as saved:
                matrix = saved["scores"]
                if matrix.shape != (len(mirnas), len(genes)):
                    raise ValueError("invalid cached shape")
            logger.info("cached prediction matrix: %s", pathlib.Path(source["path"]).name)
            return matrix
        except (ValueError, OSError, KeyError, EOFError, zipfile.BadZipFile):
            logger.warning("rebuilding unreadable prediction cache %s", path)
    rows = {name: i for i, name in enumerate(mirnas)}
    columns = {name: i for i, name in enumerate(genes)}
    matrix = np.full((len(mirnas), len(genes)), -np.inf, dtype=np.float64)
    sign = -1 if direction == "lower_is_stronger" else 1
    scanned, retained = 0, 0
    with prediction_stream(source) as handle:
        for chunk in pd.read_csv(handle, sep="\t", usecols=["Ensembl_ID", "miRNA_Name", "Score"],
                                 dtype={"Ensembl_ID": str, "miRNA_Name": str}, chunksize=250_000):
            scanned += len(chunk)
            chunk = chunk[chunk.miRNA_Name.isin(rows)]
            if len(chunk):
                gene = chunk.Ensembl_ID.str.strip().str.split(".", n=1, regex=False).str[0].map(columns)
                value = sign * pd.to_numeric(chunk.Score, errors="coerce")
                good = gene.notna() & np.isfinite(value)
                retained += int(good.sum())
                np.maximum.at(matrix, (chunk.miRNA_Name[good].map(rows).to_numpy(int), gene[good].to_numpy(int)),
                              value[good].to_numpy(float))
            if scanned % 2_000_000 == 0:
                logger.info("read %s rows; retained %s", f"{scanned:,}", f"{retained:,}")
    matrix[~np.isfinite(matrix)] = np.nan
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as handle:
        tmp = pathlib.Path(handle.name)
        np.savez_compressed(handle, scores=matrix)
    tmp.replace(path)
    logger.info("parsed %s rows into %s miRNAs x %s genes", f"{scanned:,}", len(mirnas), len(genes))
    return matrix


class PredictionPanel:
    """Cache ranked predictions for repeated gene masks, with a bounded memory budget."""
    def __init__(self, matrix, max_cache_bytes=64 << 20):
        self.matrix = matrix
        self.cache = OrderedDict()
        self.cache_bytes = 0
        self.max_cache_bytes = max_cache_bytes

    def ranked(self, mask):
        key = np.packbits(mask).tobytes()
        if key in self.cache:
            self.cache.move_to_end(key)
            return self.cache[key]
        values = self.matrix[:, mask]
        if values.shape[1] < 3:
            result = np.full(values.shape, np.nan)
        else:
            # Missing scores form a bottom tie. No arbitrary numeric floor is needed.
            values = np.where(np.isfinite(values), values, -np.inf)
            ranks = rankdata(values, axis=1, method="average")
            ranks -= ranks.mean(axis=1, keepdims=True)
            norms = np.linalg.norm(ranks, axis=1, keepdims=True)
            result = np.divide(ranks, norms, out=np.full_like(ranks, np.nan), where=norms > 0)
        while self.cache and self.cache_bytes + result.nbytes > self.max_cache_bytes:
            _, old = self.cache.popitem(last=False)
            self.cache_bytes -= old.nbytes
        if result.nbytes <= self.max_cache_bytes:
            self.cache[key] = result
            self.cache_bytes += result.nbytes
        return result


def response_view(table):
    response = table.oriented_logFC_shrunk.to_numpy(float)
    usable = table.eligible.to_numpy(bool) & np.isfinite(response)
    groups = []
    bins = table.expression_tenth.to_numpy(float)
    for b in np.unique(bins[np.isfinite(bins)]):
        mask = usable & (bins == b)
        unit = unit_rank_vector(response[mask])
        groups.append((int(b), mask, unit))
    return usable, unit_rank_vector(response[usable]), groups


def batch_alignment(panel, view):
    usable, observed, groups = view
    n = len(panel.matrix)
    # These small matrix-vector products are memory-bound. Avoid starting a large
    # BLAS thread pool on shared servers for each experiment/expression bin.
    plain = (np.einsum("ij,j->i", panel.ranked(usable), observed, optimize=False)
             if usable.sum() >= 3 else np.full(n, np.nan))
    values = np.full((n, len(groups)), np.nan)
    for i, (_, mask, unit) in enumerate(groups):
        if mask.sum() >= 20 and np.isfinite(unit).all():
            values[:, i] = np.einsum("ij,j->i", panel.ranked(mask), unit, optimize=False)
    weights = np.array([mask.sum() for _, mask, _ in groups], dtype=float)
    finite = np.isfinite(values)
    den = finite @ weights
    within = np.divide(np.nansum(values * weights, axis=1), den,
                       out=np.full(n, np.nan), where=den > 0)
    support = [",".join(str(groups[i][0]) for i in np.flatnonzero(row)) for row in finite]
    return plain, within, support, den.astype(int), values


def analyse(tables, comps, tools, mirnas, sources, cache, *, input_files, export_vectors=False):
    genes = next(iter(tables.values())).index
    names = mirnas.index
    signature = {"engine": ENGINE_VERSION, "inputs": [fingerprint(p) for p in input_files],
                 "tools": tools.to_dict("records"), "mirnas": mirnas.to_dict("index"),
                 "comparisons": comps.to_dict("records"), "export_vectors": export_vectors,
                 "sources": {k: {**v, "file": fingerprint(v["path"])} for k, v in sources.items()}}
    key = digest(signature)
    destination = cache / "analyses" / key
    required = ("candidate_alignment.tsv", "specificity.tsv", "expression_bins.tsv", "analysis.json", "vector_files.tsv")
    if all((destination / name).is_file() for name in required):
        try:
            saved = json.loads((destination / "analysis.json").read_text())
            intact = bool(saved.get("outputs")) and all(
                (destination / name).is_file() and
                {k: v for k, v in fingerprint(destination / name).items() if k != "path"} == meta
                for name, meta in saved.get("outputs", {}).items())
            if intact:
                logger.info("cached analysis: %s", key[:12])
                return destination, key, True
        except (ValueError, OSError):
            logger.warning("rebuilding unreadable analysis cache")
    destination.parent.mkdir(parents=True, exist_ok=True)
    start = time.perf_counter()
    views = {exp: response_view(t) for exp, t in tables.items()}
    seeds = mirnas.mirna_sequence.str.strip().str.upper().str.replace("T", "U", regex=False).str[1:8]
    candidate_rows, summary_rows, bin_rows, files = [], [], [], []
    with tempfile.TemporaryDirectory(dir=destination.parent) as tmp:
        tmp = pathlib.Path(tmp)
        for tool_number, tool in enumerate(tools.itertuples(index=False)):
            matrix = load_matrix(sources[tool.tool_id], tool.score_direction, names, genes, cache)
            panel = PredictionPanel(matrix)
            n_scored = np.isfinite(matrix).sum(axis=1)
            for exp_number, c in enumerate(comps.itertuples(index=False)):
                usable, observed, groups = views[c.experiment_id]
                plain, within, support, n_within, bins = batch_alignment(panel, views[c.experiment_id])
                n_scored_usable = np.isfinite(matrix[:, usable]).sum(axis=1)
                rows = []
                for i, candidate in enumerate(names):
                    reason = ("fewer_than_three_usable_genes" if usable.sum() < 3 else
                              "constant_response" if not np.isfinite(observed).all() else
                              "constant_prediction" if not np.isfinite(plain[i]) else "")
                    rows.append({"experiment_id": c.experiment_id, "predictor": tool.tool_id,
                                 "candidate_mirna": candidate, "candidate_seed": seeds[candidate],
                                 "is_cognate": candidate == c.mirna_name, "n_genes": len(genes),
                                 "n_usable": int(usable.sum()), "n_scored": int(n_scored[i]),
                                 "n_scored_usable": int(n_scored_usable[i]), "spearman": plain[i],
                                 "spearman_within_tenths": within[i], "within_tenths_bins": support[i],
                                 "n_within_tenths_genes": n_within[i], "alignment_na_reason": reason})
                candidates = pd.DataFrame(rows)
                candidate_rows.extend(rows)
                own = candidates[candidates.is_cognate].iloc[0]
                cognate = names.get_loc(c.mirna_name)
                for metric in ("spearman", "spearman_within_tenths"):
                    summary_rows.append({"experiment_id": c.experiment_id, "predictor": tool.tool_id,
                                         "mirna_name": c.mirna_name, "action": c.action,
                                         "dependency_group": c.dependency_group, "metric": metric,
                                         "alignment": own[metric], "n_genes": len(genes), "n_usable": own.n_usable,
                                         "n_scored_usable": own.n_scored_usable,
                                         "coverage": own.n_scored_usable / own.n_usable if own.n_usable else np.nan,
                                         "alignment_na_reason": own.alignment_na_reason if metric == "spearman" else
                                         ("" if np.isfinite(own[metric]) else "no_usable_expression_bins"),
                                         **specificity(candidates, c.mirna_name, metric)})
                for b, (bin_number, mask, unit) in enumerate(groups):
                    value = bins[cognate, b]
                    reason = ("fewer_than_20_usable_genes" if mask.sum() < 20 else
                              "constant_response" if not np.isfinite(unit).all() else
                              "constant_prediction" if not np.isfinite(value) else "")
                    bin_rows.append({"experiment_id": c.experiment_id, "predictor": tool.tool_id,
                                     "expression_tenth": bin_number, "n_usable": int(mask.sum()),
                                     "spearman": value, "na_reason": reason})
                if export_vectors:
                    scores = pd.Series(matrix[cognate], index=genes)
                    vectors = build_vectors(tables[c.experiment_id], scores)
                    vectors["expression_tenth"] = tables[c.experiment_id].expression_tenth
                    filename = f"vectors/{exp_number:04d}_{tool_number:03d}.tsv.gz"
                    (tmp / "vectors").mkdir(exist_ok=True)
                    vectors.to_csv(tmp / filename, sep="\t", na_rep="NA")
                    files.append({"experiment_id": c.experiment_id, "predictor": tool.tool_id, "file": filename})
            logger.info("analysed %s across %d experiments", tool.tool_id, len(tables))
        pd.DataFrame(candidate_rows).to_csv(tmp / "candidate_alignment.tsv", sep="\t", index=False, na_rep="NA")
        pd.DataFrame(summary_rows).to_csv(tmp / "specificity.tsv", sep="\t", index=False, na_rep="NA")
        pd.DataFrame(bin_rows, columns=["experiment_id", "predictor", "expression_tenth", "n_usable", "spearman",
                                       "na_reason"]).to_csv(tmp / "expression_bins.tsv", sep="\t", index=False, na_rep="NA")
        pd.DataFrame(files, columns=["experiment_id", "predictor", "file"]).to_csv(
            tmp / "vector_files.tsv", sep="\t", index=False)
        outputs = {str(p.relative_to(tmp)): {k: v for k, v in fingerprint(p).items() if k != "path"}
                   for p in tmp.rglob("*") if p.is_file()}
        atomic_json(tmp / "analysis.json", {"signature": signature, "key": key, "outputs": outputs,
                                           "analysis_seconds": time.perf_counter() - start})
        if destination.exists():
            shutil.rmtree(destination)
        shutil.move(str(tmp), str(destination))
    return destination, key, False
