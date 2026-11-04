"""Inspect gene-level response/prediction vectors without changing benchmark labels.

Uses existing OE benchmark tables and the same predictor registry as funmirbench.
Prediction scores are rankings, not calibrated predictions of log2 fold changes.
"""
from __future__ import annotations

import argparse
import logging
import pathlib
import shutil

import numpy as np
import pandas as pd
import yaml
from scipy.stats import rankdata

from funmirbench.benchmark import binned_spearman, expression_tenths, load_predictions, spearman
from funmirbench.benchmark_config import (
    build_run_dir_name, filter_df, resolve_predictor_output_path, validate_predictor_output_files,
)
from funmirbench.logger import parse_log_level, setup_logging

logger = logging.getLogger(__name__)
RESPONSE_COLUMNS = ["gene_id", "eligible", "oriented_logFC_shrunk", "control_mean_normalized_count"]


def unit_rank_vector(values: np.ndarray) -> np.ndarray:
    """Centered average ranks with unit norm; undefined for fewer than 3 or constant values."""
    values = np.asarray(values, dtype=float)
    if len(values) < 3 or not np.isfinite(values).all():
        return np.full(len(values), np.nan)
    ranks = rankdata(values, method="average")
    ranks -= ranks.mean()
    norm = np.linalg.norm(ranks)
    return ranks / norm if norm > 0 else np.full(len(values), np.nan)


def build_vectors(table: pd.DataFrame, scores: pd.Series) -> pd.DataFrame:
    """Align to the table gene order, preserving missing responses and score coverage.

    An ordinal score puts all unscored genes in one bottom tie without inventing a
    numeric effect size. Ranks are calculated on the usable response genes only.
    """
    if table.index.has_duplicates or table.index.hasnans or len(table) == 0:
        raise ValueError("response table needs a nonempty, unique gene_id index without missing IDs")
    if scores.index.has_duplicates:
        raise ValueError("prediction scores need unique gene IDs")
    eligible = table.eligible
    if eligible.isna().any() or not eligible.isin([True, False]).all():
        raise ValueError("eligible must contain boolean values without missing values")
    response = pd.to_numeric(table.oriented_logFC_shrunk, errors="raise")
    usable = eligible.astype(bool) & np.isfinite(response)
    prediction = pd.to_numeric(scores.reindex(table.index), errors="raise")
    prediction = prediction.where(np.isfinite(prediction))
    scored = prediction.notna()
    # Finite scores have ranks >= 1; 0 denotes the unscored bottom tie.
    ordinal = prediction.rank(method="average").fillna(0.0)
    result = pd.DataFrame({
        "observed_response": response.where(usable), "eligible": eligible.astype(bool),
        "response_usable": usable, "control_mean_normalized_count": table.control_mean_normalized_count,
        "prediction_score": prediction, "predictor_scored": scored, "prediction_ordinal": ordinal,
    }, index=table.index)
    result.index.name = "gene_id"
    result["response_rank_unit"] = np.nan
    result["prediction_rank_unit"] = np.nan
    result.loc[usable, "response_rank_unit"] = unit_rank_vector(response[usable].to_numpy())
    result.loc[usable, "prediction_rank_unit"] = unit_rank_vector(ordinal[usable].to_numpy())
    return result


def alignment(vectors: pd.DataFrame, bins: pd.Series) -> dict:
    """Rank alignment on the common response mask, including unscored genes."""
    usable = vectors.response_usable
    x = vectors.loc[usable, "prediction_ordinal"].to_numpy(float)
    y = vectors.loc[usable, "observed_response"].to_numpy(float)
    if len(y) < 3:
        reason = "fewer_than_three_usable_genes"
    elif np.all(x == x[0]):
        reason = "constant_prediction"
    elif np.all(y == y[0]):
        reason = "constant_response"
    else:
        reason = ""
    return {
        "n_genes": len(vectors), "n_usable": int(usable.sum()),
        "n_scored": int(vectors.predictor_scored.sum()),
        "n_scored_usable": int((usable & vectors.predictor_scored).sum()),
        "spearman": spearman(x, y),
        "spearman_within_tenths": binned_spearman(x, y, bins.loc[usable].to_numpy()),
        "alignment_na_reason": reason,
    }


def specificity(candidates: pd.DataFrame, cognate: str, metric: str) -> dict:
    """Cognate alignment versus finite different-seed miRNAs, with midrank ties.

    This is a descriptive panel comparison, not a calibrated significance test or
    a search over all human seed families. Same-seed members are not decoys.
    """
    own = candidates.loc[candidates.candidate_mirna == cognate].iloc[0]
    decoys = candidates.loc[(candidates.candidate_seed != own.candidate_seed) &
                            np.isfinite(candidates[metric]), metric].to_numpy(float)
    value = own[metric]
    result = {"n_decoys": len(decoys), "decoy_mean": np.nan, "decoy_sd": np.nan,
              "cognate_rank": np.nan, "fraction_decoys_better": np.nan, "decoy_z": np.nan,
              "specificity_na_reason": ""}
    if not np.isfinite(value):
        result["specificity_na_reason"] = "cognate_alignment_unavailable"
        return result
    if not len(decoys):
        result["specificity_na_reason"] = "no_usable_different_seed_decoys"
        return result
    # Near-equal correlations count as ties, avoiding float-roundoff rank changes.
    ties = np.isclose(decoys, value, rtol=1e-12, atol=1e-15)
    better = ((decoys > value) & ~ties).sum() + 0.5 * ties.sum()
    result.update(decoy_mean=float(decoys.mean()), cognate_rank=float(1 + better),
                  fraction_decoys_better=float(better / len(decoys)))
    if len(decoys) < 2:
        result["specificity_na_reason"] = "fewer_than_two_usable_decoys"
    else:
        sd = float(decoys.std(ddof=1))
        result["decoy_sd"] = sd
        if np.isclose(sd, 0.0, rtol=0, atol=1e-15):
            result["specificity_na_reason"] = "no_decoy_variation"
        else:
            result["decoy_z"] = float((value - decoys.mean()) / sd)
    return result


def run_analysis(config_path: pathlib.Path) -> pathlib.Path:
    config_path = config_path.expanduser().resolve()
    with config_path.open(encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    root = config_path.parent
    resolve = lambda p: resolve_predictor_output_path(p, root)
    input_dir = resolve(config["input_dir"])
    if not (input_dir / "comparisons.tsv").is_file():
        raise FileNotFoundError(f"no benchmark tables in {input_dir}; run `uv run funmirbench-download-data`")
    comps = pd.read_csv(input_dir / "comparisons.tsv", sep="\t", keep_default_na=False)
    if config.get("experiments"):
        comps = filter_df(comps, config["experiments"])
    skipped = comps[(comps.perturbation != "OE") | (comps.action == "set_aside")].copy()
    comps = comps[(comps.perturbation == "OE") & (comps.action != "set_aside")]
    if comps.empty or comps.experiment_id.duplicated().any():
        raise ValueError("select at least one OE comparison, with unique experiment IDs")
    mirnas = pd.read_csv(resolve(config["mirnas_tsv"]), sep="\t").set_index("mirna_name")
    if mirnas.index.has_duplicates or mirnas.index.hasnans:
        raise ValueError("miRNA registry needs unique, nonmissing names")
    sequences = mirnas.mirna_sequence.str.upper().str.replace("T", "U", regex=False)
    if not sequences.str.fullmatch("[ACGU]{8,}", na=False).all():
        raise ValueError("every registry miRNA needs an RNA/DNA sequence of at least eight bases")
    if not set(comps.mirna_name).issubset(mirnas.index):
        raise ValueError("every comparison miRNA needs a sequence in mirnas_tsv")
    seeds = sequences.str[1:8]
    tools = pd.read_csv(resolve(config["predictions_tsv"]), sep="\t")
    if config.get("predictors"):
        tools = filter_df(tools, config["predictors"])
    if tools.empty or tools.tool_id.duplicated().any():
        raise ValueError("select at least one predictor, with unique tool IDs")
    if not tools.score_direction.isin(["higher_is_stronger", "lower_is_stronger"]).all():
        raise ValueError("invalid predictor score_direction")
    validate_predictor_output_files({r["tool_id"]: r for r in tools.to_dict("records")}, root)
    final = pd.read_csv(input_dir / "experiments.tsv", sep="\t").set_index("experiment_id")
    if final.index.has_duplicates:
        raise ValueError("experiments.tsv needs unique experiment IDs")
    missing = set(comps.experiment_id) - set(final.index)
    if missing:
        raise ValueError(f"comparisons without a benchmark table: {sorted(missing)}")
    tables = {}
    for exp in comps.experiment_id:
        tables[exp] = pd.read_csv(input_dir / final.loc[exp, "table"], sep="\t",
                                  usecols=RESPONSE_COLUMNS).set_index("gene_id")
    universe = next(iter(tables.values())).index
    if any(not t.index.equals(universe) for t in tables.values()):
        raise ValueError("benchmark tables must share the same gene list in the same order")
    predictions = {}
    for tool in tools.itertuples(index=False):
        logger.info("loading %s", tool.tool_id)
        predictions[tool.tool_id] = load_predictions(
            resolve(tool.predictor_output_path), tool.score_direction, set(mirnas.index), set(universe))

    # A separate subfolder prevents collisions with the normal benchmark outputs.
    out = resolve(config.get("out_dir", "results")) / build_run_dir_name(
        experiments=None, tool_ids=None, eval_cfg=None) / "response_vectors"
    out.mkdir(parents=True, exist_ok=False)
    shutil.copy2(config_path, out / "analysis_config.yaml")
    comps.to_csv(out / "comparisons.tsv", sep="\t", index=False)
    skipped.to_csv(out / "skipped_comparisons.tsv", sep="\t", index=False)
    mirnas.to_csv(out / "mirnas.tsv", sep="\t")
    tools.to_csv(out / "predictors.tsv", sep="\t", index=False)
    candidate_rows, summary_rows, files = [], [], []
    for exp_number, comparison in enumerate(comps.itertuples(index=False)):
        table = tables[comparison.experiment_id]
        bins = expression_tenths(table.control_mean_normalized_count.fillna(0.0))
        for tool_number, (tool, by_mirna) in enumerate(predictions.items()):
            rows = []
            for candidate in mirnas.index:
                vectors = build_vectors(table, by_mirna.get(candidate, pd.Series(dtype=float)))
                rows.append({"experiment_id": comparison.experiment_id, "predictor": tool,
                             "candidate_mirna": candidate, "candidate_seed": seeds[candidate],
                             "is_cognate": candidate == comparison.mirna_name, **alignment(vectors, bins)})
                if candidate == comparison.mirna_name:
                    # Numeric filenames remain safe for user-supplied experiment/tool IDs.
                    filename = f"vectors/{exp_number:04d}_{tool_number:03d}.tsv.gz"
                    (out / "vectors").mkdir(exist_ok=True)
                    vectors.to_csv(out / filename, sep="\t", na_rep="NA")
                    files.append({"experiment_id": comparison.experiment_id, "predictor": tool, "file": filename})
            candidates = pd.DataFrame(rows)
            candidate_rows.extend(rows)
            own = candidates[candidates.is_cognate].iloc[0]
            for metric in ("spearman", "spearman_within_tenths"):
                summary_rows.append({"experiment_id": comparison.experiment_id, "predictor": tool,
                                     "mirna_name": comparison.mirna_name, "action": comparison.action,
                                     "dependency_group": comparison.dependency_group,
                                     "metric": metric, "alignment": own[metric],
                                     "alignment_na_reason": own.alignment_na_reason if metric == "spearman"
                                     else ("" if np.isfinite(own[metric]) else "no_usable_expression_bins"),
                                     **specificity(candidates, comparison.mirna_name, metric)})
        logger.info("analysed %s", comparison.experiment_id)
    pd.DataFrame(candidate_rows).to_csv(out / "candidate_alignment.tsv", sep="\t", index=False, na_rep="NA")
    pd.DataFrame(summary_rows).to_csv(out / "specificity.tsv", sep="\t", index=False, na_rep="NA")
    pd.DataFrame(files).to_csv(out / "vector_files.tsv", sep="\t", index=False)
    (out / "README.txt").write_text(
        "Gene-level response vector analysis\n\n"
        "Observed: oriented shrunken log2FC from the OE benchmark tables (positive = repression).\n"
        "Unusable responses remain NA; unscored predictions form one bottom tie.\n"
        "Unit centered average-rank vectors use the same response mask; their dot product equals Spearman.\n"
        "Scores are not calibrated effect sizes; vector lengths cannot be compared as predicted repression.\n"
        "Within-expression-tenth alignment follows the benchmark (minimum 20 genes per bin).\n"
        "specificity.tsv compares the cognate miRNA against finite different-seed registry candidates.\n"
        "cognate_rank = 1 + better decoys + half tied decoys; 1 is best.\n"
        "fraction_decoys_better = (better + half tied) / n_decoys; 0 is best.\n"
        "decoy_z uses sample SD and is NA for fewer than two decoys or no decoy variation.\n"
        "These are descriptive panel comparisons, not p-values or all-human seed-family ranks.\n"
        "Candidates with unavailable alignment are listed but omitted from the decoy denominator.\n"
        "Decoy families can have multiple members; counts refer to miRNAs, not independent families.\n"
        "Coverage and registry composition affect comparisons. No experiment is excluded by alignment.\n"
        "Replicate stability, controls-only contrasts and a seed-count baseline require additional inputs\n"
        "and are not computed by this command.\n", encoding="utf-8")
    logger.info("wrote %s", out)
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description="Inspect OE response vectors and miRNA specificity.")
    parser.add_argument("--config", type=pathlib.Path, required=True)
    parser.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"])
    args = parser.parse_args()
    setup_logging(parse_log_level(args.log_level))
    run_analysis(args.config)


if __name__ == "__main__":
    main()
