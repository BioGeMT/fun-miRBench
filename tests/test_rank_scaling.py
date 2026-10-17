import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from funmirbench.cross_dataset import _rank_values_for_spec
from funmirbench.evaluate import evaluate_joined_dataframe
from funmirbench.evaluate_common import _rank_scale_scores
from funmirbench.predictor_combinations import _prepare_combo_frame


def _single_valued_joined():
    return pd.DataFrame({
        "gene_id": ["ENSG000001", "ENSG000002", "ENSG000003", "ENSG000004"],
        "logFC": [-2.0, -2.0, 0.0, 0.0],
        "FDR": [0.01, 0.01, 0.5, 0.5],
        "perturbation": "Overexpression",
        "score_binary": [1.0, np.nan, 1.0, np.nan],
    })


class RankScalingTests(unittest.TestCase):
    def test_single_valued_scores_keep_missing_rows_missing(self):
        ranks = _rank_scale_scores(pd.Series([2.0, np.nan, 2.0, np.nan]))

        self.assertEqual(ranks.tolist()[0], 1.0)
        self.assertEqual(ranks.tolist()[2], 1.0)
        self.assertTrue(ranks.iloc[[1, 3]].isna().all())

    def test_multi_valued_scores_keep_missing_rows_missing(self):
        ranks = _rank_scale_scores(pd.Series([1.0, np.nan, 3.0]))

        self.assertEqual(ranks.iloc[0], 0.0)
        self.assertTrue(pd.isna(ranks.iloc[1]))
        self.assertEqual(ranks.iloc[2], 1.0)

    def test_evaluator_local_rank_keeps_missing_scores_missing(self):
        joined = _single_valued_joined()
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            evaluate_joined_dataframe(
                joined,
                plots_dir=root / "plots",
                reports_dir=root / "reports",
                fdr_threshold=0.05,
                effect_threshold=1.0,
                predictor_top_fraction=0.1,
                dataset_id="dataset",
                perturbation="Overexpression",
            )

        self.assertEqual(joined["local_rank_binary"].notna().tolist(), [True, False, True, False])

    def test_cross_dataset_score_ranks_keep_missing_scores_missing(self):
        values = _rank_values_for_spec(_single_valued_joined(), "score_binary", "score")

        self.assertEqual(values.notna().tolist(), [True, False, True, False])

    def test_combination_coverage_is_taken_from_raw_scores(self):
        prepared = _prepare_combo_frame(
            _single_valued_joined(),
            ("binary",),
            fdr_threshold=0.05,
            effect_threshold=1.0,
        )

        self.assertIsNotNone(prepared)
        work, coverage_info = prepared
        self.assertEqual(coverage_info["rows_scored"], 2)
        self.assertEqual(coverage_info["coverage"], 0.5)
        self.assertEqual(sorted(work["gene_id"]), ["ENSG000001", "ENSG000003"])


if __name__ == "__main__":
    unittest.main()
