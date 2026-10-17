import tempfile
import unittest
from pathlib import Path

import pandas as pd

from funmirbench.evaluate_common import _annotate_ground_truth, _normalize_perturbation
from funmirbench.validate_experiments import (
    _normalize_perturbation as _normalize_registry_perturbation,
    validate_experiments,
)


ALIASES = [
    ("OE", "Overexpression"),
    ("oe", "Overexpression"),
    ("KO", "Knockout"),
    ("ko", "Knockout"),
    ("KD", "Knockdown"),
    ("Kd", "Knockdown"),
    (" overexpression ", "Overexpression"),
]


class PerturbationAliasTests(unittest.TestCase):
    def test_evaluator_canonicalizes_abbreviations(self):
        for value, expected in ALIASES:
            with self.subTest(value=value):
                self.assertEqual(_normalize_perturbation(value), expected)

    def test_validator_canonicalizes_abbreviations(self):
        for value, expected in ALIASES:
            with self.subTest(value=value):
                self.assertEqual(_normalize_registry_perturbation(value), expected)

    def test_abbreviations_drive_expected_effect_sign(self):
        frame = pd.DataFrame({"logFC": [-2.0, 2.0], "FDR": [0.01, 0.01]})

        overexpression = _annotate_ground_truth(frame, perturbation="OE")
        knockdown = _annotate_ground_truth(frame, perturbation="kd")

        self.assertEqual(overexpression["expected_effect"].tolist(), [2.0, -2.0])
        self.assertEqual(knockdown["expected_effect"].tolist(), [-2.0, 2.0])

    def test_validator_accepts_abbreviated_registry_experiment_type(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            pd.DataFrame({
                "gene_id": ["ENSG000001", "ENSG000002"],
                "logFC": [2.0, 0.0],
                "control_mean_normalized_count": [10.0, 20.0],
                "FDR": [0.01, 0.5],
            }).to_csv(root / "de.tsv", sep="\t", index=False)
            pd.DataFrame({
                "id": ["dataset"],
                "mirna_name": ["hsa-miR-1"],
                "organism": ["Homo sapiens"],
                "experiment_type": ["KO"],
                "de_table_path": ["de.tsv"],
            }).to_csv(root / "metadata.tsv", sep="\t", index=False)

            summary = validate_experiments(root / "metadata.tsv", root=root)

        self.assertEqual(summary.issues, [])
        self.assertEqual(summary.notices, [])
        self.assertEqual(summary.excluded_dataset_ids, set())


if __name__ == "__main__":
    unittest.main()
