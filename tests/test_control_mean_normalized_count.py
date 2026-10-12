import tempfile
import unittest
from pathlib import Path

import pandas as pd

from funmirbench.validate_experiments import validate_experiments


class ControlMeanNormalizedCountTests(unittest.TestCase):
    def test_validator_requires_control_mean_normalized_count(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            de_table = root / "de.tsv"
            metadata = root / "metadata.tsv"

            pd.DataFrame({
                "gene_id": ["ENSG000001"],
                "logFC": [-2.0],
                "FDR": [0.01],
            }).to_csv(de_table, sep="\t", index=False)
            pd.DataFrame({
                "id": ["dataset"],
                "mirna_name": ["hsa-miR-1"],
                "organism": ["Homo sapiens"],
                "experiment_type": ["Overexpression"],
                "de_table_path": [de_table.name],
            }).to_csv(metadata, sep="\t", index=False)

            summary = validate_experiments(metadata, root=root)

            self.assertEqual([issue.check for issue in summary.issues], ["required_de_columns"])
            self.assertIn("control_mean_normalized_count", summary.issues[0].message)


if __name__ == "__main__":
    unittest.main()
