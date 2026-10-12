import tempfile
import unittest
from pathlib import Path

import pandas as pd

from funmirbench.sync_metadata import sync_metadata


class SyncMetadataTests(unittest.TestCase):
    def test_old_candidate_schema_fallbacks_are_per_row_and_nan_safe(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            registry = root / "registry.tsv"
            columns = [
                "id",
                "geo_accession",
                "mirna_name",
                "experiment_type",
                "pubmed_id",
                "control_samples",
                "condition_samples",
                "de_table_path",
            ]
            pd.DataFrame([{
                "id": "dataset",
                "geo_accession": "GSEOLD",
                "mirna_name": "hsa-miR-1",
                "experiment_type": "Overexpression",
                "pubmed_id": "11111",
                "control_samples": "GSM1",
                "condition_samples": "GSM2",
                "de_table_path": "old.tsv",
            }], columns=columns).to_csv(registry, sep="\t", index=False)

            old_run = root / "old"
            new_run = root / "new"
            old_run.mkdir()
            new_run.mkdir()
            pd.DataFrame([{
                "id": "dataset",
                "gse_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123",
                "article_pubmed_id": "https://pubmed.ncbi.nlm.nih.gov/31101765/",
                "mirna_name": "hsa-miR-1",
                "experiment_type": "Overexpression",
                "de_table_path": "new.tsv",
            }]).to_csv(old_run / "candidate_metadata.tsv", sep="\t", index=False)
            pd.DataFrame([{
                "id": "other",
                "geo_accession": "GSE456",
                "pubmed_id": "22222",
                "mirna_name": "hsa-miR-2",
                "experiment_type": "Inhibition",
                "de_table_path": "other.tsv",
            }]).to_csv(new_run / "candidate_metadata.tsv", sep="\t", index=False)

            sync_metadata(inputs=[root], registry_path=registry)

            synced = pd.read_csv(registry, sep="\t", dtype=str).fillna("")
            row = synced.set_index("id").loc["dataset"]
            self.assertEqual(row["geo_accession"], "GSE123")
            self.assertEqual(row["pubmed_id"], "31101765")
            self.assertEqual(row["control_samples"], "GSM1")
            self.assertEqual(row["condition_samples"], "GSM2")


if __name__ == "__main__":
    unittest.main()
