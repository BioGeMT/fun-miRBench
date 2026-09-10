import importlib
import sys
import tempfile
import types
import unittest
from pathlib import Path

bio_module = types.ModuleType("Bio")
bio_module.Entrez = types.SimpleNamespace()
sys.modules.setdefault("Bio", bio_module)

geo_download = importlib.import_module("pipelines.geo.geo_download")
DEFAULT_GENOME_FASTA = geo_download.DEFAULT_GENOME_FASTA
DEFAULT_GTF = geo_download.DEFAULT_GTF
SRRInfo = geo_download.SRRInfo
build_geo_sample_entries = geo_download.build_geo_sample_entries
generate_yaml_config = geo_download.generate_yaml_config


class GeoDownloadTests(unittest.TestCase):
    def test_multiple_runs_remain_one_gsm_sample(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            output_directory = Path(temporary_directory)
            runs = [
                SRRInfo(srr="SRR2", gsm="GSM1", layout="PAIRED"),
                SRRInfo(srr="SRR1", gsm="GSM1", layout="PAIRED"),
            ]

            control, treated = build_geo_sample_entries(
                runs, ["GSM1"], [], output_directory
            )

            self.assertEqual(treated, [])
            self.assertEqual(control, [{
                "sample_id": "GSM1",
                "reads_1": [
                    str(output_directory / "SRR1_1.fastq.gz"),
                    str(output_directory / "SRR2_1.fastq.gz"),
                ],
                "reads_2": [
                    str(output_directory / "SRR1_2.fastq.gz"),
                    str(output_directory / "SRR2_2.fastq.gz"),
                ],
            }])

    def test_generated_config_preserves_multiple_run_paths(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            output_directory = Path(temporary_directory)
            entry = {
                "sample_id": "GSM1",
                "reads_1": ["run1_1.fastq.gz", "run2_1.fastq.gz"],
                "reads_2": ["run1_2.fastq.gz", "run2_2.fastq.gz"],
            }
            experiment = {
                "id": "dataset",
                "mirna": "hsa-miR-1-3p",
                "experiment_type": "Overexpression",
                "gse": "GSE1",
            }

            config_path = generate_yaml_config(
                experiment, [entry], [dict(entry, sample_id="GSM2")], output_directory
            )

            config_text = config_path.read_text(encoding="utf-8")
            self.assertIn("- run1_1.fastq.gz", config_text)
            self.assertIn("- run2_2.fastq.gz", config_text)

    def test_generated_configs_use_downloaded_ensembl_v115_references(self):
        self.assertIn("ensembl_v115", DEFAULT_GENOME_FASTA)
        self.assertIn("GRCh38.115.gtf.gz", DEFAULT_GTF)


if __name__ == "__main__":
    unittest.main()
