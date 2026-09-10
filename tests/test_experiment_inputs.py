import gzip
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from funmirbench.experiments_pipeline import (
    combine_count_matrices,
    group_sample_ids_by_layout,
    infer_library_layout,
    materialize_sample_reads,
    normalize_sample_entry,
    run_reads_mode,
)


class ExperimentInputTests(unittest.TestCase):
    def test_single_path_config_remains_supported(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            reads = root / "sample.fastq.gz"
            reads.write_bytes(b"reads")

            sample = normalize_sample_entry(
                {"sample_id": "GSM1", "reads_1": str(reads)},
                group_name="control",
                root=root,
                repo=root,
            )
            materialized = materialize_sample_reads(sample, run_dir=root / "run")

            self.assertEqual(materialized["reads_1"], str(reads))
            self.assertEqual(materialized["reads_2"], "")

    def test_multiple_runs_are_materialized_as_one_biological_sample(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            run_1 = root / "run_1.fastq.gz"
            run_2 = root / "run_2.fastq.gz"
            with gzip.open(run_1, "wb") as handle:
                handle.write(b"@run1\nAC\n+\nII\n")
            with gzip.open(run_2, "wb") as handle:
                handle.write(b"@run2\nGT\n+\nII\n")

            sample = normalize_sample_entry(
                {"sample_id": "GSM1", "reads_1": [str(run_1), str(run_2)]},
                group_name="control",
                root=root,
                repo=root,
            )
            materialized = materialize_sample_reads(sample, run_dir=root / "run")

            self.assertEqual(materialized["sample_id"], "GSM1")
            self.assertEqual(materialized["source_reads_1"], [str(run_1), str(run_2)])
            with gzip.open(materialized["reads_1"], "rb") as handle:
                self.assertEqual(
                    handle.read(),
                    b"@run1\nAC\n+\nII\n@run2\nGT\n+\nII\n",
                )

    def test_paired_run_lists_must_have_equal_lengths(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            paths = [root / name for name in ("r1a.fastq", "r1b.fastq", "r2a.fastq")]
            for path in paths:
                path.write_text("reads", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "same number"):
                normalize_sample_entry(
                    {
                        "sample_id": "GSM1",
                        "reads_1": [str(paths[0]), str(paths[1])],
                        "reads_2": [str(paths[2])],
                    },
                    group_name="control",
                    root=root,
                    repo=root,
                )

    def test_mixed_library_layout_is_supported(self):
        samples = [
            {"sample_id": "single", "reads_2": ""},
            {"sample_id": "paired", "reads_2": "mate.fastq.gz"},
        ]

        self.assertEqual(infer_library_layout(samples), "mixed")
        self.assertEqual(
            group_sample_ids_by_layout(samples, ["paired", "single"]),
            {"single": ["single"], "paired": ["paired"]},
        )

    def test_separate_featurecounts_matrices_are_combined_in_sample_order(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            single = root / "single.tsv"
            paired = root / "paired.tsv"
            output = root / "combined.tsv"
            pd.DataFrame({"gene_id": ["ENSG1", "ENSG2"], "single": [1, 2]}).to_csv(
                single, sep="\t", index=False
            )
            pd.DataFrame({"gene_id": ["ENSG1", "ENSG2"], "paired": [3, 4]}).to_csv(
                paired, sep="\t", index=False
            )

            combine_count_matrices([single, paired], ["paired", "single"], output)

            combined = pd.read_csv(output, sep="\t")
            self.assertEqual(list(combined.columns), ["gene_id", "paired", "single"])
            self.assertEqual(combined.to_dict(orient="list"), {
                "gene_id": ["ENSG1", "ENSG2"],
                "paired": [3, 4],
                "single": [1, 2],
            })

    def test_reads_mode_splits_mixed_layouts_before_featurecounts(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            single_reads = root / "single.fastq.gz"
            paired_reads_1 = root / "paired_1.fastq.gz"
            paired_reads_2 = root / "paired_2.fastq.gz"
            for path in (single_reads, paired_reads_1, paired_reads_2):
                path.write_bytes(b"reads")
            (root / "run").mkdir()

            samples = (
                [{
                    "sample_id": "single",
                    "group": "control",
                    "reads_1": [str(single_reads)],
                    "reads_2": [],
                    "count_matrix_column": "single",
                    "accession": "single",
                    "title": "single",
                    "group_label": "control",
                }],
                [{
                    "sample_id": "paired",
                    "group": "treated",
                    "reads_1": [str(paired_reads_1)],
                    "reads_2": [str(paired_reads_2)],
                    "count_matrix_column": "paired",
                    "accession": "paired",
                    "title": "paired",
                    "group_label": "treated",
                }],
            )
            featurecounts_calls = []

            def fake_fastqc(**kwargs):
                return {}, [], root / "stdout", root / "stderr"

            def fake_fastp(**kwargs):
                return kwargs["sample"], [], root / "stdout", root / "stderr"

            def fake_star(**kwargs):
                sample_id = kwargs["sample"]["sample_id"]
                return root / f"{sample_id}.bam", [], root / "stdout", root / "stderr"

            def fake_featurecounts(**kwargs):
                featurecounts_calls.append(
                    (kwargs["sample_order"], kwargs["paired_end"])
                )
                return root / "raw.tsv", [], root / "stdout", root / "stderr"

            def fake_build_featurecounts_matrix(**kwargs):
                sample_order = kwargs["sample_order"]
                pd.DataFrame({
                    "gene_id": ["ENSG1"],
                    **{sample_id: [1] for sample_id in sample_order},
                }).to_csv(kwargs["out_path"], sep="\t", index=False)
                return kwargs["out_path"]

            captured = {}

            def fake_run_de_from_counts(**kwargs):
                captured.update(kwargs)
                return {"run_dir": str(root), "de_table_path": "output.tsv"}

            with (
                patch("funmirbench.experiments_pipeline.load_reads_samples", return_value=samples),
                patch(
                    "funmirbench.experiments_pipeline.prepare_reads_reference_assets",
                    return_value={
                        "star_index": root / "star",
                        "gtf_path": root / "genes.gtf",
                        "genome_fasta": str(root / "genome.fa"),
                        "generated_star_index": False,
                        "reused_star_index": True,
                    },
                ),
                patch("funmirbench.experiments_pipeline.run_fastqc", side_effect=fake_fastqc),
                patch("funmirbench.experiments_pipeline.run_fastp", side_effect=fake_fastp),
                patch("funmirbench.experiments_pipeline.run_star_alignment", side_effect=fake_star),
                patch("funmirbench.experiments_pipeline.run_featurecounts", side_effect=fake_featurecounts),
                patch(
                    "funmirbench.experiments_pipeline.build_featurecounts_matrix",
                    side_effect=fake_build_featurecounts_matrix,
                ),
                patch(
                    "funmirbench.experiments_pipeline.run_de_from_counts",
                    side_effect=fake_run_de_from_counts,
                ),
            ):
                run_reads_mode(
                    {"dataset_id": "dataset", "source": {}},
                    config_path=root / "config.yaml",
                    repo=root,
                    run_dir=root / "run",
                    force=False,
                )

            self.assertEqual(
                featurecounts_calls,
                [(["single"], False), (["paired"], True)],
            )
            self.assertEqual(
                list(captured["counts_df"].columns),
                ["gene_id", "single", "paired"],
            )
            self.assertEqual(
                captured["extra_manifest"]["library_layout"], "mixed"
            )


if __name__ == "__main__":
    unittest.main()
