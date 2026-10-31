import importlib
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
STANDARDIZED_PREDICTORS = ROOT / "pipelines" / "standardized_predictors"
if str(STANDARDIZED_PREDICTORS) not in sys.path:
    sys.path.insert(0, str(STANDARDIZED_PREDICTORS))
targetscan_utils = importlib.import_module("pipelines.standardized_predictors.targetscan.utils")


FAMILY_INFO_COLUMNS = [
    "miR family",
    "Seed+m8",
    "Species ID",
    "MiRBase ID",
    "Mature sequence",
    "Family Conservation?",
    "MiRBase Accession",
]
FAMILY_INFO_ROWS = [
    ["let-7-5p/98-5p", "GAGGUAG", "9606", "hsa-let-7a-5p", "UGAGGUAGUAGGUUGUAUAGUU", "2", "MIMAT0000062"],
    ["let-7-5p/98-5p", "GAGGUAG", "9606", "hsa-let-7b-5p", "UGAGGUAGUAGGUUGUGUGGUU", "2", "MIMAT0000063"],
    ["let-7-5p/98-5p", "GAGGUAG", "9606", "hsa-miR-98-5p", "UGAGGUAGUAAGUUGUAUUGUU", "2", "MIMAT0000096"],
    ["miR-1-3p/206", "GGAAUGU", "9606", "hsa-miR-1-3p", "UGGAAUGUAAAGAAGUAUGUAU", "2", "MIMAT0000416"],
    ["let-7/98", "GAGGUAG", "10090", "mmu-let-7a-5p", "UGAGGUAGUAGGUUGUAUAGUU", "2", "MIMAT0000521"],
]
MIRBASE_ACC2NAME = {
    "MIMAT0000062": "hsa-let-7a-5p",
    "MIMAT0000063": "hsa-let-7b-5p",
    "MIMAT0000096": "hsa-miR-98-5p",
    "MIMAT0000416": "hsa-miR-1-3p",
    "MIMAT0000521": "mmu-let-7a-5p",
}
SUMMARY_COLUMNS = [
    "Transcript ID",
    "Gene Symbol",
    "miRNA family",
    "Species ID",
    "Representative miRNA",
    "Cumulative weighted context++ score",
]
SUMMARY_ROWS = [
    ["ENST00000000001.1", "GENE1", "GAGGUAG", "9606", "hsa-let-7a-5p", "-0.5"],
    ["ENST00000000002.1", "GENE2", "GAGGUAG", "9606", "hsa-let-7b-5p", "-0.3"],
    ["ENST00000000001.1", "GENE1", "GGAAUGU", "9606", "hsa-miR-1-3p", "-0.2"],
    ["ENST00000000001.1", "GENE1", "GAGGUAG", "10090", "mmu-let-7a-5p", "-0.9"],
]
TX_INDEX = {
    "rep_txs_by_gene_id": {
        "ENSG00000000001": [("ENST00000000001.1", "", "GENE1")],
        "ENSG00000000002": [("ENST00000000002.1", "", "GENE2")],
    },
    "targetscan_gene_id_by_tx": {
        "ENST00000000001.1": "ENSG00000000001",
        "ENST00000000002.1": "ENSG00000000002",
    },
}
ENSEMBL_TABLES = {
    "tx_to_gene": {
        "ENST00000000001": "ENSG00000000001",
        "ENST00000000002": "ENSG00000000002",
    },
    "gene_to_name": {
        "ENSG00000000001": "GENE1",
        "ENSG00000000002": "GENE2",
    },
}


def _write_tsv(path, columns, rows):
    path.write_text(
        "\t".join(columns) + "\n" + "".join("\t".join(row) + "\n" for row in rows),
        encoding="utf-8",
    )


class TargetScanFamilyExpansionTests(unittest.TestCase):
    def setUp(self):
        self._temporary_directory = tempfile.TemporaryDirectory()
        self.root = Path(self._temporary_directory.name)
        self.family_info = self.root / "miR_Family_Info.txt"
        self.summary_counts = self.root / "Summary_Counts.all_predictions.txt"
        _write_tsv(self.family_info, FAMILY_INFO_COLUMNS, FAMILY_INFO_ROWS)
        _write_tsv(self.summary_counts, SUMMARY_COLUMNS, SUMMARY_ROWS)
        self.annotations = targetscan_utils.step_build_human_mirna_annotations(
            self.family_info,
            mirbase_acc2name=MIRBASE_ACC2NAME,
        )

    def tearDown(self):
        self._temporary_directory.cleanup()

    def _write(self, output_name, **kwargs):
        output_path = self.root / output_name
        targetscan_utils.step_write_standardized_predictions(
            self.summary_counts,
            tx_index=TX_INDEX,
            ensembl_tables=ENSEMBL_TABLES,
            mirna_annotations=self.annotations,
            output_path=output_path,
            **kwargs,
        )
        return pd.read_csv(output_path, sep="\t", dtype={"Score": str})

    def test_family_members_are_read_per_seed_for_human_only(self):
        members = targetscan_utils.step_build_human_mirna_family_members(self.family_info)

        self.assertEqual(
            members,
            {
                "GAGGUAG": ["hsa-let-7a-5p", "hsa-let-7b-5p", "hsa-miR-98-5p"],
                "GGAAUGU": ["hsa-miR-1-3p"],
            },
        )

    def test_default_output_keeps_representative_miRNA_only(self):
        out = self._write("legacy.tsv")

        self.assertEqual(list(out.columns), ["Ensembl_ID", "Gene_Name", "miRNA_ID", "miRNA_Name", "Score"])
        self.assertEqual(
            sorted(zip(out["Ensembl_ID"], out["miRNA_Name"])),
            [
                ("ENSG00000000001", "hsa-let-7a-5p"),
                ("ENSG00000000001", "hsa-miR-1-3p"),
                ("ENSG00000000002", "hsa-let-7b-5p"),
            ],
        )

    def test_opt_in_expansion_maps_family_scores_to_every_member(self):
        members = targetscan_utils.step_build_human_mirna_family_members(self.family_info)
        out = self._write("expanded.tsv", family_members=members)

        self.assertEqual(
            list(out.columns),
            [
                "Ensembl_ID",
                "Gene_Name",
                "miRNA_ID",
                "miRNA_Name",
                "Score",
                "TargetScan_miRNA_Family",
                "TargetScan_Representative_miRNA",
            ],
        )
        self.assertEqual(len(out), 7)
        gene1_let7 = out[
            (out["Ensembl_ID"] == "ENSG00000000001")
            & (out["TargetScan_miRNA_Family"] == "GAGGUAG")
        ]
        self.assertEqual(
            sorted(gene1_let7["miRNA_Name"]),
            ["hsa-let-7a-5p", "hsa-let-7b-5p", "hsa-miR-98-5p"],
        )
        self.assertEqual(set(gene1_let7["Score"]), {"-0.5"})
        self.assertEqual(set(gene1_let7["TargetScan_Representative_miRNA"]), {"hsa-let-7a-5p"})
        gene2_mir98 = out[
            (out["Ensembl_ID"] == "ENSG00000000002") & (out["miRNA_Name"] == "hsa-miR-98-5p")
        ]
        self.assertEqual(gene2_mir98["miRNA_ID"].tolist(), ["MIMAT0000096"])
        self.assertEqual(gene2_mir98["Score"].tolist(), ["-0.3"])
        self.assertEqual(gene2_mir98["TargetScan_Representative_miRNA"].tolist(), ["hsa-let-7b-5p"])


if __name__ == "__main__":
    unittest.main()
