import unittest

from funmirbench.run_report_pdf import _paginate_table_rows


class RunReportPdfTests(unittest.TestCase):
    def test_inventory_rows_paginate_without_losing_order(self):
        rows = [
            [
                f"Experiment {index} with a descriptive label that wraps",
                f"GSE{index:06d}_OE_miR_example_cell_line",
                "No selected predictor had evaluable scored rows after filtering.",
            ]
            for index in range(30)
        ]

        pages = _paginate_table_rows(
            rows,
            col_widths=[0.36, 0.27, 0.37],
        )

        self.assertGreater(len(pages), 1)
        self.assertEqual([row for page in pages for row in page], rows)


if __name__ == "__main__":
    unittest.main()
