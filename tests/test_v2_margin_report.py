import unittest

import numpy as np
import pandas as pd

from funmirbench import v2_margin_report as mr


class MarginRowsTests(unittest.TestCase):
    def test_counts_ratio_and_site_shares(self):
        # Genes 0-3 moved like targets, 4-5 the other way, 6-7 barely; gene 8 has no usable estimate.
        o = np.array([2.0, 2.0, 2.0, 2.0, -2.0, -2.0, 0.0, 0.0, np.nan])
        sval = np.array([0.01] * 6 + [0.9, 0.9, np.nan])
        state = ["responsive"] * 4 + ["no_response_above_margin"] * 4 + ["unresolved"]
        table = pd.DataFrame({"eligible": [True] * 8 + [False], "logFC_shrunk": -o, "lfcSE_shrunk": 0.1,
                              "oriented_logFC_shrunk": o, "oriented_ci95_low_shrunk": o - 0.2,
                              "oriented_ci95_high_shrunk": o + 0.2, "state_d05": state, "svalue_d05": sval})
        site = np.array([True, True, True, False, False, False, False, False, False])
        (row,) = mr.margin_rows(table, {"d05": 0.5}, site)
        self.assertEqual((row["n_responsive"], row["n_no_response"], row["n_unresolved"]), (4, 4, 1))
        self.assertEqual(row["n_opposite"], 2)
        self.assertAlmostEqual(row["opposite_to_target_ratio"], 0.5)
        self.assertAlmostEqual(row["opposite_median_change"], -2.0)
        self.assertAlmostEqual(row["site_share_target"], 0.75)
        self.assertAlmostEqual(row["site_share_opposite"], 0.0)
        self.assertAlmostEqual(row["site_share_all"], 3 / 8)
        self.assertLess(row["expected_wrong_share_responsive"], 1e-6)


if __name__ == "__main__":
    unittest.main()
