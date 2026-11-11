"""Render representative layouts and check geometry rather than image pixels."""
import pathlib
import unittest
from unittest.mock import patch

import matplotlib.pyplot as plt
import pandas as pd

from funmirbench import plots, report_pdf


class LayoutTests(unittest.TestCase):
    def test_ten_predictor_legends_do_not_overlap_or_cover_panels(self):
        ids = [f"tool_{i}" for i in range(10)]
        tools = pd.DataFrame({"tool_id": ids, "official_name": [f"Synthetic predictor {i}" for i in range(10)]})
        per = pd.DataFrame([{"experiment_id": f"E{e}", "predictor": t,
                             "auroc_d05": 0.5 + i / 30, "spearman_shrunk": i / 30,
                             "spearman_above_decoys_shrunk": i / 60}
                            for e in range(3) for i, t in enumerate(ids)])
        coverage = pd.DataFrame({"predictor": ids, "genes_covered": 100, "genes_measured": 100,
                                 "experiments_retained": 1000, "experiments_total": 1000,
                                 "responsive_scored": 8, "responsive_total": 10,
                                 "no_response_scored": 70, "no_response_total": 90})
        recovery = pd.DataFrame([{"experiment_id": "E1", "predictor": t, "budget": b, "recall": b / 100}
                                 for t in ids for b in [10, 50, 100]])
        precision = pd.DataFrame({"experiment_id": "E1", "predictor": ids, "precision": 0.2,
                                  "random_precision": 0.1, "precision_over_random": 2.,
                                  "false_positives_per_true_positive": 4.})

        def check(fig, out, stem):
            try:
                fig.canvas.draw()
                renderer = fig.canvas.get_renderer()
                legend = fig.legends[0]
                boxes = [t.get_window_extent(renderer) for t in legend.get_texts()]
                for i, box in enumerate(boxes):
                    for other in boxes[i + 1:]:
                        self.assertFalse(box.overlaps(other), stem)
                bounds = legend.get_window_extent(renderer)
                for ax in fig.axes:
                    self.assertFalse(bounds.overlaps(ax.get_window_extent(renderer)), stem)
                    for label in [ax.xaxis.label, *ax.get_xticklabels()]:
                        if label.get_text():
                            self.assertFalse(bounds.overlaps(label.get_window_extent(renderer)), stem)
                    if stem == "coverage":
                        percentages = [text.get_window_extent(renderer) for text in ax.texts]
                        for i, box in enumerate(percentages):
                            for other in percentages[i + 1:]:
                                self.assertFalse(box.overlaps(other), "coverage percentage labels")
                self.assertEqual(len(boxes), 10)
                self.assertLessEqual(fig.get_figwidth(), 16)
            finally:
                plt.close(fig)

        with patch.object(plots, "save", check):
            plots.plot_coverage(coverage, tools, pathlib.Path("unused"))
            plots.plot_all_genes(per, tools, pathlib.Path("unused"))
            plots.plot_recovery(recovery, precision, tools, pathlib.Path("unused"))

    def test_predictor_tables_paginate_with_headers_inside_page(self):
        rows = [[f"Synthetic predictor {i}", "0.500"] for i in range(40)]
        pages = []

        def check(pdf, fig):
            try:
                fig.canvas.draw()
                renderer = fig.canvas.get_renderer()
                for table in fig.axes[0].tables:
                    box = table.get_window_extent(renderer).transformed(fig.transFigure.inverted())
                    self.assertGreaterEqual(box.y0, report_pdf.BOTTOM - 1e-9)
                    self.assertLessEqual(box.y1, report_pdf.TOP)
                    pages.append(table)
            finally:
                plt.close(fig)

        with patch.object(report_pdf, "_save_page", check):
            report_pdf.draw_tables(None, [("Performance", ["Predictor", "Mean AUROC"], [0.7, 0.3], rows)])
        self.assertGreater(len(pages), 1)

    def test_small_ranked_first_counts_are_separated(self):
        from matplotlib.text import Text
        ids = [f"tool_{i}" for i in range(10)]
        per = pd.DataFrame([{"experiment_id": f"E{e}", "predictor": tool,
                             "auroc_d05": float(i == (e if e < 9 else 9)),
                             "spearman_shrunk": float(i == (e if e < 9 else 9))}
                            for e in range(72) for i, tool in enumerate(ids)])
        fig, ax = plt.subplots(figsize=(7, 4))
        try:
            plots.draw_won(ax, per, ids, plots.tool_colors(ids))
            fig.canvas.draw()
            renderer = fig.canvas.get_renderer()
            labels = [text for text in ax.texts if text.get_text().isdigit()]
            self.assertEqual(len(labels), 20)
            boxes = [Text.get_window_extent(text, renderer) for text in labels]
            for i, box in enumerate(boxes):
                for other in boxes[i + 1:]:
                    self.assertFalse(box.overlaps(other))
        finally:
            plt.close(fig)


if __name__ == "__main__":
    unittest.main()
