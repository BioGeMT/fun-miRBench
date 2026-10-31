"""REPORT.pdf of the benchmark."""
from __future__ import annotations

import pathlib
import textwrap

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

from funmirbench.plots import best_per_experiment, drawn_tools, experiments_won

REPORT_PAGE_SIZE = (8.27, 11.69)

PUBLICATION_BLUE = "#17324D"
PUBLICATION_RULE = "#D8DEE9"
PUBLICATION_TABLE_HEADER = "#E9F1FB"
PUBLICATION_TABLE_ALT = "#F9FBFD"
REPORT_TITLE_SIZE = 23
REPORT_SECTION_SIZE = 13.4
REPORT_BODY_SIZE = 11.2
REPORT_TABLE_SIZE = 10.2
INVENTORY_WRAP_SCALE = 75
INVENTORY_ROWS, ROW_HEIGHT = 24.0, 0.0345  # rows per page and their height, as a share of the page height
MARGIN_INCHES, GAP_INCHES = 0.45, 0.3  # around every page, and between the blocks on a page
LEFT = MARGIN_INCHES / REPORT_PAGE_SIZE[0]
WIDTH = 1.0 - 2 * LEFT
TOP = 1.0 - MARGIN_INCHES / REPORT_PAGE_SIZE[1]
BOTTOM = MARGIN_INCHES / REPORT_PAGE_SIZE[1]
BLOCK_GAP = GAP_INCHES / REPORT_PAGE_SIZE[1]
HEADER_RULE_DROP = 0.04  # from the top of the page title to the rule under it
TITLE_DROP, LINE_STEP, PARAGRAPH_STEP = 0.038, 0.030, 0.007  # block title to first line, lines, extra after each item
PLOT_BORDER_PIXELS = 45  # the white border of the plot PNGs (0.15 inches at 300 dpi), cropped so that the margins hold
OUTPUT_TABLES = ["per_comparison.tsv: every metric, for each experiment and tool.",
                 "summary.tsv: means and 95% intervals, on all experiments and on unflagged ones.",
                 "coverage.tsv: genes and pairs covered, and experiments retained, per tool.",
                 "recovery.tsv: recall per experiment and tool, at each prediction budget.",
                 "precision_at_recall.tsv: precision and false positives at 50% recall."]
COVERAGE = ("coverage", "Coverage",
          ["Gene coverage: share of measured genes the tool scores in at least one experiment.",
           "Experiments retained: experiments where it scores responsive and no-response genes.",
           "Pair coverage: share of responsive and of no-response miRNA-gene pairs it scores."])
PERFORMANCE = ("oe_all_genes", "Performance",
          ["AUROC: chance that a responsive gene scores above a no-response gene; 0.5 is chance.",
           "Spearman rho: agreement of the score's ranking with the shrunken log2FC; 0 is none.",
           "Decoy z-score: Spearman for the tested miRNA against different-seed miRNAs, in SDs.",
           "Experiments won: experiments where the tool has the highest AUROC or Spearman rho."])
RECOVERY = ("recovery", "Recovery",
          ["Recall: share of responsive genes among the top-ranked genes, mean over experiments.",
           "Best precision / random: best tool's precision at 50% recall over responsive share."])


def _format_summary_value(value, *, percent=False):
    if pd.isna(value):
        return "NA"
    number = float(value)
    if percent:
        return f"{number:.1%}"
    return f"{number:.3f}"

def _new_page():
    fig = plt.figure(figsize=REPORT_PAGE_SIZE)
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    ax.axis("off")
    fig.patch.set_facecolor("white")
    return fig, ax

def _save_page(pdf, fig):
    fig.patch.set_facecolor("white")
    pdf.savefig(fig, facecolor="white")
    plt.close(fig)

def _header(ax, title) -> float:
    """Draw the page title and its rule; returns where the next block starts."""
    ax.text(LEFT, TOP, title, fontsize=REPORT_TITLE_SIZE, fontweight="bold", color=PUBLICATION_BLUE, va="top", ha="left", family="DejaVu Sans")
    rule = TOP - HEADER_RULE_DROP
    ax.add_line(plt.Line2D([LEFT, LEFT + WIDTH], [rule, rule], color=PUBLICATION_RULE, linewidth=1.2))
    return rule - BLOCK_GAP

def _bullet_block(ax, title, lines, *, y, body_size=REPORT_BODY_SIZE, title_size=REPORT_SECTION_SIZE) -> float:
    """Draw a titled list from y down; returns where the next block starts."""
    ax.text(LEFT, y, title, fontsize=title_size, fontweight="bold", color="#2F5D8C", va="top", ha="left", family="DejaVu Sans")
    current_y = y - TITLE_DROP
    wrap_width = int(WIDTH * 96)
    for line in lines:
        chunks = textwrap.wrap(str(line), width=wrap_width) or [""]
        for index, chunk in enumerate(chunks):
            prefix = "- " if index == 0 else "  "
            ax.text(
                LEFT,
                current_y,
                prefix + chunk,
                fontsize=body_size,
                color="#22303C",
                va="top",
                ha="left",
                family="DejaVu Sans",
            )
            current_y -= LINE_STEP
        current_y -= PARAGRAPH_STEP
    last_line_bottom = current_y + LINE_STEP + PARAGRAPH_STEP - body_size / 72 / REPORT_PAGE_SIZE[1]
    return last_line_bottom - BLOCK_GAP

def _wrap_table_cell(value, width_hint, *, wrap_scale):
    max_chars = max(8, int(width_hint * wrap_scale))
    return "\n".join(
        textwrap.wrap(
            str(value),
            width=max_chars,
            break_long_words=True,
            break_on_hyphens=False,
        )
    )

def _table_row_weight(row, *, col_widths, wrap_scale):
    wrapped = [
        _wrap_table_cell(value, width, wrap_scale=wrap_scale)
        for value, width in zip(row, col_widths)
    ]
    return max(1.0, max(value.count("\n") + 1 for value in wrapped) * 0.82)

def _paginate_table_rows(rows, *, col_widths, wrap_scale=96, max_weight=18.0):
    pages = []
    current_page = []
    current_weight = 0.0
    for row in rows:
        row_weight = _table_row_weight(
            row,
            col_widths=col_widths,
            wrap_scale=wrap_scale,
        )
        if current_page and current_weight + row_weight > max_weight:
            pages.append(current_page)
            current_page = []
            current_weight = 0.0
        current_page.append(row)
        current_weight += row_weight
    if current_page:
        pages.append(current_page)
    return pages

def _draw_basic_table(ax, rows, *, columns, col_widths, bbox, font_size=REPORT_TABLE_SIZE, wrap_scale=96):
    wrapped_rows = [
        [
            _wrap_table_cell(value, width, wrap_scale=wrap_scale)
            for value, width in zip(row, col_widths)
        ]
        for row in rows
    ]
    wrapped_columns = [
        _wrap_table_cell(value, width, wrap_scale=wrap_scale)
        for value, width in zip(columns, col_widths)
    ]
    table = ax.table(
        cellText=wrapped_rows,
        colLabels=wrapped_columns,
        colWidths=col_widths,
        cellLoc="left",
        colLoc="left",
        bbox=bbox,
    )
    table.auto_set_font_size(False)
    table.set_fontsize(font_size)
    table.scale(1.0, 1.35)
    row_weights = [1.0]
    for row in wrapped_rows:
        row_weights.append(max(1.0, max(str(value).count("\n") + 1 for value in row) * 0.82))
    weight_total = sum(row_weights)
    for (row, col), cell in table.get_celld().items():
        cell.PAD = 0.095
        cell.set_height(bbox[3] * row_weights[row] / weight_total)
        cell.set_edgecolor("#D9E2EC")
        cell.set_linewidth(0.7)
        if row == 0:
            cell.set_facecolor(PUBLICATION_TABLE_HEADER)
            cell.set_text_props(weight="bold", color=PUBLICATION_BLUE)
        else:
            cell.set_facecolor("#FFFFFF" if row % 2 else PUBLICATION_TABLE_ALT)
    return table


def tool_means(summary: pd.DataFrame, tools: pd.DataFrame) -> pd.DataFrame:
    """Mean AUROC, Spearman and decoy z-score over all comparisons; one row per drawn tool, indexed by its label."""
    labels = drawn_tools(tools)
    rows = summary[(summary.subset == "all") & summary.predictor.isin(labels)]
    wide = rows.pivot(index="predictor", columns="metric", values="mean")
    return wide.loc[list(labels), ["auroc_d05", "spearman_shrunk", "decoy_z_shrunk"]].rename(index=labels)


def key_results(means: pd.DataFrame, responsive_pairs: pd.Series) -> list[str]:
    auroc, spearman = means.auroc_d05.idxmax(), means.spearman_shrunk.idxmax()
    covered, sparse = responsive_pairs.idxmax(), responsive_pairs.idxmin()
    return [f"Highest mean AUROC: {auroc} ({_format_summary_value(means.auroc_d05[auroc])}).",
            f"Highest mean Spearman: {spearman} ({_format_summary_value(means.spearman_shrunk[spearman])}).",
            f"Highest responsive pair coverage: {covered} ({_format_summary_value(responsive_pairs[covered], percent=True)}).",
            f"Lowest responsive pair coverage: {sparse} ({_format_summary_value(responsive_pairs[sparse], percent=True)})."]


def draw_inventory(pdf, rows: list[list[str]]) -> None:
    widths = [0.56, 0.44]
    pages = _paginate_table_rows(rows, col_widths=widths, wrap_scale=INVENTORY_WRAP_SCALE, max_weight=INVENTORY_ROWS)
    for page_rows in pages:
        fig, ax = _new_page()
        y = _bullet_block(ax, "Experiment inventory", [], y=TOP)
        weight = 1 + sum(_table_row_weight(row, col_widths=widths, wrap_scale=INVENTORY_WRAP_SCALE) for row in page_rows)
        height = min(y - BOTTOM, ROW_HEIGHT * weight)
        _draw_basic_table(ax, page_rows, columns=["Experiment", "Experiment ID"], col_widths=widths,
                          bbox=[LEFT, y - height, WIDTH, height], font_size=REPORT_TABLE_SIZE,
                          wrap_scale=INVENTORY_WRAP_SCALE)
        _save_page(pdf, fig)


def draw_plot(fig, page_ax, plots: pathlib.Path, stem: str, title: str, lines: list[str], y: float) -> float:
    """Draw a plot at full width under its title and one line per metric; returns where the next block starts."""
    page_w, page_h = REPORT_PAGE_SIZE
    y = _bullet_block(page_ax, title, lines, y=y)
    border = PLOT_BORDER_PIXELS
    image = plt.imread(plots / f"{stem}.png")[border:-border, border:-border]
    height = image.shape[0] / image.shape[1] * WIDTH * page_w / page_h
    ax = fig.add_axes([LEFT, y - height, WIDTH, height])
    ax.imshow(image, interpolation="antialiased")
    ax.axis("off")
    return y - height - BLOCK_GAP


def plot_tables(per: pd.DataFrame, coverage: pd.DataFrame, recovery: pd.DataFrame, precision: pd.DataFrame,
                means: pd.DataFrame, tools: pd.DataFrame) -> list[tuple]:
    """One (title, columns, column widths, rows) per plot, with the numbers that plot draws."""
    labels = drawn_tools(tools)
    tool_ids = list(labels)
    pct = lambda value: _format_summary_value(value, percent=True)
    c = coverage.set_index("predictor").loc[tool_ids]
    won_auroc, won_spearman = experiments_won(per, tool_ids, "auroc_d05"), experiments_won(per, tool_ids, "spearman_shrunk")
    budget = recovery.budget.max()
    recall = recovery[recovery.budget == budget].groupby("predictor").recall.mean()
    best = best_per_experiment(precision, labels).predictor.value_counts()
    return [
        ("Coverage", ["Predictor", "Gene coverage", "Experiments retained", "Responsive pairs", "No-response pairs"],
         [0.24, 0.19, 0.19, 0.19, 0.19],
         [[labels[t], pct(c.genes_covered[t] / c.genes_measured[t]), f"{c.experiments_retained[t]}/{c.experiments_total[t]}",
           pct(c.responsive_scored[t] / c.responsive_total[t]), pct(c.no_response_scored[t] / c.no_response_total[t])]
          for t in tool_ids]),
        ("Performance", ["Predictor", "Mean AUROC", "Mean Spearman", "Mean decoy z-score", "Won, AUROC", "Won, Spearman"],
         [0.22, 0.156, 0.156, 0.156, 0.156, 0.156],
         [[labels[t], _format_summary_value(means.auroc_d05[labels[t]]), _format_summary_value(means.spearman_shrunk[labels[t]]),
           _format_summary_value(means.decoy_z_shrunk[labels[t]]), str(won_auroc[t]), str(won_spearman[t])] for t in tool_ids]),
        ("Recovery", ["Predictor", f"Mean recall at {budget:,} genes", "Experiments won at 50% recall"],
         [0.24, 0.38, 0.38],
         [[labels[t], _format_summary_value(recall[t]), str(best.get(t, 0))] for t in tool_ids]),
    ]


def draw_tables(pdf, tables: list[tuple]) -> None:
    fig, ax = _new_page()
    y = TOP
    for title, columns, widths, rows in tables:
        height = ROW_HEIGHT * (len(rows) + 1)
        if y < TOP and y - 2 * BLOCK_GAP - height < BOTTOM:  # many tools: the title and table go to the next page
            _save_page(pdf, fig)
            fig, ax = _new_page()
            y = TOP
        y = _bullet_block(ax, title, [], y=y)
        _draw_basic_table(ax, rows, columns=columns, col_widths=widths, bbox=[LEFT, y - height, WIDTH, height],
                          font_size=REPORT_TABLE_SIZE, wrap_scale=70)
        y -= height + BLOCK_GAP
    _save_page(pdf, fig)


def write_report(path: pathlib.Path, *, comps: pd.DataFrame, set_aside: pd.DataFrame, summary: pd.DataFrame,
                 per: pd.DataFrame, coverage: pd.DataFrame, recovery: pd.DataFrame, precision: pd.DataFrame,
                 tools: pd.DataFrame, n_genes: int, min_per_state: int, plots: pathlib.Path) -> None:
    means = tool_means(summary, tools)
    pairs = coverage.set_index("predictor").loc[list(drawn_tools(tools))].rename(index=drawn_tools(tools))
    responsive_pairs = pairs.responsive_scored / pairs.responsive_total
    auroc_rows = summary[(summary.subset == "all") & (summary.metric == "auroc_d05")]
    n_auroc = int(auroc_rows.n_comparisons.max())
    rows = [[f"{c.mirna_name} · {c.cell_line}", c.experiment_id] for c in comps.itertuples(index=False)]
    rows += [[f"{c.mirna_name} · {c.cell_line} (set aside: {c.flags.replace(';', ', ').replace('_', ' ')})",
              c.experiment_id] for c in set_aside.itertuples(index=False)]
    with PdfPages(path) as pdf:
        fig, ax = _new_page()
        y = _header(ax, "fun-miRBench Benchmark Report")
        settings = [["Experiments evaluated", str(len(comps))],
                    ["Experiments skipped", str(len(set_aside))],
                    ["Predictors evaluated", str(len(means))],
                    ["Responsive genes", "oriented shrunken log2FC > 0; FSOS s-value < 0.05 against d = 0.5"],
                    ["No-response genes", "whole shrunken 95% interval of the oriented change below d = 0.5"]]
        settings_height = ROW_HEIGHT * (len(settings) + 1)
        _draw_basic_table(
            ax, settings,
            columns=["Run setting", "Value"], col_widths=[0.27, 0.73], bbox=[LEFT, y - settings_height, WIDTH, settings_height],
            font_size=REPORT_TABLE_SIZE)
        y = _bullet_block(
            ax, "Evaluation design",
            [f"Every predictor is scored on the same {n_genes:,} genes; unscored genes rank last.",
             "Unresolved genes are left out of the AUROC.",
             f"AUROC on {n_auroc} of {len(comps)} experiments (at least {min_per_state} responsive and "
             f"{min_per_state} no-response genes).",
             f"Spearman on all {len(comps)} experiments."],
            y=y - settings_height - BLOCK_GAP)
        y = _bullet_block(ax, "Key results", key_results(means, responsive_pairs), y=y)
        _bullet_block(ax, "Output tables", OUTPUT_TABLES, y=y)
        _save_page(pdf, fig)

        with plt.rc_context({"savefig.dpi": 300}):  # the plots are PNGs; the default would blur them in the PDF
            fig, ax = _new_page()
            y = draw_plot(fig, ax, plots, *COVERAGE, TOP)
            draw_plot(fig, ax, plots, *RECOVERY, y)
            _save_page(pdf, fig)

            fig, ax = _new_page()
            draw_plot(fig, ax, plots, *PERFORMANCE, TOP)
            _save_page(pdf, fig)

        draw_tables(pdf, plot_tables(per, coverage, recovery, precision, means, tools))

        draw_inventory(pdf, rows)
