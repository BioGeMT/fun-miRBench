"""Small, cached figures and standalone reports for vector exploration."""
from __future__ import annotations

import base64
import html
import json
import logging
import pathlib
import shutil
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D

from funmirbench.plots import tool_colors
from funmirbench.vector_sources import atomic_json

logger = logging.getLogger(__name__)
PLOT_VERSION = "explore-figures-3"


def read_results(directory):
    return {name: pd.read_csv(directory / f"{name}.tsv", sep="\t", keep_default_na=False,
                             na_values=["NA"]) for name in
            ("specificity", "candidate_alignment", "expression_bins", "comparisons", "predictors")}


def figures(data, directory, *, dpi=130, details=3):
    directory.mkdir(parents=True, exist_ok=True)
    summary, candidates, bins, comps, tools = (data[n] for n in
        ("specificity", "candidate_alignment", "expression_bins", "comparisons", "predictors"))
    ids = tools.tool_id.tolist()
    labels = dict(zip(ids, tools.get("official_name", tools.tool_id)))
    colors = tool_colors(ids)
    tick_labels = ["\n".join(textwrap.wrap(labels[t], 18)) for t in ids]
    specs = []

    def save(fig, stem, title, caption):
        fig.savefig(directory / f"{stem}.png", dpi=dpi, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        specs.append({"file": f"{stem}.png", "title": title, "caption": caption})

    fig, ax = plt.subplots(figsize=(10, 4.8), layout="constrained")
    for x, tool in enumerate(ids):
        sub = summary[summary.predictor == tool]
        pivot = sub.pivot(index="experiment_id", columns="metric", values="alignment")
        if not pivot.notna().any().any():
            ax.text(x, .04, "NA", transform=ax.get_xaxis_transform(), ha="center", color="grey")
        for _, row in pivot.iterrows():
            ax.plot([x - .13, x + .13], row.reindex(["spearman", "spearman_within_tenths"]),
                    color=colors[tool], alpha=.15, linewidth=.8)
        for dx, metric, marker in ((-.13, "spearman", "o"), (.13, "spearman_within_tenths", "s")):
            values = pivot[metric].dropna()
            ax.scatter(np.full(len(values), x + dx), values, color=colors[tool], alpha=.45, s=18, marker=marker)
            if len(values):
                ax.scatter(x + dx, values.median(), color="black", s=65, marker=marker, zorder=4)
    ax.axhline(0, color="grey", linewidth=.8)
    ax.set(xticks=range(len(ids)), xticklabels=tick_labels, ylabel="Rank-vector cosine / Spearman")
    ax.legend(handles=[Line2D([], [], marker="o", color="black", linestyle="", label="Overall"),
                       Line2D([], [], marker="s", color="black", linestyle="", label="Within expression bins")], frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    save(fig, "alignment", "Does agreement survive expression stratification?",
         "Small points are experiments; paired lines join each experiment's two estimates. Black symbols are descriptive medians. "
         "The within-bin measure is a gene-count-weighted average; it need not equal the overall correlation.")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), layout="constrained")
    for ax, metric, title in zip(axes, ("spearman", "spearman_within_tenths"), ("Overall", "Within expression bins")):
        valid = summary[(summary.metric == metric) & summary.alignment.notna() & summary.fraction_decoys_better.notna()]
        for tool in ids:
            sub = summary[(summary.predictor == tool) & (summary.metric == metric)]
            for flagged in (False, True):
                part = sub[(sub.action != "keep") == flagged]
                ax.scatter(part.alignment, part.fraction_decoys_better, label=labels[tool] if not flagged else None,
                           facecolors="none" if flagged else colors[tool], edgecolors=colors[tool], s=38, alpha=.8)
        ax.set(title=title, xlabel="Cognate alignment", ylabel="Fraction of decoys doing better", ylim=(-.03, 1.03))
        ax.axvline(0, color="grey", linewidth=.6)
        ax.spines[["top", "right"]].set_visible(False)
        if valid.empty:
            ax.text(.5, .5, "No comparable specificity results", transform=ax.transAxes, ha="center", color="grey")
    handles = [Line2D([], [], marker="o", color=colors[t], linestyle="", label=labels[t]) for t in ids]
    axes[0].legend(handles=handles, fontsize=8, frameon=False)
    save(fig, "specificity", "Is the agreement miRNA specific?",
         "Each point is an experiment and predictor. Hollow markers identify flagged comparisons. "
         "Higher alignment and a lower fraction of better decoys indicate stronger, more specific agreement. "
         "Panel-dependent decoy counts are in specificity.tsv; these fractions are not p-values.")

    for start in range(0, len(comps), 30):
        subset = comps.iloc[start:start + 30]
        fig, axes = plt.subplots(1, 2, figsize=(14, max(4, .28 * len(subset) + 1.6)), layout="constrained")
        plain = summary[summary.metric == "spearman"]
        for ax, value, title, cmap, limits in (
                (axes[0], "alignment", "Overall alignment", "coolwarm", (-1, 1)),
                (axes[1], "fraction_decoys_better", "Fraction of better decoys", "viridis_r", (0, 1))):
            wide = plain.pivot(index="experiment_id", columns="predictor", values=value).reindex(
                index=subset.experiment_id, columns=ids)
            cm = plt.get_cmap(cmap).copy()
            cm.set_bad("#D9D9D9")
            image = ax.imshow(np.ma.masked_invalid(wide.to_numpy(float)), aspect="auto", cmap=cm,
                              vmin=limits[0], vmax=limits[1])
            ax.set_xticks(range(len(ids)), tick_labels, fontsize=8, rotation=30, ha="right")
            row_labels = [f"[{start + i + 1:02d}] {r.study} {r.mirna_name.removeprefix('hsa-')}" +
                          (" *" if r.action != "keep" else "") if hasattr(r, "study") else r.experiment_id
                          for i, r in enumerate(subset.itertuples(index=False))]
            ax.set_yticks(range(len(subset)), row_labels, fontsize=8)
            ax.set_title(title)
            fig.colorbar(image, ax=ax, shrink=.7)
        save(fig, f"experiments_{start // 30 + 1}", "Which experiments carry signal?",
             "Bracketed numbers follow the experiment table's row order, disambiguating repeated study/miRNA labels. "
             "Asterisk = flagged; grey = unavailable, never zero. These are descriptive experiment-level results.")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), layout="constrained")
    plain = summary[summary.metric == "spearman"]
    for x, tool in enumerate(ids):
        part = plain[plain.predictor == tool]
        values = part.coverage.dropna() * 100
        axes[0].scatter(np.full(len(values), x), values, color=colors[tool], alpha=.5, s=20)
        if len(values):
            axes[0].plot([x - .2, x + .2], [values.median()] * 2, color="black", linewidth=2)
        for dx, metric in ((-.12, "spearman"), (.12, "spearman_within_tenths")):
            decoys = summary[(summary.predictor == tool) & (summary.metric == metric)].n_decoys
            axes[1].scatter(np.full(len(decoys), x + dx), decoys, color=colors[tool],
                            marker="o" if dx < 0 else "s", alpha=.6, s=20)
    for ax in axes:
        ax.set_xticks(range(len(ids)), tick_labels)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set(ylabel="Scored genes / usable-response genes (%)", ylim=(-2, 102))
    axes[1].set(ylabel="Comparable decoy miRNAs")
    save(fig, "coverage", "How much information supports each result?",
         "Left: prediction coverage among genes with a usable response; black marks show medians. "
         "Right: decoy counts for overall (circles) and within-bin (squares) specificity. Coverage is not accuracy.")

    fig, ax = plt.subplots(figsize=(10, 4.5), layout="constrained")
    for tool in ids:
        group = bins[bins.predictor == tool].groupby("expression_tenth").spearman.median()
        ax.plot(group.index, group.values, "o-", color=colors[tool], label=labels[tool])
    if not bins.spearman.notna().any():
        ax.text(.5, .5, "No usable expression-bin correlations", transform=ax.transAxes, ha="center", color="grey")
    ax.set(xlabel="Control-expression tenth (low to high)", ylabel="Median per-bin Spearman", xticks=range(1, 11))
    ax.axhline(0, color="grey", linewidth=.8)
    ax.legend(frameon=False, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    save(fig, "expression_bins", "Where does alignment occur along baseline expression?",
         "Each point is the median of defined experiment-level correlations in that bin. Missing bins are unavailable. "
         "This is not a pooled-gene correlation or a confidence interval; expression_bins.tsv contains every value and gene count.")

    for number, comparison in enumerate(comps.head(details).itertuples(index=False)):
        fig, ax = plt.subplots(figsize=(10, 4.5), layout="constrained")
        for x, tool in enumerate(ids):
            sub = candidates[(candidates.experiment_id == comparison.experiment_id) & (candidates.predictor == tool)]
            own = sub[sub.is_cognate].iloc[0]
            values = sub.loc[(sub.candidate_seed != own.candidate_seed), "spearman"].dropna().sort_values()
            # Bounded plotting cost even for large registries; numerical results always retain all candidates.
            if len(values) > 2000:
                values = values.iloc[np.linspace(0, len(values) - 1, 2000).astype(int)]
            ax.scatter(x + np.linspace(-.16, .16, len(values)), values, color=colors[tool], alpha=.5, s=13)
            if np.isfinite(own.spearman):
                ax.scatter(x, own.spearman, marker="D", color="black", s=65, zorder=4)
        ax.set(xticks=range(len(ids)), xticklabels=tick_labels, ylabel="Overall Spearman")
        ax.axhline(0, color="grey", linewidth=.7)
        ax.spines[["top", "right"]].set_visible(False)
        save(fig, f"decoys_{number + 1}", f"Cognate and decoys: {comparison.experiment_id}",
             "Colored points are different-seed decoys; black diamonds are the cognate. "
             "Detail examples are the first comparisons in registry order, not selected by performance. "
             "At most 2,000 evenly spaced sorted decoys per tool are drawn; full values remain in candidate_alignment.tsv.")
    atomic_json(directory / "figures.json", specs)
    return specs


def prepare_figures(data, result_cache, output, *, dpi=130, details=3):
    directory = result_cache / f"{PLOT_VERSION}_dpi{dpi}_details{details}"
    manifest = directory / "figures.json"
    if manifest.is_file():
        specs = json.loads(manifest.read_text())
        if not all((directory / s["file"]).is_file() for s in specs):
            specs = figures(data, directory, dpi=dpi, details=details)
        else:
            logger.info("reusing %d cached figures", len(specs))
    else:
        specs = figures(data, directory, dpi=dpi, details=details)
    output.mkdir(parents=True, exist_ok=True)
    for spec in specs:
        shutil.copy2(directory / spec["file"], output / spec["file"])
    atomic_json(output / "figures.json", specs)
    return specs


METHODS = [
    "Observed response: negative shrunken log2FC for OE, positive for repression. Both response and predictor scores "
    "are converted to average ranks, centered and normalized. Their cosine equals Spearman; score magnitudes are not predicted fold changes.",
    "Unusable responses are excluded; all unscored genes form a bottom tie. Overall and within-bin alignment are shown separately. "
    "Expression ties stay together, unknown baselines have no bin, and bins require at least 20 usable genes.",
    "Specificity compares the cognate with finite different-seed registry candidates (nt 2-8). Ties count half. "
    "Within-bin decoys must use the same contributing bins as the cognate. Z-scores need at least two decoys and nonzero SD.",
    "Decoy ranks, fractions and z-scores are descriptive, not p-values. Families may contribute several members. "
    "Different predictors can have different decoy panels. Shared-control comparisons are not independent replicates.",
    "Flagged experiments are retained and identified. Results summarize a published data snapshot, not a new run of the final RNA-seq pipeline. "
    "Direct binding cannot be established from transcriptomic agreement alone.",
]
PENDING = ("Replicate stability and controls-only contrasts require sample-level inputs and additional fits. "
           "The seed-site baseline needs matching UTR/site annotations and a specified scoring rule. "
           "These checks are not computed or displayed as results here. The final metric still needs agreement with Steph.")


def write_reports(directory, specs, provenance):
    data = read_results(directory)
    summary, comps = data["specificity"], data["comparisons"]
    grouped = pd.concat([summary.assign(subset="all"), summary[summary.action == "keep"].assign(subset="unflagged")])
    medians = grouped.groupby(["subset", "predictor", "metric"], sort=False).agg(
        median_alignment=("alignment", "median"), valid_experiments=("alignment", "count"),
        median_fraction_decoys_better=("fraction_decoys_better", "median"),
        median_decoy_count=("n_decoys", "median"), median_coverage=("coverage", "median")).reset_index()
    medians.to_csv(directory / "overview.tsv", sep="\t", index=False, na_rep="NA")
    intro = (f"Exploratory vector analysis of {len(comps)} OE comparisons in {comps.dependency_group.nunique()} dependency groups; "
             f"{int((comps.action != 'keep').sum())} flagged comparisons retained. "
             f"Gene universe: {int(summary.n_genes.iloc[0]):,}. "
             f"Experimental tables: Zenodo record {provenance['benchmark_record']}.")
    prediction_sources = "; ".join(f"{name}: " + (f"Zenodo {source['record']}" if "record" in source else source["origin"])
                                   for name, source in provenance["prediction_sources"].items())
    esc = html.escape
    images = []
    for spec in specs:
        encoded = base64.b64encode((directory / "plots" / spec["file"]).read_bytes()).decode()
        images.append(f"<section><h2>{esc(spec['title'])}</h2><img alt='{esc(spec['title'])}' "
                      f"src='data:image/png;base64,{encoded}'><p>{esc(spec['caption'])}</p></section>")
    methods = "".join(f"<p>{esc(line)}</p>" for line in METHODS)
    tables = []
    columns = ["predictor", "metric", "alignment", "coverage", "n_usable", "n_decoys", "cognate_rank",
               "fraction_decoys_better", "decoy_z", "n_decoys_excluded_support", "alignment_na_reason", "specificity_na_reason"]
    for i, c in enumerate(comps.itertuples(index=False)):
        rows = summary[summary.experiment_id == c.experiment_id]
        tables.append(f"<details><summary>[{i + 1:02d}] {esc(c.experiment_id)} — {esc(c.action)}</summary>" +
                      rows[columns].to_html(index=False, na_rep="NA", float_format=lambda x: f"{x:.4g}") + "</details>")
    page = f"""<!doctype html><html lang="en"><meta charset="utf-8"><title>Exploratory response vectors</title>
<style>body{{font:16px/1.5 Arial,sans-serif;color:#22303c;max-width:1200px;margin:32px auto;padding:0 24px}}
h1,h2{{color:#17324d}}nav a{{margin-right:20px}}section{{margin:32px 0}}img{{width:100%;height:auto}}
table{{border-collapse:collapse;font-size:13px;display:block;overflow:auto}}th,td{{border:1px solid #d9d9d9;padding:8px}}
th{{background:#e9f1fb}}tr:nth-child(even){{background:#f9fbfd}}details{{margin:15px 0}}summary{{cursor:pointer;font-weight:bold}}</style>
<h1>Exploratory response vectors</h1><p>{esc(intro)}</p><p>{esc(prediction_sources)}</p>
<nav><a href="#overview">Overview</a><a href="#figures">Figures</a><a href="#experiments">Experiments</a><a href="#methods">Methods and limits</a></nav>
<section id="overview"><h2>Descriptive overview</h2>{medians.to_html(index=False, na_rep='NA', float_format=lambda x: f'{x:.4g}')}</section>
<div id="figures">{''.join(images)}</div><section id="experiments"><h2>Per experiment values and missing results</h2>{''.join(tables)}</section>
<section id="methods"><h2>Methods and interpretation</h2>{methods}<h2>Checks still pending</h2><p>{esc(PENDING)}</p></section></html>"""
    (directory / "REPORT.html").write_text(page, encoding="utf-8")
    temporary_pdf = directory / ".REPORT.pdf.tmp"
    with PdfPages(temporary_pdf) as pdf:
        fig = plt.figure(figsize=(8.27, 11.69))
        fig.text(.08, .93, "Exploratory response vectors", fontsize=19, weight="bold", color="#17324D")
        y = .87
        for title, text in [("Scope", intro), ("Prediction sources", prediction_sources),
                            *[("Method" if i == 0 else "", line) for i, line in enumerate(METHODS)],
                            ("Still pending", PENDING)]:
            if title:
                fig.text(.08, y, title, fontsize=11, weight="bold")
                y -= .025
            lines = textwrap.wrap(text, 93)
            fig.text(.08, y, "\n".join(lines), va="top", fontsize=9.2, linespacing=1.35)
            y -= len(lines) * .017 + .028
        pdf.savefig(fig); plt.close(fig)
        for metric, title in (("spearman", "Overall"), ("spearman_within_tenths", "Within expression bins")):
            rows = medians[(medians.subset == "all") & (medians.metric == metric)]
            fig, ax = plt.subplots(figsize=(11.69, 8.27))
            ax.axis("off")
            ax.set_title(f"Descriptive overview — {title}", loc="left", fontsize=17, weight="bold", color="#17324D")
            table = rows[["predictor", "median_alignment", "valid_experiments", "median_fraction_decoys_better",
                          "median_decoy_count", "median_coverage"]].copy()
            for col in table.columns[1:]:
                table[col] = table[col].map(lambda x: "NA" if pd.isna(x) else f"{x:.3g}")
            drawn = ax.table(cellText=table.values, colLabels=["Predictor", "Median\nalignment", "Valid\nexperiments",
                             "Median fraction\nbetter decoys", "Median\ndecoy count", "Median\ncoverage (0–1)"],
                             colWidths=[.26, .14, .12, .19, .14, .15], cellLoc="left", bbox=[0, .35, 1, .5])
            drawn.auto_set_font_size(False); drawn.set_fontsize(9)
            for (r, _), cell in drawn.get_celld().items():
                cell.set_edgecolor("#D9D9D9")
                if r == 0:
                    cell.set_facecolor("#E9F1FB"); cell.set_text_props(weight="bold")
            ax.text(0, .24, "Medians over defined experiment-level values; no independence or significance is implied.\n"
                    "All and unflagged subsets, exact per-experiment values and NA reasons are in the HTML report and TSVs.",
                    fontsize=10, va="top")
            pdf.savefig(fig); plt.close(fig)
        for spec in specs:
            fig = plt.figure(figsize=(11.69, 8.27))
            title = "\n".join(textwrap.wrap(spec["title"], 95))
            fig.text(.055, .95, title, fontsize=15, weight="bold", va="top", color="#17324D")
            image = plt.imread(directory / "plots" / spec["file"])
            ax = fig.add_axes([.04, .15, .92, .70])
            ax.imshow(image); ax.axis("off")
            fig.text(.055, .11, "\n".join(textwrap.wrap(spec["caption"], 145)), va="top", fontsize=9)
            pdf.savefig(fig); plt.close(fig)
    with temporary_pdf.open("rb") as handle:
        handle.seek(max(0, temporary_pdf.stat().st_size - 1024))
        if not handle.read().rstrip().endswith(b"%%EOF"):
            raise OSError("PDF did not finalize correctly")
    temporary_pdf.replace(directory / "REPORT.pdf")
    return medians
