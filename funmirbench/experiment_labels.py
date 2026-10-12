"""Human-readable experiment labels for reports and plots."""

from __future__ import annotations


_MISSING_VALUES = {"", "NA", "N/A", "NONE", "NAN"}


def _display_value(value) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    return "" if text.upper() in _MISSING_VALUES else text


def format_experiment_label(
    *,
    geo_accession=None,
    mirna=None,
    cell_line=None,
    perturbation=None,
    separator=" · ",
) -> str:
    """Return a compact human-readable experiment label."""
    parts = [
        _display_value(geo_accession),
        _display_value(mirna),
        _display_value(cell_line),
        _display_value(perturbation),
    ]
    return separator.join(part for part in parts if part)


def format_experiment_id_label(dataset_id) -> str:
    """Return the stable machine identifier as a provenance label."""
    return f"Experiment ID: {dataset_id}"
