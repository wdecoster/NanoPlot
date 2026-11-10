"""Tests for NanoPlot --percentqual handling of quality-less input (e.g. fasta)."""
import pandas as pd

from nanoplot.filteroptions import filter_and_transform_data


def _settings(**overrides):
    settings = {
        "alength": False,
        "bam": False,
        "drop_outliers": False,
        "maxlength": None,
        "minlength": None,
        "minqual": None,
        "loglength": False,
        "runtime_until": None,
        "downsample": None,
        "percentqual": True,
    }
    settings.update(overrides)
    return settings


def test_percentqual_without_quality_column_is_skipped():
    """--percentqual on data without qualities must not crash with KeyError."""
    df = pd.DataFrame({"lengths": [100, 200, 300]})
    out_df, _ = filter_and_transform_data(df, _settings())
    assert list(out_df["lengths"]) == [100, 200, 300]


def test_percentqual_with_quality_column_converts():
    """--percentqual on data with qualities converts phred to percent identity."""
    df = pd.DataFrame({"lengths": [100, 200], "quals": [20, 30]})
    out_df, _ = filter_and_transform_data(df, _settings())
    assert list(out_df["quals"]) == [99.0, 99.9]
