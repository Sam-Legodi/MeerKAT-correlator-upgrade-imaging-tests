from __future__ import annotations

import numpy as np
import pandas as pd

from meerkat_corr_imaging.vis_amp_analyze import (
    acceptance_summary,
    derive_rfi_free_channel_mask,
    spectrum_metrics,
)


def test_rfi_mask_uses_strict_twenty_percent_limit() -> None:
    flagged = np.array([0, 1, 2, 3])
    total = np.array([0, 10, 10, 10])

    mask = derive_rfi_free_channel_mask(flagged, total, threshold=0.20)

    assert mask.tolist() == [False, True, False, False]


def test_spectrum_metrics_exclude_masked_outlier() -> None:
    clean = np.linspace(10.0, 11.0, 101)
    amplitude = clean.copy()
    amplitude[50] = 1_000_000.0
    valid = np.ones(amplitude.shape, dtype=bool)
    valid[50] = False

    metrics = spectrum_metrics(amplitude, valid, window=51, order=3)

    assert metrics["N_VALID"] == 100
    assert np.isclose(metrics["MEAN"], clean[valid].mean())
    assert metrics["OSC_FRAC"] < 1e-10


def test_spectrum_metrics_detect_one_percent_oscillation_concern() -> None:
    channel = np.arange(101)
    amplitude = 10.0 + 0.2 * np.where(channel % 2, 1.0, -1.0)

    metrics = spectrum_metrics(
        amplitude, np.ones(amplitude.shape, dtype=bool), window=51, order=3
    )

    assert metrics["OSC_FRAC"] > 0.01


def test_acceptance_summary_applies_strict_limits_and_worst_spectrum() -> None:
    scan_stats = pd.DataFrame(
        [
            {"CLASS": "auto", "POL": "XX", "OSC_FRAC": 0.005},
            {"CLASS": "auto", "POL": "XX", "OSC_FRAC": 0.009},
            {"CLASS": "cross", "POL": "XX", "OSC_FRAC": 0.010},
        ]
    )
    flagging = pd.DataFrame(
        [
            {"CLASS": "auto", "POL": "XX", "FLAGGED": 1, "TOTAL": 10},
            {"CLASS": "cross", "POL": "XX", "FLAGGED": 2, "TOTAL": 10},
        ]
    )

    scan_stats["BASE"] = ["0-0", "0-0", "0-1"]
    scan_stats["SCAN"] = [1, 2, 1]

    summary = acceptance_summary(scan_stats, flagging)
    auto = summary[(summary.CLASS == "auto") & (summary.POL == "XX")].iloc[0]
    cross = summary[(summary.CLASS == "cross") & (summary.POL == "XX")].iloc[0]

    assert auto.FLAG_STATUS == "Pass"
    assert auto.OSC_STATUS == "Pass"
    assert cross.FLAG_STATUS == "Concern"
    assert cross.OSC_STATUS == "Concern"
    assert cross.OSC_MAX_PCT == 1.0


def test_oscillation_eligibility_and_robust_scatter():
    values = 10 + 0.001 * np.sin(np.arange(401))
    valid = np.ones(401, dtype=bool)
    clean = spectrum_metrics(values, valid)["OSC_FRAC"]
    values[200] = 10000
    assert spectrum_metrics(values, valid)["OSC_FRAC"] < 0.01
    assert clean < 0.01
    assert np.isnan(spectrum_metrics(values, valid, flag_fraction=0.20)["OSC_FRAC"])
    assert np.isnan(spectrum_metrics(values, valid, flag_fraction=np.nan)["OSC_FRAC"])
    assert np.isnan(spectrum_metrics(np.zeros(401), valid)["OSC_FRAC"])
    valid[:] = False
    valid[:50] = True
    assert np.isfinite(spectrum_metrics(values, valid, expected_channels=100)["OSC_FRAC"])
    valid[49] = False
    assert np.isnan(spectrum_metrics(values, valid, expected_channels=100)["OSC_FRAC"])


def test_baseline_medians_resist_one_scan_outlier_and_keep_warning():
    frame = pd.DataFrame([
        dict(CLASS="cross", POL="XX", BASE=f"0-{b}", SCAN=s, OSC_FRAC=.005)
        for b in range(1, 11) for s in range(3)
    ])
    frame.loc[0, "OSC_FRAC"] = 10.0
    flags = pd.DataFrame([dict(CLASS="cross", POL="XX", FLAGGED=1, TOTAL=100)])
    result = acceptance_summary(frame, flags).iloc[0]
    assert result.OSC_STATUS == "Pass"
    assert result.OSC_OUTLIER_WARNING
    assert result.OSC_MAX_FRAC == 10.0
    assert result.OSC_BASELINE_P95_FRAC == .005
    frame.loc[1, "OSC_FRAC"] = 10.0
    assert acceptance_summary(frame, flags).iloc[0].OSC_STATUS == "Concern"
    frame["OSC_FRAC"] = np.nan
    result = acceptance_summary(frame, flags).iloc[0]
    assert result.OSC_STATUS == "Not assessed"
    assert result.N_UNASSESSED_SPECTRA == 30


def test_single_panel_skips_absent_classes_and_nonfinite_products(tmp_path, monkeypatch):
    from meerkat_corr_imaging import vis_amp_analyze as vis
    frame = pd.DataFrame([
        {'CLASS': 'cross', 'POL': 'XX', 'CHAN': 1, 'FLAG_FRAC': .1, 'TOTAL': 10},
        {'CLASS': 'cross', 'POL': 'YY', 'CHAN': 1, 'FLAG_FRAC': .2, 'TOTAL': 10},
        {'CLASS': 'cross', 'POL': 'XY', 'CHAN': 1, 'FLAG_FRAC': np.nan, 'TOTAL': 0},
    ])
    captured = []
    original = vis.plt.close
    def capture(fig):
        if hasattr(fig, 'axes'):
            captured.append([(axis.get_title(), [line.get_label() for line in axis.lines]) for axis in fig.axes])
        original(fig)
    monkeypatch.setattr(vis.plt, 'close', capture)
    paths = vis._plot_two_classes(frame, 'CHAN', 'FLAG_FRAC', 'Flagging', 'Channel', 'Fraction', tmp_path / 'flagging_by_channel.png')
    assert len(paths) == 1
    assert paths[0].name == 'flagging_by_channel.png'
    assert captured[-1] == [('Flagging (cross baselines)', ['XX', 'YY'])]
    frame.loc[len(frame)] = ['auto', 'XX', 1, .1, 10]
    assert len(vis._plot_two_classes(frame, 'CHAN', 'FLAG_FRAC', 'Flagging', 'Channel', 'Fraction', tmp_path / 'both.png')) == 2
    assert all(len(axes) == 1 for axes in captured)


def test_physical_baseline_lengths():
    from meerkat_corr_imaging.vis_amp_analyze import _add_baseline_lengths
    frame = pd.DataFrame({'ANT1': [0, 0, 1], 'ANT2': [0, 1, 0]})
    _add_baseline_lengths(frame, np.array([[0, 0, 0], [3, 4, 0]]))
    assert frame.BASELINE_LENGTH_M.tolist() == [0, 5, 5]
