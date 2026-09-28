import numpy as np
import pandas as pd

from meerkat_corr_imaging.reference_comparison import (
    degradation_decision, parity_decision, position_reference_test,
    visibility_exposure,
)
from meerkat_corr_imaging.uncertainty import UncertaintyConfig


def test_confidence_intervals_use_reference_parity():
    assert parity_decision({"ci_low": .97, "ci_high": 1.01})["status"] == "Pass"
    assert parity_decision({"ci_low": 1.02, "ci_high": 1.08})["status"] == "Concern"
    assert degradation_decision({"ci_low": -.03, "ci_high": .08})["status"] == "Pass"
    assert degradation_decision({"ci_low": .01, "ci_high": .08})["status"] == "Concern"


def test_position_compares_with_catalogue_noise_and_joint_shift():
    covariance = np.tile(np.diag([.1**2, .1**2]), (30, 1, 1))
    config = UncertaintyConfig(200, .95, 42)
    normal = position_reference_test(.25, covariance, [0., 0.], np.diag([.1, .1]), config)
    shifted = position_reference_test(.25, covariance, [.6, .0], np.diag([.01, .01]), config)
    extended = position_reference_test(1.0, covariance, None, None, config)
    assert normal["status"] == "Pass"
    assert shifted["translation_concern"] and shifted["status"] == "Concern"
    assert extended["p95_concern"] and extended["status"] == "Concern"


def test_visibility_exposure_uses_unflagged_cross_samples(tmp_path):
    pd.DataFrame([dict(CLASS="cross", POL="XX", SCAN=1, FLAGGED=10, TOTAL=100),
                  dict(CLASS="cross", POL="XX", SCAN=2, FLAGGED=20, TOTAL=100),
                  dict(CLASS="auto", POL="XX", SCAN=1, FLAGGED=0, TOTAL=100)]).to_csv(
        tmp_path / "flagging_vs_scan.csv", index=False)
    pd.DataFrame({"FREQ_HZ": [100., 200., 300.]}).to_csv(tmp_path / "rfi_free_channel_mask.csv", index=False)
    pd.DataFrame({"TIME": [0., 2., 4.]}).to_csv(tmp_path / "perrow_amp_stats.csv", index=False)
    result = visibility_exposure(tmp_path)
    assert result["unflagged_cross_samples"] == 170
    assert result["effective_exposure"] == 34_000
