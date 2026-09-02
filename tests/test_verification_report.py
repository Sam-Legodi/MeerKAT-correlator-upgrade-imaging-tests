import numpy as np
import pytest

from meerkat_corr_imaging.verification_report import _correction_state, fit_rigid_transform


def test_fit_rigid_transform_recovers_rotation_and_translation():
    rng = np.random.default_rng(12)
    reference = rng.normal(size=(80, 2)) * 500.0
    angle_deg = 0.013
    angle_rad = np.deg2rad(angle_deg)
    rotation = np.array(
        [[np.cos(angle_rad), -np.sin(angle_rad)], [np.sin(angle_rad), np.cos(angle_rad)]]
    )
    translation = np.array([0.42, -0.27])
    test = (rotation @ reference.T).T + translation

    result = fit_rigid_transform(reference, test, bootstrap_samples=20)

    assert result.rotation_deg == pytest.approx(angle_deg, abs=1e-8)
    assert result.translation_east_arcsec == pytest.approx(translation[0], abs=1e-8)
    assert result.translation_north_arcsec == pytest.approx(translation[1], abs=1e-8)
    assert result.postfit_p95_arcsec < 1e-8


def test_fit_rigid_transform_rejects_too_few_matches():
    reference = np.zeros((9, 2))
    with pytest.raises(ValueError, match="at least 10"):
        fit_rigid_transform(reference, reference)


@pytest.mark.parametrize(
    ("header_value", "expected"),
    [
        (False, "non-PB"),
        ("F", "non-PB"),
        ("false", "non-PB"),
        ("0", "non-PB"),
        (True, "PB-corrected"),
        ("T", "PB-corrected"),
        ("true", "PB-corrected"),
        ("1", "PB-corrected"),
    ],
)
def test_correction_state_parses_fits_boolean_representations(
    monkeypatch, header_value, expected
):
    monkeypatch.setattr(
        "meerkat_corr_imaging.verification_report.fits.getheader",
        lambda path: {"PBCOR": header_value},
    )

    assert _correction_state("image.fits") == expected


@pytest.mark.parametrize(
    ("filename", "expected"),
    [
        ("image_IClean_PB.fits", "PB-corrected"),
        ("image_non_pb.fits", "non-PB"),
        ("image_lowband.fits", "non-PB"),
    ],
)
def test_correction_state_filename_fallback(monkeypatch, filename, expected):
    monkeypatch.setattr(
        "meerkat_corr_imaging.verification_report.fits.getheader",
        lambda path: {},
    )

    assert _correction_state(filename) == expected
