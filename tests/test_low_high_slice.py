from __future__ import annotations

from argparse import Namespace
from pathlib import Path

import numpy as np
from astropy.io import fits

from meerkat_corr_imaging.config import _dict_to_dataclass, load_config
from meerkat_corr_imaging.low_high_slice import Subband, extract_slices, select_subband
from meerkat_corr_imaging.steps import step4_low_high_slice


def write_cuboid(path: Path) -> np.ndarray:
    data = np.zeros((1, 5, 4, 3), dtype=np.float32)
    for plane_index in range(data.shape[1]):
        data[0, plane_index] = plane_index

    header = fits.Header()
    header["CTYPE1"] = "RA---SIN"
    header["CTYPE2"] = "DEC--SIN"
    header["CTYPE3"] = "SPECLNMF"
    header["CTYPE4"] = "STOKES"
    header["CRVAL1"] = 10.0
    header["CRVAL2"] = -30.0
    header["CRVAL3"] = 200.0
    header["CRVAL4"] = 1.0
    header["CRPIX1"] = 2.0
    header["CRPIX2"] = 2.0
    header["CRPIX3"] = 1.0
    header["CRPIX4"] = 1.0
    header["CDELT1"] = -0.01
    header["CDELT2"] = 0.01
    header["CDELT3"] = 1.0
    header["CDELT4"] = 1.0
    header["CUNIT1"] = "deg"
    header["CUNIT2"] = "deg"
    header["NTERM"] = 2
    header["NSPEC"] = 3
    header["RFALPHA"] = 200.0
    header["BUNIT"] = "Jy/beam"
    header["CLEANBMJ"] = 0.01
    header["CLEANBMN"] = 0.008
    header["CLEANBPA"] = 12.0
    for number, (low, effective, high) in enumerate(
        ((100.0, 130.0, 160.0), (160.0, 200.0, 240.0), (240.0, 270.0, 300.0)),
        start=1,
    ):
        header[f"FREL{number:04d}"] = low
        header[f"FEFF{number:04d}"] = effective
        header[f"FREH{number:04d}"] = high
    fits.PrimaryHDU(data=data, header=header).writeto(path)
    return data


def test_selection_maximises_overlap_and_uses_directional_ties() -> None:
    subbands = [
        Subband(1, 2, 100.0, 130.0, 160.0),
        Subband(2, 3, 160.0, 200.0, 240.0),
        Subband(3, 4, 240.0, 270.0, 300.0),
    ]
    assert select_subband(subbands, (110.0, 220.0), "lowband").subband.number == 2
    assert select_subband(subbands, (200.0, 290.0), "highband").subband.number == 3

    tied_low = select_subband(subbands, (100.0, 300.0), "lowband")
    tied_high = select_subband(subbands, (100.0, 300.0), "highband")
    assert tied_low.subband.number == 2  # largest full coverage (80 Hz)
    assert tied_high.subband.number == 2

    equal_width = [
        Subband(1, 2, 100.0, 130.0, 160.0),
        Subband(2, 3, 160.0, 190.0, 220.0),
    ]
    assert select_subband(equal_width, (100.0, 220.0), "lowband").subband.number == 1
    assert select_subband(equal_width, (100.0, 220.0), "highband").subband.number == 2


def test_enabled_config_requires_slice_frequency_ranges() -> None:
    try:
        _dict_to_dataclass(
            {
                "project_name": "test",
                "low_high_slice": {"enabled": True},
            }
        )
    except ValueError as exc:
        assert "lowband_hz" in str(exc)
        assert "highband_hz" in str(exc)
    else:
        raise AssertionError("Enabled low_high_slice accepted missing frequency ranges")


def test_extracts_single_2d_planes_with_frequency_and_beam_metadata(tmp_path: Path) -> None:
    cuboid = tmp_path / "observation_IClean.fits"
    source_data = write_cuboid(cuboid)

    low_path, high_path = extract_slices(
        cuboid,
        lowband_hz=(110.0, 220.0),
        highband_hz=(200.0, 290.0),
    )

    with fits.open(low_path) as low_hdul:
        assert low_hdul[0].data.shape == (4, 3)
        np.testing.assert_array_equal(low_hdul[0].data, source_data[0, 3])
        assert low_hdul[0].header["SUBBAND"] == 2
        assert low_hdul[0].header["CUBPLANE"] == 3
        assert low_hdul[0].header["RESTFRQ"] == 200.0
        assert low_hdul[0].header["PBCOR"] is False
        assert low_hdul[0].header["BMAJ"] == 0.01
        assert "CTYPE3" not in low_hdul[0].header

    with fits.open(high_path) as high_hdul:
        np.testing.assert_array_equal(high_hdul[0].data, source_data[0, 4])
        assert high_hdul[0].header["SUBBAND"] == 3
        assert high_hdul[0].header["RESTFRQ"] == 270.0


def test_pipeline_step_creates_default_outputs_and_resolves_source_images(
    tmp_path: Path,
    monkeypatch,
) -> None:
    cuboid = tmp_path / "test_IClean.fits"
    write_cuboid(cuboid)
    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        f"""
project_name: test
paths:
  raw_dir: {tmp_path / 'raw'}
  interim_dir: {tmp_path / 'interim'}
  processed_dir: {tmp_path / 'processed'}
  reports_dir: {tmp_path / 'reports'}
  sky_xmatches_dir: {tmp_path / 'xmatches'}
reference: {{name: none, ms_paths: [], images: []}}
tests:
  - name: test_target
    ms_paths: []
    images: []
    cuboid: {cuboid}
frequency_ranges:
  lowband_hz: [110, 220]
  highband_hz: [200, 290]
casa: {{}}
low_high_slice:
  enabled: true
"""
    )
    cfg = load_config(config_path)

    # The packaged extractor must remain runnable outside the repository root.
    monkeypatch.chdir(tmp_path)
    step4_low_high_slice.run(cfg)

    expected = [
        str((tmp_path / "test_IClean_lowband.fits").resolve()),
        str((tmp_path / "test_IClean_highband.fits").resolve()),
    ]
    assert step4_low_high_slice.source_images(cfg) == expected

    # Existing images are reused even when the cuboid is no longer available.
    cuboid.unlink()
    step4_low_high_slice.run(cfg)


def test_pybdsf_global_frequency_is_only_a_fallback(tmp_path: Path) -> None:
    from meerkat_corr_imaging.pybdsf_srcfind import determine_frequency

    image = tmp_path / "image.fits"
    header = fits.Header()
    header["RESTFRQ"] = 123.0
    fits.PrimaryHDU(data=np.ones((2, 2), dtype=np.float32), header=header).writeto(image)
    args = Namespace(freq_hz=999.0, freq_mhz=None)

    assert determine_frequency(str(image.resolve()), {}, args) is None
    assert determine_frequency(str(image.resolve()), {str(image.resolve()): 456.0}, args) == 456.0


def test_legacy_slice_config_is_normalised_without_duplicate_runtime_state() -> None:
    cfg = _dict_to_dataclass(
        {
            "project_name": "legacy",
            "reference": {"name": "reference"},
            "tests": [{"name": "test"}],
            "casa": {
                "lowband_hz": [100.0, 200.0],
                "highband_hz": [300.0, 400.0],
            },
            "low_high_slice": {
                "enabled": True,
                "lowband_hz": [100.0, 200.0],
                "highband_hz": [300.0, 400.0],
                "reference": {"name": "reference", "cuboid": "reference.fits"},
                "tests": [{"name": "test", "cuboid": "test.fits"}],
            },
        }
    )

    assert cfg.frequency_ranges.lowband_hz == [100.0, 200.0]
    assert cfg.frequency_ranges.highband_hz == [300.0, 400.0]
    assert cfg.casa.lowband_hz is None
    assert cfg.low_high_slice.lowband_hz is None
    assert cfg.reference.cuboid == "reference.fits"
    assert cfg.tests[0].cuboid == "test.fits"


def test_pybdsf_uses_unambiguous_frequency_parameter(tmp_path: Path, monkeypatch) -> None:
    import sys
    from types import SimpleNamespace

    from meerkat_corr_imaging.pybdsf_srcfind import find_sources

    image = tmp_path / "image.fits"
    header = fits.Header()
    header["RESTFRQ"] = 123.0e6
    fits.PrimaryHDU(data=np.ones((2, 2), dtype=np.float32), header=header).writeto(image)

    captured: dict[str, object] = {}

    class FakeImage:
        def export_image(self, **kwargs) -> None:
            pass

        def write_catalog(self, **kwargs) -> None:
            pass

    def process_image(input_name: str, **kwargs) -> FakeImage:
        captured["input"] = input_name
        captured["kwargs"] = kwargs
        return FakeImage()

    monkeypatch.setitem(sys.modules, "bdsf", SimpleNamespace(process_image=process_image))
    find_sources(str(image), base="image", adaptive_rms_box=False)

    kwargs = captured["kwargs"]
    assert isinstance(kwargs, dict)
    assert kwargs["frequency"] == 123.0e6
    assert "freq" not in kwargs
    assert kwargs["adaptive_rms_box"] is False
