from __future__ import annotations

from astropy.table import Table

from meerkat_corr_imaging.xmatch_pybdsf import sanitize_fits_meta


def _history_text(table: Table) -> str:
    history = table.meta.get("HISTORY", [])
    if isinstance(history, str):
        history = [history]
    return "".join(history)


def test_sanitize_fits_meta_keeps_complete_card_that_fits(tmp_path) -> None:
    table = Table({"source_id": [1]})
    keyword = "HIERARCH T1_INIMAGE"
    value = "1785941476_continuum_image_J2147-8132_IClean_lowband.fits"
    table.meta[keyword] = value

    sanitized = sanitize_fits_meta(table)

    assert sanitized.meta[keyword] == value
    sanitized.write(tmp_path / "lowband.fits", format="fits")


def test_sanitize_fits_meta_moves_oversized_inimage_to_history(tmp_path) -> None:
    table = Table({"source_id": [1]})
    keyword = "HIERARCH T1_INIMAGE"
    value = "1785941476_continuum_image_J2147-8132_IClean_highband.fits"
    table.meta[keyword] = value

    sanitized = sanitize_fits_meta(table)

    assert keyword not in sanitized.meta
    assert f"T1_INIMAGE={value}" in _history_text(sanitized)
    sanitized.write(tmp_path / "highband.fits", format="fits")


def test_sanitize_fits_meta_preserves_existing_history() -> None:
    table = Table({"source_id": [1]})
    table.meta["HISTORY"] = ["existing provenance"]
    table.meta["HIERARCH T1_INIMAGE"] = (
        "1785941476_continuum_image_J2147-8132_IClean_highband.fits"
    )

    sanitized = sanitize_fits_meta(table)

    history = sanitized.meta["HISTORY"]
    assert history[0] == "existing provenance"
    assert "T1_INIMAGE=" in "".join(history[1:])
