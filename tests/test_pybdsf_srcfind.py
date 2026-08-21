from __future__ import annotations

from pathlib import Path

from meerkat_corr_imaging import pybdsf_srcfind
from meerkat_corr_imaging.config import Config, PyBDSFCfg, Target
from meerkat_corr_imaging.steps import step5_srcfind


def _write_completed_catalogues(image: Path) -> None:
    base = pybdsf_srcfind.compute_base_name(image.stem, None)
    catalogues = pybdsf_srcfind.catalogue_paths(str(image.resolve()), base)
    catalogues[0].parent.mkdir(parents=True)
    for catalogue in catalogues:
        catalogue.touch()


def test_main_skips_only_inputs_with_completed_catalogues(
    tmp_path: Path,
    monkeypatch,
    capsys,
) -> None:
    completed = tmp_path / "completed.fits"
    pending = tmp_path / "pending.fits"
    completed.touch()
    pending.touch()
    _write_completed_catalogues(completed)

    processed: list[str] = []

    def fake_find_sources(image_path: str, **kwargs):
        processed.append(image_path)
        return str(tmp_path / "outputs"), str(tmp_path / "pybdsf.log")

    monkeypatch.setattr(pybdsf_srcfind, "find_sources", fake_find_sources)

    assert pybdsf_srcfind.main(["--images", str(completed), str(pending)]) == 0
    assert processed == [str(pending.resolve())]
    assert f"Skipping {completed.resolve()}" in capsys.readouterr().out


def test_main_overwrite_reprocesses_completed_input(
    tmp_path: Path,
    monkeypatch,
) -> None:
    image = tmp_path / "completed.fits"
    image.touch()
    _write_completed_catalogues(image)

    processed: list[str] = []

    def fake_find_sources(image_path: str, **kwargs):
        processed.append(image_path)
        return str(tmp_path / "outputs"), str(tmp_path / "pybdsf.log")

    monkeypatch.setattr(pybdsf_srcfind, "find_sources", fake_find_sources)

    assert pybdsf_srcfind.main(["--images", str(image), "--overwrite"]) == 0
    assert processed == [str(image.resolve())]


def test_step_passes_configured_overwrite_flag(monkeypatch) -> None:
    cfg = Config(
        project_name="overwrite-test",
        reference=Target(name="reference", images=["image.fits"]),
        pybdsf=PyBDSFCfg(overwrite=True),
    )
    commands: list[list[str]] = []
    monkeypatch.setattr(
        step5_srcfind,
        "_run",
        lambda command, **kwargs: commands.append(list(command)),
    )

    step5_srcfind.run(cfg)

    assert commands[0][-1] == "--overwrite"


def test_step_can_disable_adaptive_rms_boxes(monkeypatch) -> None:
    cfg = Config(
        project_name="fixed-rms-box-test",
        reference=Target(name="reference", images=["image.fits"]),
        pybdsf=PyBDSFCfg(adaptive_rms_box=False),
    )
    commands: list[list[str]] = []
    monkeypatch.setattr(
        step5_srcfind,
        "_run",
        lambda command, **kwargs: commands.append(list(command)),
    )

    step5_srcfind.run(cfg)

    assert "--no-adaptive-rms-box" in commands[0]


def test_main_passes_disabled_adaptive_rms_setting(
    tmp_path: Path,
    monkeypatch,
) -> None:
    image = tmp_path / "image.fits"
    image.touch()
    captured: dict[str, object] = {}

    def fake_find_sources(image_path: str, **kwargs):
        captured.update(kwargs)
        return str(tmp_path / "outputs"), str(tmp_path / "pybdsf.log")

    monkeypatch.setattr(pybdsf_srcfind, "find_sources", fake_find_sources)

    assert (
        pybdsf_srcfind.main(
            ["--images", str(image), "--no-adaptive-rms-box"]
        )
        == 0
    )
    assert captured["adaptive_rms_box"] is False
