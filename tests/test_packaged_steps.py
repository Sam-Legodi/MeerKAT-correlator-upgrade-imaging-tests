from __future__ import annotations

import sys
from pathlib import Path

from meerkat_corr_imaging.config import Config, PathsCfg, Target
from meerkat_corr_imaging.steps import (
    step2_vis_analysis,
    step3_calibrate_image,
    step5_srcfind,
    step6_xmatch,
    step7_flux,
    step7_positions,
)


def _config(tmp_path: Path) -> Config:
    cfg = Config(
        project_name="packaged-step-test",
        paths=PathsCfg(
            raw_dir=str(tmp_path / "raw"),
            interim_dir=str(tmp_path / "interim"),
            processed_dir=str(tmp_path / "processed"),
            reports_dir=str(tmp_path / "reports"),
            sky_xmatches_dir=str(tmp_path / "xmatches"),
        ),
        reference=Target(
            name="reference",
            ms_paths=["reference.ms"],
            images=["reference.fits"],
        ),
        tests=[
            Target(
                name="test",
                ms_paths=["test.ms"],
                images=["test.fits"],
            )
        ],
    )
    cfg.extra = {
        "xmatch_pairs": [["reference-cat.fits", "test-cat.fits", "matched.fits"]],
        "positions": [
            {
                "xmatch_table": "matched.fits",
                "ref_fits": "reference.fits",
                "other_fits": "test.fits",
            }
        ],
        "flux": {
            "ref_low_xmatch": "low.fits",
            "ref_high_xmatch": "high.fits",
            "ref_mfs_xmatch": "mfs.fits",
        },
    }
    return cfg


def test_python_steps_launch_packaged_modules_from_any_directory(
    tmp_path: Path,
    monkeypatch,
) -> None:
    cfg = _config(tmp_path)
    monkeypatch.chdir(tmp_path)

    commands: list[list[str]] = []
    monkeypatch.setattr(
        step2_vis_analysis,
        "_run_cmd",
        lambda cmd, **kwargs: commands.append(list(cmd)),
    )
    monkeypatch.setattr(
        step5_srcfind,
        "_run",
        lambda cmd, **kwargs: commands.append(list(cmd)),
    )
    monkeypatch.setattr(
        step6_xmatch,
        "_run",
        lambda cmd, **kwargs: commands.append(list(cmd)),
    )
    monkeypatch.setattr(
        step7_positions,
        "_run",
        lambda cmd, **kwargs: commands.append(list(cmd)),
    )
    monkeypatch.setattr(
        step7_flux,
        "_run",
        lambda cmd, **kwargs: commands.append(list(cmd)),
    )

    step2_vis_analysis.run(cfg)
    step5_srcfind.run(cfg)
    step6_xmatch.run(cfg)
    step7_positions.run(cfg)
    step7_flux.run(cfg)

    modules = [
        command[2]
        for command in commands
        if command[:2] == [sys.executable, "-m"]
    ]
    assert modules == [
        "meerkat_corr_imaging.vis_amp_analyze",
        "meerkat_corr_imaging.vis_amp_analyze",
        "meerkat_corr_imaging.pybdsf_srcfind",
        "meerkat_corr_imaging.xmatch_pybdsf",
        "meerkat_corr_imaging.positions_analysis",
        "meerkat_corr_imaging.flux_analysis",
    ]
    visibility_commands = [
        command
        for command in commands
        if command[:3]
        == [sys.executable, "-m", "meerkat_corr_imaging.vis_amp_analyze"]
    ]
    assert all("--exact-outdir" in command for command in visibility_commands)
    assert "--ms-ref-results" not in visibility_commands[0]
    assert "--ms-ref-results" in visibility_commands[1]
    assert "--ms-ref" not in visibility_commands[1]


def test_casa_step_uses_absolute_packaged_scripts(tmp_path: Path, monkeypatch) -> None:
    cfg = _config(tmp_path)
    cfg.extra["force_calibrate"] = True
    monkeypatch.chdir(tmp_path)

    for name in ("reference.ms", "test.ms"):
        (tmp_path / name).mkdir()
    monkeypatch.setenv("MCI_TCLEAN_MSFILE", "stale.ms")
    calls = []
    monkeypatch.setattr(
        step3_calibrate_image, "_run_cmd",
        lambda cmd, env=None, **kwargs: calls.append((list(cmd), dict(env))),
    )
    step3_calibrate_image.run(cfg)

    assert len(calls) == 4
    for index, (command, env) in enumerate(calls):
        script = Path(command[command.index("-c") + 1])
        assert script.is_absolute() and script.is_file()
        assert script.parent.name == "meerkat_corr_imaging"
        assert env["MCI_CASA_SCRIPT_DIR"] == str(script.parent)
        ms = str(tmp_path / ("reference.ms" if index < 2 else "test.ms"))
        if index % 2 == 0:
            assert script.name == "standalone_xxyy_solve.py"
            assert env["MCI_CAL_MSFILE"] == ms
            assert env["MCI_CAL_REFANT"] == cfg.casa.refant
            assert env["MCI_CAL_FLUX_FIELD"] == cfg.casa.flux_field
        else:
            assert script.name == "tclean_two_bands.py"
            assert command[-1] == env["MCI_TCLEAN_MSFILE"] == ms


def test_casa_missing_input_fails_before_launch(tmp_path, monkeypatch):
    import pytest
    cfg = _config(tmp_path)
    cfg.extra["force_calibrate"] = True
    monkeypatch.chdir(tmp_path)
    calls = []
    monkeypatch.setattr(step3_calibrate_image, "_run_cmd", lambda *a, **k: calls.append(a))
    with pytest.raises(FileNotFoundError, match="reference.ms.*test.ms"):
        step3_calibrate_image.run(cfg)
    assert calls == []


def test_solver_requires_explicit_existing_input(monkeypatch, tmp_path):
    import pytest
    # Execute the input validation in isolation without importing CASA tasks.
    script = Path(step3_calibrate_image.__file__).resolve().parents[1] / "standalone_xxyy_solve.py"
    source = script.read_text().split("# MeasurementSet to calibrate\n", 1)[1].split('antennas =', 1)[0]
    import os
    monkeypatch.delenv("MCI_CAL_MSFILE", raising=False)
    with pytest.raises(ValueError, match="MCI_CAL_MSFILE"):
        exec(source, {"os": os})
    monkeypatch.setenv("MCI_CAL_MSFILE", str(tmp_path / "missing.ms"))
    with pytest.raises(IOError, match="missing.ms"):
        exec(source, {"os": os})


def test_flux_step_runs_each_chained_analysis(tmp_path: Path, monkeypatch) -> None:
    cfg = _config(tmp_path)
    cfg.extra["flux"] = [
        {
            "ref_low_xmatch": "test-a-low.fits",
            "ref_high_xmatch": "test-a-high.fits",
            "ref_mfs_xmatch": "test-a-mfs.fits",
        },
        {
            "ref_low_xmatch": "test-b-low.fits",
            "ref_high_xmatch": "test-b-high.fits",
            "ref_mfs_xmatch": "test-b-mfs.fits",
        },
    ]
    commands: list[list[str]] = []
    monkeypatch.setattr(
        step7_flux,
        "_run",
        lambda cmd, **kwargs: commands.append(list(cmd)),
    )

    step7_flux.run(cfg)

    assert len(commands) == 2
    assert commands[0][2] == "meerkat_corr_imaging.flux_analysis"
    assert "test-a-mfs.fits" in commands[0]
    assert "test-b-mfs.fits" in commands[1]


def test_calibration_only_passes_exclusions(tmp_path, monkeypatch):
    import json
    cfg = _config(tmp_path)
    cfg.extra['force_calibrate'] = True
    cfg.casa.imaging_enabled = False
    cfg.casa.exclude_fields = ['target']
    cfg.casa.calibration_exclude_fields = ['2']
    monkeypatch.chdir(tmp_path)
    for ms in ('reference.ms', 'test.ms'):
        (tmp_path / ms).mkdir()
    calls = []
    monkeypatch.setattr(step3_calibrate_image, '_run_cmd',
                        lambda cmd, env=None, **kw: calls.append((cmd, env)))
    step3_calibrate_image.run(cfg)
    assert len(calls) == 2
    for cmd, env in calls:
        assert cmd[-1].endswith('standalone_xxyy_solve.py')
        assert json.loads(env['MCI_CAL_EXCLUDE_FIELDS']) == ['target', '2']
