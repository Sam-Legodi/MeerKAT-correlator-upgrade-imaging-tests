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


def test_casa_step_uses_absolute_packaged_scripts(tmp_path: Path, monkeypatch) -> None:
    cfg = _config(tmp_path)
    cfg.extra["force_calibrate"] = True
    monkeypatch.chdir(tmp_path)

    commands: list[list[str]] = []
    monkeypatch.setattr(
        step3_calibrate_image,
        "_run_cmd",
        lambda cmd, env=None, **kwargs: commands.append(list(cmd)),
    )
    step3_calibrate_image.run(cfg)

    assert len(commands) == 3
    for command in commands:
        script = Path(command[command.index("-c") + 1])
        assert script.is_absolute()
        assert script.is_file()
        assert script.parent.name == "meerkat_corr_imaging"
    assert commands[0][commands[0].index("-c") + 1].endswith("standalone_xxyy_solve.py")
    assert commands[1][commands[1].index("-c") + 1].endswith("tclean_two_bands.py")
    assert commands[2][commands[2].index("-c") + 1].endswith("tclean_two_bands.py")
    assert commands[1][-1] == "reference.ms"
    assert commands[2][-1] == "test.ms"
