from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

from meerkat_corr_imaging import cli
from meerkat_corr_imaging.config import (
    Config,
    LowHighSliceCfg,
    PathsCfg,
    Target,
)
from meerkat_corr_imaging.image_pipeline import (
    build_image_pipeline_plan,
    wire_image_pipeline,
)


def _config(tmp_path: Path) -> Config:
    return Config(
        project_name="image-pipeline-test",
        paths=PathsCfg(
            reports_dir=str(tmp_path / "reports"),
            sky_xmatches_dir=str(tmp_path / "Sky-CrossMatches"),
        ),
        reference=Target(
            name="reference",
            images=[str(tmp_path / "reference_mfs.fits")],
            cuboid=str(tmp_path / "reference_cuboid.fits"),
        ),
        tests=[
            Target(
                name="test_32k",
                images=[str(tmp_path / "test_mfs.fits")],
                cuboid=str(tmp_path / "test_cuboid.fits"),
            )
        ],
        low_high_slice=LowHighSliceCfg(
            enabled=True,
            add_to_source_finding=True,
        ),
    )


def test_image_plan_chains_catalogues_xmatches_positions_and_flux(tmp_path: Path) -> None:
    cfg = _config(tmp_path)

    plan = build_image_pipeline_plan(cfg)

    assert {product.band for product in plan.products} == {"mfs", "low", "high"}
    assert len(plan.products) == 6
    assert len(plan.xmatch_pairs) == 3
    assert {match.band for match in plan.matches} == {"mfs", "low", "high"}
    assert len(plan.positions) == 3
    assert len(plan.flux) == 1

    catalogues = {product.catalogue for product in plan.products}
    for input1, input2, output in plan.xmatch_pairs:
        assert input1 in catalogues
        assert input2 in catalogues
        assert Path(output).parent == tmp_path / "Sky-CrossMatches"

    flux = plan.flux[0]
    tables_by_band = {match.band: match.table for match in plan.matches}
    assert flux == {
        "ref_low_xmatch": tables_by_band["low"],
        "ref_high_xmatch": tables_by_band["high"],
        "ref_mfs_xmatch": tables_by_band["mfs"],
    }


def test_image_plan_keeps_explicit_handoffs_and_fills_missing_bands(tmp_path: Path) -> None:
    cfg = _config(tmp_path)
    initial = build_image_pipeline_plan(cfg)
    mfs = next(match for match in initial.matches if match.band == "mfs")
    custom_mfs = str(tmp_path / "Sky-CrossMatches" / "custom_mfs.fits")
    cfg.extra = {
        "xmatch_pairs": [[mfs.reference.catalogue, mfs.test.catalogue, custom_mfs]],
        "positions": [
            {
                "xmatch_table": custom_mfs,
                "ref_fits": mfs.reference.image,
                "other_fits": mfs.test.image,
            }
        ],
        "flux": {},
    }

    plan = wire_image_pipeline(cfg)

    assert len(plan.xmatch_pairs) == 3
    assert sum(pair[2] == custom_mfs for pair in plan.xmatch_pairs) == 1
    assert len(plan.positions) == 3
    assert sum(entry["xmatch_table"] == custom_mfs for entry in plan.positions) == 1
    assert plan.flux[0]["ref_mfs_xmatch"] == custom_mfs
    assert cfg.extra["xmatch_pairs"] == [list(pair) for pair in plan.xmatch_pairs]


def test_images_command_runs_only_image_steps_in_order(monkeypatch, tmp_path: Path) -> None:
    cfg = _config(tmp_path)
    events: list[str] = []
    runners = {
        command: (command, lambda config, command=command: events.append(command))
        for command in cli.IMAGE_STEPS
    }

    monkeypatch.setattr(cli, "load_config", lambda path: cfg)
    monkeypatch.setattr(cli, "STEP_RUNNERS", runners)
    monkeypatch.setattr(
        cli,
        "wire_image_pipeline",
        lambda config: SimpleNamespace(
            products=(), xmatch_pairs=(), positions=(), flux=()
        ),
    )
    monkeypatch.setattr(
        cli,
        "run_step_with_audit",
        lambda step_name, reports_dir, runner: runner(),
    )

    cli.main(["--config", "unused.yaml", "images"])

    assert events == ["low_high_slice", "src", "xm", "pos", "flux", "report"]
