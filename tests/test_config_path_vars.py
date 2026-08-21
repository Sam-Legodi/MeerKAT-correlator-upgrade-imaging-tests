from __future__ import annotations

from pathlib import Path

import pytest

from meerkat_corr_imaging.config import _dict_to_dataclass, load_config


def test_path_vars_expand_recursively_in_nested_config_values(tmp_path: Path) -> None:
    config_path = tmp_path / "config.yaml"
    raw_config = """
project_name: path-vars-test
path_vars:
  root: __ROOT__
  observations: ${root}/Obsevations
  image_dir: ${observations}/123/images
  stem: target_image
paths:
  processed_dir: ${root}/processed
reference:
  images: ["${image_dir}/${stem}_PB.fits"]
  cuboid: "${image_dir}/${stem}.fits"
extra:
  xmatch_pairs:
    - ["${image_dir}/ref.fits", "${image_dir}/test.fits"]
""".strip()
    root = tmp_path / "commissioning"
    config_path.write_text(raw_config.replace("__ROOT__", str(root)))

    cfg = load_config(config_path)

    assert cfg.paths.processed_dir == str(root / "processed")
    assert cfg.reference.images == [
        str(root / "Obsevations/123/images/target_image_PB.fits")
    ]
    assert cfg.reference.cuboid == str(
        root / "Obsevations/123/images/target_image.fits"
    )
    assert cfg.extra["xmatch_pairs"][0][0] == str(
        root / "Obsevations/123/images/ref.fits"
    )
    assert "path_vars" not in cfg.extra


def test_path_vars_reject_cycles() -> None:
    with pytest.raises(ValueError, match="Cyclic path_vars reference: first -> second -> first"):
        _dict_to_dataclass(
            {
                "project_name": "cycle-test",
                "path_vars": {
                    "first": "${second}",
                    "second": "${first}",
                },
            }
        )


def test_path_vars_reject_undefined_references() -> None:
    with pytest.raises(ValueError, match="Undefined path variable: missing"):
        _dict_to_dataclass(
            {
                "project_name": "undefined-test",
                "path_vars": {"root": "/data"},
                "reference": {"images": ["${missing}/image.fits"]},
            }
        )
