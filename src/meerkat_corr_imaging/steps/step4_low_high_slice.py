from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Optional

from astropy.io import fits

from ..audit import (
    raise_for_failures,
    record_failure,
    record_skip,
    record_success,
    register_inputs,
    run_logged_command,
)
from ..config import Config, Target


@dataclass(frozen=True)
class ResolvedSliceTarget:
    role: str
    name: str
    cuboid: Optional[Path]
    low_image: Path
    high_image: Path


def _absolute(path: str) -> Path:
    return Path(path).expanduser().resolve()


def _default_output(cuboid: Path, band_name: str) -> Path:
    return cuboid.with_name(f"{cuboid.stem}_{band_name}.fits")


def _resolve_target(role: str, spec: Target) -> ResolvedSliceTarget:
    cuboid = _absolute(spec.cuboid) if spec.cuboid else None
    if spec.low_image:
        low_image = _absolute(spec.low_image)
    elif cuboid is not None:
        low_image = _default_output(cuboid, "lowband")
    else:
        raise ValueError(f"low_high_slice {role} '{spec.name}' needs cuboid or low_image")

    if spec.high_image:
        high_image = _absolute(spec.high_image)
    elif cuboid is not None:
        high_image = _default_output(cuboid, "highband")
    else:
        raise ValueError(f"low_high_slice {role} '{spec.name}' needs cuboid or high_image")

    return ResolvedSliceTarget(
        role=role,
        name=spec.name or role,
        cuboid=cuboid,
        low_image=low_image,
        high_image=high_image,
    )


def resolved_targets(cfg: Config) -> List[ResolvedSliceTarget]:
    """Resolve configured or deterministic low/high paths without creating files."""
    section = cfg.low_high_slice
    if not section.enabled:
        return []

    targets: List[ResolvedSliceTarget] = []
    reference = cfg.reference
    if reference.cuboid or reference.low_image or reference.high_image:
        targets.append(_resolve_target("reference", reference))

    seen_names = set()
    for spec in cfg.tests:
        if not spec.name:
            raise ValueError("Each tests entry requires a name")
        if spec.name in seen_names:
            raise ValueError(f"Duplicate test name: {spec.name}")
        seen_names.add(spec.name)
        if spec.cuboid or spec.low_image or spec.high_image:
            targets.append(_resolve_target("test", spec))
    return targets


def _validate_image(path: Path, label: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(
            f"Missing {label}: {path}. Run the low_high_slice step before this stage."
        )
    with fits.open(path, memmap=True, do_not_scale_image_data=True) as hdul:
        if not hdul or hdul[0].data is None:
            raise ValueError(f"Configured {label} has no primary image data: {path}")
        if hdul[0].data.squeeze().ndim != 2:
            raise ValueError(f"Configured {label} is not a two-dimensional image: {path}")


def source_images(cfg: Config) -> List[str]:
    """Return validated low/high images for automatic PyBDSF inclusion."""
    section = cfg.low_high_slice
    if not section.enabled or not section.add_to_source_finding:
        return []
    images: List[str] = []
    for target in resolved_targets(cfg):
        _validate_image(target.low_image, f"{target.name} low-band image")
        _validate_image(target.high_image, f"{target.name} high-band image")
        images.extend((str(target.low_image), str(target.high_image)))
    return images


def _run(
    command: Iterable[str],
    *,
    inputs: Iterable[str] = (),
    mark_success: bool = True,
) -> None:
    run_logged_command(
        command,
        prefix="[LOW/HIGH SLICE]",
        inputs=inputs,
        mark_success=mark_success,
    )


def _audit_label(target: ResolvedSliceTarget) -> str:
    if target.cuboid is not None:
        return str(target.cuboid)
    return f"{target.low_image} + {target.high_image}"


def _process_target(
    cfg: Config,
    target: ResolvedSliceTarget,
    lowband_hz: Iterable[float],
    highband_hz: Iterable[float],
) -> None:
    label = _audit_label(target)
    low_exists = target.low_image.is_file()
    high_exists = target.high_image.is_file()
    if low_exists and high_exists and (not cfg.low_high_slice.overwrite or target.cuboid is None):
        _validate_image(target.low_image, f"{target.name} low-band image")
        _validate_image(target.high_image, f"{target.name} high-band image")
        print(f"[LOW/HIGH SLICE] Reusing existing images for {target.name}.")
        record_success(label, "reused existing low/high images")
        return
    if target.cuboid is None:
        missing = [
            str(path)
            for path in (target.low_image, target.high_image)
            if not path.is_file()
        ]
        raise FileNotFoundError(
            f"Cannot create missing image(s) for {target.name} without a cuboid: {missing}"
        )

    command = [
        sys.executable,
        "-m",
        "meerkat_corr_imaging.low_high_slice",
        "--cuboid",
        str(target.cuboid),
        "--lowband-hz",
        *(str(value) for value in lowband_hz),
        "--highband-hz",
        *(str(value) for value in highband_hz),
        "--low-output",
        str(target.low_image),
        "--high-output",
        str(target.high_image),
    ]
    if cfg.low_high_slice.overwrite:
        command.append("--overwrite")
    _run(command, inputs=[label], mark_success=False)
    _validate_image(target.low_image, f"{target.name} low-band image")
    _validate_image(target.high_image, f"{target.name} high-band image")
    record_success(label)


def run(cfg: Config) -> None:
    section = cfg.low_high_slice
    if not section.enabled:
        message = "Disabled in config; skipping."
        print(f"[LOW/HIGH SLICE] {message}")
        record_skip(message)
        return

    targets = resolved_targets(cfg)
    if not targets:
        message = "No reference/test cuboids or images configured; skipping."
        print(f"[LOW/HIGH SLICE] {message}")
        record_skip(message)
        return

    register_inputs(_audit_label(target) for target in targets)

    lowband_hz = section.lowband_hz or cfg.frequency_ranges.lowband_hz
    highband_hz = section.highband_hz or cfg.frequency_ranges.highband_hz

    failures: list[BaseException] = []
    for target in targets:
        try:
            _process_target(cfg, target, lowband_hz, highband_hz)
        except Exception as exc:
            label = _audit_label(target)
            record_failure(label, exc)
            print(f"[LOW/HIGH SLICE] FAILED {label}: {exc}")
            failures.append(exc)

    raise_for_failures("low/high slicing", failures)
