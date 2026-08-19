"""Extract low- and high-band image planes from an MFImage FITS cuboid.

Each output is one existing independent subband plane, not an average or a
re-image.  Among planes that overlap the requested frequency range, selection
maximises frequency overlap.  Equal-overlap ties choose the lowest effective
frequency for the low-band output and the highest effective frequency for the
high-band output.
"""

from __future__ import annotations

import argparse
import math
import os
import re
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
from astropy.io import fits


FITS_SUFFIXES = {".fits", ".fit", ".fts"}


@dataclass(frozen=True)
class Subband:
    """Metadata for one independent MFImage subband plane."""

    number: int
    plane_index: int
    low_hz: float
    effective_hz: float
    high_hz: float

    @property
    def width_hz(self) -> float:
        return self.high_hz - self.low_hz

    def overlap_hz(self, band_hz: Sequence[float]) -> float:
        low_hz, high_hz = validate_band_range(band_hz, "requested")
        return max(0.0, min(self.high_hz, high_hz) - max(self.low_hz, low_hz))


@dataclass(frozen=True)
class Selection:
    """A selected cuboid plane and its overlap with the requested range."""

    band_name: str
    requested_low_hz: float
    requested_high_hz: float
    subband: Subband
    overlap_hz: float


def validate_band_range(values: Sequence[float], label: str) -> tuple[float, float]:
    if len(values) != 2:
        raise ValueError(f"{label} range must contain exactly two frequencies")
    low_hz, high_hz = (float(values[0]), float(values[1]))
    if not (math.isfinite(low_hz) and math.isfinite(high_hz)):
        raise ValueError(f"{label} range must contain finite frequencies")
    if low_hz >= high_hz:
        raise ValueError(f"{label} range must be increasing: {low_hz}, {high_hz}")
    return low_hz, high_hz


def subbands_from_header(header: fits.Header) -> list[Subband]:
    """Read the notebook-verified MFImage plane and frequency mapping."""
    if str(header.get("CTYPE3", "")).strip() != "SPECLNMF":
        raise ValueError("Input is not an MFImage cuboid: CTYPE3 must be SPECLNMF")

    try:
        nterm = int(header["NTERM"])
        nspec = int(header["NSPEC"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("MFImage cuboid must contain integer NTERM and NSPEC") from exc
    if nterm < 1 or nspec < 1:
        raise ValueError(f"Invalid MFImage dimensions: NTERM={nterm}, NSPEC={nspec}")

    subbands: list[Subband] = []
    for number in range(1, nspec + 1):
        keys = {
            "low": f"FREL{number:04d}",
            "effective": f"FEFF{number:04d}",
            "high": f"FREH{number:04d}",
        }
        missing = [key for key in keys.values() if key not in header]
        if missing:
            raise ValueError(
                "MFImage cuboid lacks subband frequency metadata: " + ", ".join(missing)
            )
        low_hz = float(header[keys["low"]])
        effective_hz = float(header[keys["effective"]])
        high_hz = float(header[keys["high"]])
        if not all(math.isfinite(value) for value in (low_hz, effective_hz, high_hz)):
            raise ValueError(f"Subband {number} contains non-finite frequency metadata")
        if not low_hz < high_hz:
            raise ValueError(f"Subband {number} has invalid edges: {low_hz}, {high_hz}")
        if not low_hz <= effective_hz <= high_hz:
            raise ValueError(
                f"Subband {number} effective frequency {effective_hz} is outside its edges"
            )
        subbands.append(
            Subband(
                number=number,
                plane_index=nterm + number - 1,
                low_hz=low_hz,
                effective_hz=effective_hz,
                high_hz=high_hz,
            )
        )
    return subbands


def select_subband(
    subbands: Iterable[Subband],
    band_hz: Sequence[float],
    band_name: str,
) -> Selection:
    """Select the maximum-overlap plane, with a frequency-directed tie-break."""
    requested_low_hz, requested_high_hz = validate_band_range(band_hz, band_name)
    candidates: list[tuple[float, Subband]] = []
    for subband in subbands:
        overlap_hz = subband.overlap_hz((requested_low_hz, requested_high_hz))
        if overlap_hz > 0.0:
            candidates.append((overlap_hz, subband))
    if not candidates:
        raise ValueError(
            f"No cuboid subband overlaps the {band_name} range "
            f"{requested_low_hz:.6g}-{requested_high_hz:.6g} Hz"
        )

    if band_name == "lowband":
        overlap_hz, selected = max(
            candidates,
            key=lambda item: (item[0], -item[1].effective_hz),
        )
    elif band_name == "highband":
        overlap_hz, selected = max(
            candidates,
            key=lambda item: (item[0], item[1].effective_hz),
        )
    else:
        raise ValueError(f"Unsupported band name: {band_name}")

    return Selection(
        band_name=band_name,
        requested_low_hz=requested_low_hz,
        requested_high_hz=requested_high_hz,
        subband=selected,
        overlap_hz=overlap_hz,
    )


def default_output_path(cuboid_path: Path, band_name: str) -> Path:
    return cuboid_path.with_name(f"{cuboid_path.stem}_{band_name}.fits")


def _remove_non_spatial_wcs(header: fits.Header) -> None:
    """Reduce the cuboid WCS to its two celestial axes."""
    scalar_axis = re.compile(r"^(?:CTYPE|CUNIT|CRVAL|CRPIX|CDELT|CROTA)(\d+)[A-Z]?$", re.I)
    matrix_axis = re.compile(r"^(?:PC|CD)(\d+)_(\d+)[A-Z]?$", re.I)
    parameter_axis = re.compile(r"^(?:PV|PS)(\d+)_(\d+)[A-Z]?$", re.I)
    for key in list(header.keys()):
        if not key:
            continue
        scalar_match = scalar_axis.match(key)
        matrix_match = matrix_axis.match(key)
        parameter_match = parameter_axis.match(key)
        if scalar_match and int(scalar_match.group(1)) > 2:
            del header[key]
        elif matrix_match and (
            int(matrix_match.group(1)) > 2 or int(matrix_match.group(2)) > 2
        ):
            del header[key]
        elif parameter_match and int(parameter_match.group(1)) > 2:
            del header[key]

    for key in (
        "NAXIS3", "NAXIS4", "NTERM", "NSPEC", "RFALPHA",
        "SPECSYS", "SSYSOBS", "VELOSYS",
    ):
        header.pop(key, None)
    for key in list(header.keys()):
        if re.match(r"^(?:FREQ|FREL|FEFF|FREH)\d{4}$", key, re.I):
            del header[key]
    header["WCSAXES"] = (2, "Number of coordinate axes")


def output_header(source: fits.Header, selection: Selection, cuboid_path: Path) -> fits.Header:
    """Build a 2-D, PyBDSF-readable FITS header for a selected plane."""
    header = source.copy()
    _remove_non_spatial_wcs(header)

    effective_hz = selection.subband.effective_hz
    header["RESTFRQ"] = (effective_hz, "Selected subband effective frequency (Hz)")
    header["FREQ"] = (effective_hz, "Selected subband effective frequency (Hz)")
    header["CFREQ"] = (effective_hz, "Selected subband effective frequency (Hz)")
    header["REQFLO"] = (selection.requested_low_hz, "Requested band lower edge (Hz)")
    header["REQFHI"] = (selection.requested_high_hz, "Requested band upper edge (Hz)")
    header["SUBFLO"] = (selection.subband.low_hz, "Selected subband lower edge (Hz)")
    header["SUBFEFF"] = (effective_hz, "Selected subband effective frequency (Hz)")
    header["SUBFHI"] = (selection.subband.high_hz, "Selected subband upper edge (Hz)")
    header["OVLPHZ"] = (selection.overlap_hz, "Overlap with requested band (Hz)")
    header["SUBBAND"] = (selection.subband.number, "MFImage subband number (1-based)")
    header["CUBPLANE"] = (selection.subband.plane_index, "Cuboid plane index (0-based)")
    header["MCI_BAND"] = (selection.band_name, "MCI low/high slice role")
    header["PBCOR"] = (False, "Image is not primary-beam corrected")
    header["BTYPE"] = (header.get("BTYPE", "Intensity"), "Image type")

    if "BMAJ" not in header and "CLEANBMJ" in source:
        header["BMAJ"] = (float(source["CLEANBMJ"]), "Restoring beam major axis (deg)")
    if "BMIN" not in header and "CLEANBMN" in source:
        header["BMIN"] = (float(source["CLEANBMN"]), "Restoring beam minor axis (deg)")
    if "BPA" not in header and "CLEANBPA" in source:
        header["BPA"] = (float(source["CLEANBPA"]), "Restoring beam position angle (deg)")

    header.add_history("MCI low_high_slice extracted one non-PB MFImage subband plane.")
    header.add_history(f"Input cuboid: {cuboid_path}")
    header.add_history(
        f"{selection.band_name}: selected subband {selection.subband.number}, "
        f"cuboid plane {selection.subband.plane_index}."
    )
    return header


def validate_existing_image(path: Path) -> None:
    if not path.is_file():
        raise FileNotFoundError(path)
    if path.suffix.lower() not in FITS_SUFFIXES:
        raise ValueError(f"Expected a FITS image path: {path}")
    with fits.open(path, memmap=True, do_not_scale_image_data=True) as hdul:
        if not hdul or hdul[0].data is None:
            raise ValueError(f"FITS image has no primary image data: {path}")
        if np.asarray(hdul[0].data).squeeze().ndim != 2:
            raise ValueError(f"Expected a two-dimensional image product: {path}")


def _write_atomic(output_path: Path, data: np.ndarray, header: fits.Header) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{output_path.stem}.", suffix=".tmp.fits", dir=output_path.parent
    )
    os.close(descriptor)
    temporary_path = Path(temporary_name)
    try:
        fits.PrimaryHDU(data=np.asarray(data, dtype=np.float32), header=header).writeto(
            temporary_path,
            overwrite=True,
            checksum=True,
        )
        os.replace(temporary_path, output_path)
    finally:
        if temporary_path.exists():
            temporary_path.unlink()


def extract_slices(
    cuboid_path: Path,
    lowband_hz: Sequence[float],
    highband_hz: Sequence[float],
    low_output: Path | None = None,
    high_output: Path | None = None,
    overwrite: bool = False,
) -> tuple[Path, Path]:
    """Extract selected low/high cuboid planes and return their output paths."""
    cuboid_path = cuboid_path.expanduser().resolve()
    low_output = (low_output or default_output_path(cuboid_path, "lowband")).expanduser().resolve()
    high_output = (high_output or default_output_path(cuboid_path, "highband")).expanduser().resolve()
    if low_output == high_output:
        raise ValueError("Low- and high-band output paths must differ")
    if cuboid_path in (low_output, high_output):
        raise ValueError("An output path must not overwrite the input cuboid")

    outputs = {"lowband": low_output, "highband": high_output}
    pending = {
        name: path
        for name, path in outputs.items()
        if overwrite or not path.exists()
    }
    for name, path in pending.items():
        if path.parent != cuboid_path.parent:
            raise ValueError(
                f"A newly created {name} image must be beside its cuboid: "
                f"{path.parent} != {cuboid_path.parent}"
            )
    for name, path in outputs.items():
        if name not in pending:
            validate_existing_image(path)
            print(f"[{name}] Reusing existing FITS image: {path}")
    if not pending:
        return low_output, high_output

    if not cuboid_path.is_file():
        raise FileNotFoundError(cuboid_path)
    if cuboid_path.suffix.lower() not in FITS_SUFFIXES:
        raise ValueError(f"Expected a FITS cuboid: {cuboid_path}")

    with fits.open(cuboid_path, memmap=True, do_not_scale_image_data=True) as hdul:
        if not hdul or hdul[0].data is None:
            raise ValueError(f"FITS cuboid has no primary image data: {cuboid_path}")
        header = hdul[0].header
        raw_cube = hdul[0].data
        subbands = subbands_from_header(header)
        nterm = int(header["NTERM"])
        nspec = int(header["NSPEC"])

        if raw_cube.ndim == 4:
            if raw_cube.shape[0] != 1:
                raise ValueError(f"Expected one Stokes plane, received shape {raw_cube.shape}")
            image_planes = raw_cube[0]
        elif raw_cube.ndim == 3:
            image_planes = raw_cube
        else:
            raise ValueError(f"Expected a 3-D/4-D MFImage cuboid, received {raw_cube.shape}")
        if image_planes.shape[0] != nterm + nspec:
            raise ValueError(
                f"Cuboid has {image_planes.shape[0]} image planes; "
                f"NTERM + NSPEC requires {nterm + nspec}"
            )

        selections = {
            "lowband": select_subband(subbands, lowband_hz, "lowband"),
            "highband": select_subband(subbands, highband_hz, "highband"),
        }
        for name, output_path in pending.items():
            selection = selections[name]
            subband = selection.subband
            print(
                f"[{name}] plane={subband.plane_index} subband={subband.number} "
                f"edges={subband.low_hz:.6g}-{subband.high_hz:.6g} Hz "
                f"effective={subband.effective_hz:.6g} Hz "
                f"overlap={selection.overlap_hz:.6g} Hz"
            )
            header_2d = output_header(header, selection, cuboid_path)
            _write_atomic(output_path, image_planes[subband.plane_index], header_2d)
            validate_existing_image(output_path)
            print(f"[{name}] Wrote: {output_path}")

    return low_output, high_output


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Extract maximum-coverage low/high subband planes from an MFImage cuboid."
    )
    parser.add_argument("--cuboid", required=True, type=Path, help="Input non-PB MFImage FITS cuboid")
    parser.add_argument(
        "--lowband-hz", required=True, type=float, nargs=2, metavar=("LOW_HZ", "HIGH_HZ")
    )
    parser.add_argument(
        "--highband-hz", required=True, type=float, nargs=2, metavar=("LOW_HZ", "HIGH_HZ")
    )
    parser.add_argument("--low-output", type=Path, default=None)
    parser.add_argument("--high-output", type=Path, default=None)
    parser.add_argument("--overwrite", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    extract_slices(
        cuboid_path=args.cuboid,
        lowband_hz=args.lowband_hz,
        highband_hz=args.highband_hz,
        low_output=args.low_output,
        high_output=args.high_output,
        overwrite=args.overwrite,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
