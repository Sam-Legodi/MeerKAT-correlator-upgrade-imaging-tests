#!/usr/bin/env python3
"""Cross-match PyBDSF-style FITS catalogues on the celestial sphere.

The script reads two source catalogues, performs a spherical sky match that
mirrors TOPCAT's "sky" algorithm, and writes only the positive matches to a
third FITS file. Matching defaults to a greedy global nearest-neighbour
approach (minimum separation first) to enforce one-to-one associations, but a
one-to-many mode is also available. Input/output paths can be supplied either
positionally or via ``--paths-file`` containing a single comma-separated line
``input1,input2,output``.

Usage
-----

    python xmatch_pybdsf.py cat1.fits cat2.fits output.fits

See ``python xmatch_pybdsf.py --help`` for the full CLI documentation.
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
import tempfile
from dataclasses import dataclass
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np
from astropy.coordinates import Angle, SkyCoord, search_around_sky
from astropy.io import fits
from astropy.table import MaskedColumn, Table
import astropy.units as u


LOG = logging.getLogger("xmatch_pybdsf")


@dataclass
class PreparedCatalog:
    """Bundle table data and the subset that is valid for matching."""

    table: Table
    coords: SkyCoord
    valid_indices: np.ndarray
    header_meta: Dict[str, object]
    skipped: int


def parse_arguments(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(
        description="Cross-match two PyBDSF-style FITS catalogues on the sky.",
    )
    parser.add_argument("input1", nargs="?", help="First input FITS catalogue")
    parser.add_argument("input2", nargs="?", help="Second input FITS catalogue")
    parser.add_argument(
        "output",
        nargs="?",
        help="Output FITS filename; results are written into a sibling 'Sky-CrossMatches/' directory.",
    )
    parser.add_argument(
        "--paths-file",
        help="Optional text file containing one or more lines 'input1,input2,output'.",
    )

    parser.add_argument(
        "--max-error",
        default="1 arcsec",
        help=(
            "Maximum angular separation for a positive match (default: 1 arcsec). "
            "Units are parsed with astropy and can be arcsec, deg, etc."
        ),
    )
    parser.add_argument("--ra-col-1", default="RA", help="RA column name in catalogue 1")
    parser.add_argument("--dec-col-1", default="DEC", help="Dec column name in catalogue 1")
    parser.add_argument("--ra-col-2", default="RA", help="RA column name in catalogue 2")
    parser.add_argument("--dec-col-2", default="DEC", help="Dec column name in catalogue 2")
    parser.add_argument(
        "--coord-frame",
        default="icrs",
        help="Reference frame name understood by astropy (default: icrs)",
    )
    parser.add_argument(
        "--one-to-many",
        action="store_true",
        help="Report all matches within the radius (no one-to-one enforcement).",
    )
    parser.add_argument(
        "--progress",
        action="store_true",
        help="Print progress information (sets logging level to INFO).",
    )
    parser.add_argument(
        "--self-test",
        action="store_true",
        help="Run a built-in self-test that fabricates small catalogues.",
    )

    args = parser.parse_args(argv)

    jobs: List[Tuple[str, str, str]] = []

    if args.paths_file:
        try:
            with open(args.paths_file, "r", encoding="utf-8") as fh:
                for raw_line in fh:
                    line = raw_line.strip()
                    if not line or line.startswith("#"):
                        continue
                    parts = [p.strip() for p in line.split(",")]
                    if len(parts) < 3:
                        raise ValueError("expected three comma-separated paths")
                    jobs.append((parts[0], parts[1], parts[2]))
        except Exception as exc:
            parser.error(f"Failed to parse --paths-file '{args.paths_file}': {exc}")

    if args.input1 or args.input2 or args.output:
        if not (args.input1 and args.input2 and args.output):
            parser.error("When supplying positional arguments you must provide input1, input2, and output")
        jobs.insert(0, (args.input1, args.input2, args.output))

    if not args.self_test and not jobs:
        parser.error("No input/output paths supplied (use positional args, --paths-file, or --self-test)")

    if jobs:
        args.input1, args.input2, args.output = jobs[0]
    args.jobs = jobs

    return args


def configure_logging(progress: bool, log_path: str | None = None) -> None:
    """Configure logging to stdout and optional logfile with timestamps."""

    level = logging.INFO if progress else logging.WARNING

    root_logger = logging.getLogger()
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)

    formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S")

    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(level)
    stream_handler.setFormatter(formatter)
    root_logger.addHandler(stream_handler)

    if log_path:
        file_handler = logging.FileHandler(log_path, mode="a")
        file_handler.setLevel(level)
        file_handler.setFormatter(formatter)
        root_logger.addHandler(file_handler)

    root_logger.setLevel(level)


def parse_angle(text: str) -> Angle:
    """Parse an angular quantity using ``astropy.coordinates.Angle``."""

    try:
        ang = Angle(text)
    except Exception as exc:  # pragma: no cover - defensive fallback
        raise ValueError(f"Could not parse angle '{text}': {exc}") from exc
    if not ang.isscalar or not np.isfinite(ang.arcsec) or ang.arcsec <= 0:
        raise ValueError(f"Maximum error must be positive and finite: '{text}'")
    return ang


def _column_to_quantity(column, default_unit: u.Unit) -> u.Quantity:
    """Convert a table column to a quantity, assuming a default unit if absent."""

    data = np.asarray(column)
    unit = getattr(column, "unit", None)
    if unit is None:
        qty = u.Quantity(data, default_unit)
    else:
        qty = u.Quantity(data, unit)
    return qty


def _valid_coordinate_mask(
    ra_qty: u.Quantity, dec_qty: u.Quantity, ra_col, dec_col
) -> np.ndarray:
    """Return mask of finite (RA, Dec) rows, respecting masked values."""

    finite = np.isfinite(ra_qty.to_value(u.deg)) & np.isfinite(dec_qty.to_value(u.deg))
    if isinstance(ra_col, MaskedColumn):
        finite &= ~ra_col.mask
    if isinstance(dec_col, MaskedColumn):
        finite &= ~dec_col.mask
    ra = ra_qty.to_value(u.deg)
    dec = dec_qty.to_value(u.deg)
    return finite & (ra >= 0) & (ra <= 360) & (np.abs(dec) <= 90)


def _extract_header_metadata(path: str, prefix: str, hdu=None) -> Dict[str, object]:
    """Extract FITS header keywords and prefix them for inclusion in the output."""

    meta: Dict[str, object] = {}
    with fits.open(path) as hdul:
        header = hdul[hdu].header if hdu is not None else None
        if header is None:
            for item in hdul:
                if isinstance(item, (fits.BinTableHDU, fits.TableHDU)):
                    header = item.header
                    break
        if header is None:
            header = hdul[0].header
        for key, value in header.items():
            if not key or key in {"HISTORY", "COMMENT"}:
                continue
            meta_key = f"HIERARCH {prefix}_{key}"
            meta[meta_key] = value
    return meta


def prepare_catalog(
    path: str,
    ra_col_name: str,
    dec_col_name: str,
    frame: str,
    prefix: str,
    *,
    format: str | None = None,
    hdu: int | str | None = None,
    table_id: int | str | None = None,
    coordinate_unit: str = "deg",
    required_columns: Sequence[str] = (),
    strict_table_selection: bool = False,
) -> PreparedCatalog:
    """Read a catalogue and prepare the subset suitable for position matching."""

    LOG.info("Reading catalogue %s", path)
    if path.lower().endswith((".crdownload", ".part", ".tmp")):
        raise ValueError(f"Incomplete catalogue input: {path}")
    format = format or ("votable" if path.lower().endswith((".xml", ".vot")) else "fits")
    if format == "fits":
        if table_id is not None:
            raise ValueError("table_id applies only to VOTable inputs")
        # Survey jobs require explicit multi-HDU selection; legacy jobs keep
        # the original first-table convention.
        with fits.open(path) as hdul:
            if hdu is None:
                tables = [i for i, item in enumerate(hdul)
                          if isinstance(item, (fits.BinTableHDU, fits.TableHDU))]
                if not tables or (strict_table_selection and len(tables) != 1):
                    raise ValueError(f"Select an HDU for {path}: found {len(tables)} tables")
                hdu = tables[0]
            if not isinstance(hdul[hdu], (fits.BinTableHDU, fits.TableHDU)):
                raise ValueError(f"Selected HDU is not a table: {path} [{hdu}]")
        table = Table.read(path, format="fits", hdu=hdu)
        header_meta = _extract_header_metadata(path, prefix, hdu)
    elif format == "votable":
        if hdu is not None:
            raise ValueError("hdu applies only to FITS inputs")
        from astropy.io.votable import parse
        tables = list(parse(path).iter_tables())
        if table_id is None:
            if len(tables) != 1:
                raise ValueError(f"Select table_id for {path}: found {len(tables)} tables")
            selected = tables[0]
        elif isinstance(table_id, int):
            selected = tables[table_id]
        else:
            selected = next((t for t in tables if t.ID == table_id), None)
            if selected is None:
                raise ValueError(f"VOTable table_id {table_id!r} not found in {path}")
        table = selected.to_table(use_names_over_ids=True)
        header_meta = {}
    else:
        raise ValueError(f"Unsupported catalogue format: {format}")

    for name in required_columns:
        if name not in table.colnames:
            raise KeyError(f"Column '{name}' not found in {path}")

    if ra_col_name not in table.colnames:
        raise KeyError(f"Column '{ra_col_name}' not found in {path}")
    if dec_col_name not in table.colnames:
        raise KeyError(f"Column '{dec_col_name}' not found in {path}")

    ra_col = table[ra_col_name]
    dec_col = table[dec_col_name]

    default_unit = u.Unit(coordinate_unit)
    if not default_unit.is_equivalent(u.deg):
        raise ValueError(f"Coordinate unit must be angular: {coordinate_unit}")
    ra_qty = _column_to_quantity(ra_col, default_unit).to(u.deg)
    dec_qty = _column_to_quantity(dec_col, default_unit).to(u.deg)
    if ra_qty.ndim != 1 or dec_qty.ndim != 1:
        raise ValueError(f"Coordinate columns must be one-dimensional in {path}")

    finite_mask = _valid_coordinate_mask(ra_qty, dec_qty, ra_col, dec_col)
    skipped = int((~finite_mask).sum())
    if skipped:
        LOG.warning(
            "Skipping %d rows with invalid or masked coordinates in %s",
            skipped,
            os.path.basename(path),
        )

    valid_indices = np.nonzero(finite_mask)[0]
    if valid_indices.size:
        coords = SkyCoord(ra=ra_qty[finite_mask], dec=dec_qty[finite_mask], frame=frame)
    else:
        coords = SkyCoord([] * u.deg, [] * u.deg, frame=frame)

    return PreparedCatalog(
        table=table,
        coords=coords,
        valid_indices=valid_indices,
        header_meta=header_meta,
        skipped=skipped,
    )


def prepare_output_directory(output_path: str) -> Tuple[str, str, str]:
    """Ensure the destination resides in a ``Sky-CrossMatches`` directory."""

    abs_output = os.path.abspath(output_path)
    parent = os.path.dirname(abs_output)
    if not parent:
        parent = os.getcwd()

    if os.path.basename(parent) == "Sky-CrossMatches":
        out_dir = parent
    else:
        out_dir = os.path.join(parent, "Sky-CrossMatches")

    os.makedirs(out_dir, exist_ok=True)

    final_output = os.path.join(out_dir, os.path.basename(abs_output))
    log_path = os.path.join(out_dir, "xmatch_pybdsf.log")
    return final_output, log_path, out_dir


def _greedy_one_to_one(
    coords1: SkyCoord,
    coords2: SkyCoord,
    idx_map1: np.ndarray,
    idx_map2: np.ndarray,
    max_sep: Angle,
) -> List[Tuple[int, int, Angle]]:
    """Greedy global nearest-neighbour matching with 1-to-1 enforcement."""

    if coords1.size == 0 or coords2.size == 0:
        return []

    idx2, sep2d, _ = coords1.match_to_catalog_sky(coords2)
    candidates: List[Tuple[float, int, int, Angle]] = []
    for i1, (i2, sep) in enumerate(zip(idx2, sep2d)):
        if not np.isfinite(sep.to_value(u.arcsec)):
            continue
        if sep <= max_sep:
            candidates.append((sep.to_value(u.arcsec), idx_map1[i1], idx_map2[i2], sep))

    candidates.sort(key=lambda item: item[0])
    matched1: set[int] = set()
    matched2: set[int] = set()
    results: List[Tuple[int, int, Angle]] = []

    for _, orig1, orig2, sep in candidates:
        if orig1 in matched1 or orig2 in matched2:
            continue
        matched1.add(orig1)
        matched2.add(orig2)
        results.append((orig1, orig2, sep))

    return results


def _one_to_many(
    coords1: SkyCoord,
    coords2: SkyCoord,
    idx_map1: np.ndarray,
    idx_map2: np.ndarray,
    max_sep: Angle,
) -> List[Tuple[int, int, Angle]]:
    """Return every (i, j) pair within ``max_sep`` with no uniqueness constraint."""

    if coords1.size == 0 or coords2.size == 0:
        return []

    idx1, idx2, sep2d, _ = search_around_sky(coords1, coords2, max_sep)
    results: List[Tuple[int, int, Angle]] = []
    for i1, i2, sep in zip(idx1, idx2, sep2d):
        if not np.isfinite(sep.to_value(u.arcsec)):
            continue
        results.append((idx_map1[i1], idx_map2[i2], sep))
    return results


def build_output_table(
    table1: Table,
    table2: Table,
    matches: List[Tuple[int, int, Angle]],
    *,
    enrich: bool = True,
) -> Table:
    """Construct the combined match table with suffixed columns and separation."""

    output = Table()

    if matches:
        idx1 = [m[0] for m in matches]
        idx2 = [m[1] for m in matches]
        sub1 = table1[idx1]
        sub2 = table2[idx2]
    else:
        sub1 = table1[:0]
        sub2 = table2[:0]

    for name in sub1.colnames:
        output[f"{name}_1"] = sub1[name]

    for name in sub2.colnames:
        output[f"{name}_2"] = sub2[name]

    sep_arcsec = np.array([m[2].to_value(u.arcsec) for m in matches], dtype=float)
    output["sep_arcsec"] = sep_arcsec
    output["sep_arcsec"].unit = u.arcsec
    output["input_row_1"] = np.array([m[0] for m in matches], dtype=np.int64)
    output["input_row_2"] = np.array([m[1] for m in matches], dtype=np.int64)

    # Normalize standard PyBDSF measurements while retaining all original fields.
    if enrich:
        from .uncertainty import enrich_matches
        return enrich_matches(output)
    return output


def cross_match_catalogues(
    cat1: PreparedCatalog,
    cat2: PreparedCatalog,
    max_error: Angle,
    one_to_many: bool,
    *,
    enrich: bool = True,
) -> Tuple[Table, List[Tuple[int, int, Angle]]]:
    """Perform the cross-match and build the output table."""

    if not max_error.isscalar or not np.isfinite(max_error.arcsec) or max_error.arcsec <= 0:
        raise ValueError("Matching radius must be positive and finite")

    LOG.info(
        "Matching %d (valid %d) sources to %d (valid %d) within %.3f arcsec",
        len(cat1.table),
        cat1.coords.size,
        len(cat2.table),
        cat2.coords.size,
        max_error.to_value(u.arcsec),
    )

    if one_to_many:
        matches = _one_to_many(cat1.coords, cat2.coords, cat1.valid_indices, cat2.valid_indices, max_error)
    else:
        matches = _greedy_one_to_one(cat1.coords, cat2.coords, cat1.valid_indices, cat2.valid_indices, max_error)

    candidates = (matches if one_to_many else _one_to_many(
        cat1.coords, cat2.coords, cat1.valid_indices, cat2.valid_indices, max_error))
    n1 = np.bincount([m[0] for m in candidates], minlength=len(cat1.table))
    n2 = np.bincount([m[1] for m in candidates], minlength=len(cat2.table))
    match_table = build_output_table(cat1.table, cat2.table, matches, enrich=enrich)
    match_table["candidate_count_1"] = np.array([n1[m[0]] for m in matches], dtype=np.int64)
    match_table["candidate_count_2"] = np.array([n2[m[1]] for m in matches], dtype=np.int64)
    match_table["ambiguous"] = np.array([
        n1[m[0]] > 1 or n2[m[1]] > 1 for m in matches], dtype=bool)

    # Merge metadata while preserving existing meta from subtables.
    match_table.meta.update(cat1.table.meta)
    match_table.meta.update(cat2.table.meta)
    match_table.meta.update(cat1.header_meta)
    match_table.meta.update(cat2.header_meta)
    match_table.meta["HIERARCH XMATCH_MAX_ERROR_ARCSEC"] = float(max_error.to_value(u.arcsec))
    match_table.meta["HIERARCH XMATCH_MODE"] = "one-to-many" if one_to_many else "one-to-one"
    for side, cat in ((1, cat1), (2, cat2)):
        matched = {m[side - 1] for m in matches}
        ambiguous_rows = {m[side - 1] for m in candidates
                          if n1[m[0]] > 1 or n2[m[1]] > 1}
        for key, value in {
            "INPUT": len(cat.table), "VALID": len(cat.valid_indices),
            "INVALID": cat.skipped, "MATCHED": len(matched),
            "UNMATCHED": len(cat.valid_indices) - len(matched),
            "AMBIGUOUS": len(ambiguous_rows),
        }.items():
            match_table.meta[f"HIERARCH XM{side}_{key}"] = value

    return match_table, matches

FITS_CARD_LENGTH = 80
FITS_HISTORY_VALUE_LENGTH = 70


def _fits_meta_card_fits(keyword: str, value: object) -> bool:
    """Return whether a metadata item serializes to one FITS card."""

    try:
        image = fits.Card(keyword, value).image
    except (TypeError, ValueError):
        return False
    return len(image) <= FITS_CARD_LENGTH


def _history_chunks(keyword: str, value: object) -> List[str]:
    """Preserve a removed metadata item in FITS-safe HISTORY records."""

    label = keyword.removeprefix("HIERARCH ")
    text = f"{label}={value}"
    return [
        text[start : start + FITS_HISTORY_VALUE_LENGTH]
        for start in range(0, len(text), FITS_HISTORY_VALUE_LENGTH)
    ] or [f"{label}="]


def sanitize_fits_meta(tab: Table) -> Table:
    """Make table metadata safe for FITS output without losing provenance."""

    tab = tab.copy(copy_data=False)
    existing_history = tab.meta.get("HISTORY", [])
    if isinstance(existing_history, str):
        history = [existing_history]
    else:
        history = list(existing_history or [])

    moved_to_history = False
    for keyword, value in list(tab.meta.items()):
        if keyword in {"HISTORY", "COMMENT"}:
            continue
        if _fits_meta_card_fits(keyword, value):
            continue

        # Preserve a structured card where a long path can be shortened enough.
        if isinstance(value, str) and (
            keyword.endswith("_INIMAGE")
            or keyword.endswith("_CATALOG")
            or "/" in value
            or "\\" in value
        ):
            basename = os.path.basename(value)
            if basename != value and _fits_meta_card_fits(keyword, basename):
                tab.meta[keyword] = basename
                continue

        history.extend(_history_chunks(keyword, value))
        del tab.meta[keyword]
        moved_to_history = True

    if moved_to_history or existing_history:
        tab.meta["HISTORY"] = history
    return tab


def write_output(table: Table, path: str, *, overwrite: bool = True) -> None:
    """Write the matched table to ``path`` in FITS format."""

    LOG.info("Writing %d matched rows to %s", len(table), path)
    table = sanitize_fits_meta(table)
    table.write(path, format="fits", overwrite=overwrite)


def print_summary(
    matches: List[Tuple[int, int, Angle]],
    cat1: PreparedCatalog,
    cat2: PreparedCatalog,
    max_error: Angle,
    match_table: Table | None = None,
) -> None:
    """Print a concise summary of matching results to stdout."""

    matched1 = {m[0] for m in matches}
    matched2 = {m[1] for m in matches}
    total1 = len(cat1.table)
    total2 = len(cat2.table)
    unmatched1 = len(cat1.valid_indices) - len(matched1)
    unmatched2 = len(cat2.valid_indices) - len(matched2)

    summary = (
        f"Matches: {len(matches)} | cat1: input={total1} valid={len(cat1.valid_indices)} "
        f"invalid={cat1.skipped} matched={len(matched1)} unmatched={unmatched1} | "
        f"cat2: input={total2} valid={len(cat2.valid_indices)} invalid={cat2.skipped} "
        f"matched={len(matched2)} unmatched={unmatched2} | max_error={max_error.to_string(unit=u.arcsec)}"
    )
    if match_table is not None:
        summary += (f" | ambiguous rows: cat1={match_table.meta['HIERARCH XM1_AMBIGUOUS']} "
                    f"cat2={match_table.meta['HIERARCH XM2_AMBIGUOUS']}")
    print(summary)


def run_self_test() -> int:
    """Fabricate two small catalogues and exercise both matching modes."""

    LOG.info("Running self-test")

    cols = {
        "Gaus_id": [1, 2, 3, 4],
        "Isl_id": [1, 1, 2, 3],
        "Source_id": [1, 2, 3, 4],
        "Wave_id": [1, 1, 1, 1],
        "RA": [10.0, 10.00025, 150.0, 220.0],
        "DEC": [-20.0, -20.00020, 45.0, 10.0],
        "Total_flux": [10.0, 20.0, 30.0, 40.0],
        "Peak_flux": [5.0, 10.0, 15.0, 20.0],
    }
    cat1 = Table(cols)

    cols2 = {
        "Gaus_id": [10, 20, 30],
        "Isl_id": [1, 2, 3],
        "Source_id": [10, 20, 30],
        "Wave_id": [1, 1, 1],
        "RA": [10.0001, 150.0003, 310.0],
        "DEC": [-20.00005, 45.0001, -5.0],
        "Total_flux": [11.0, 31.0, 41.0],
        "Peak_flux": [6.0, 16.0, 26.0],
    }
    cat2 = Table(cols2)

    with tempfile.TemporaryDirectory() as tmpdir:
        path1 = os.path.join(tmpdir, "cat1.fits")
        path2 = os.path.join(tmpdir, "cat2.fits")
        out1 = os.path.join(tmpdir, "matches_one_to_one.fits")
        out2 = os.path.join(tmpdir, "matches_one_to_many.fits")

        cat1.write(path1, format="fits", overwrite=True)
        cat2.write(path2, format="fits", overwrite=True)

        out1_final, log_path, cross_dir = prepare_output_directory(out1)
        configure_logging(progress=True, log_path=log_path)

        args = argparse.Namespace(
            input1=path1,
            input2=path2,
            output=out1_final,
            max_error="1 arcsec",
            ra_col_1="RA",
            dec_col_1="DEC",
            ra_col_2="RA",
            dec_col_2="DEC",
            coord_frame="icrs",
            one_to_many=False,
            progress=True,
            paths_file=None,
        )
        max_error = parse_angle(args.max_error)
        cat1_prep = prepare_catalog(args.input1, args.ra_col_1, args.dec_col_1, args.coord_frame, "T1")
        cat2_prep = prepare_catalog(args.input2, args.ra_col_2, args.dec_col_2, args.coord_frame, "T2")
        table, matches = cross_match_catalogues(cat1_prep, cat2_prep, max_error, args.one_to_many)
        write_output(table, args.output)
        print_summary(matches, cat1_prep, cat2_prep, max_error, table)

        # Demonstrate one-to-many mode as well.
        LOG.info("Demonstrating one-to-many mode as well:")
        max_error_many = parse_angle("2 arcsec")
        table_many, matches_many = cross_match_catalogues(cat1_prep, cat2_prep, max_error_many, True)
        out2_final, _, _ = prepare_output_directory(out2)
        write_output(table_many, out2_final)
        print_summary(matches_many, cat1_prep, cat2_prep, max_error_many, table_many)

        LOG.info("Self-test outputs written to %s", cross_dir)

    return 0


def main(argv: Sequence[str] | None = None) -> int:
    """CLI entry point."""

    args = parse_arguments(argv)

    if args.self_test:
        return run_self_test()

    exit_code = 0
    total_jobs = len(args.jobs)

    for idx, (input1, input2, output) in enumerate(args.jobs, 1):
        output_path, log_path, output_dir = prepare_output_directory(output)
        configure_logging(args.progress, log_path=log_path)

        if total_jobs > 1:
            LOG.info("Processing job %d of %d", idx, total_jobs)
        if args.paths_file:
            LOG.info("Job paths: input1=%s | input2=%s | output=%s", input1, input2, output)

        try:
            max_error = parse_angle(args.max_error)
        except ValueError:
            LOG.exception("Failed to parse maximum error for job %d/%d", idx, total_jobs)
            exit_code = 2
            continue

        try:
            cat1 = prepare_catalog(input1, args.ra_col_1, args.dec_col_1, args.coord_frame, "T1")
            cat2 = prepare_catalog(input2, args.ra_col_2, args.dec_col_2, args.coord_frame, "T2")
        except Exception:  # pragma: no cover - CLI error propagation
            LOG.exception("Failed to prepare catalogues for job %d/%d", idx, total_jobs)
            exit_code = 1
            continue

        match_table, matches = cross_match_catalogues(cat1, cat2, max_error, args.one_to_many)

        try:
            write_output(match_table, output_path)
        except Exception:  # pragma: no cover - file system errors
            LOG.exception("Failed to write output for job %d/%d", idx, total_jobs)
            exit_code = 1
            continue

        print_summary(matches, cat1, cat2, max_error, match_table)
        LOG.info("Outputs written to %s (log: %s)", output_dir, log_path)

    return exit_code


if __name__ == "__main__":  # pragma: no cover - CLI behaviour
    sys.exit(main())
