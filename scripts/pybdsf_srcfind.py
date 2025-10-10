#!/usr/bin/env python3
"""
pybdsf_srcfind.py

Batch-capable PyBDSF launcher for local execution. Combines the single-image
source-finding logic from run_pybdsf.py with the flexible input discovery
provided by make_sbatch_jobs.py, but without any SLURM interaction.

The script accepts individual FITS images on the command line or a text file
listing image paths and their associated frequencies (MHz). Results for each
image are written under <image_dir>/pybdsf.results/<base>/ along with a
timestamped log file capturing the processing details.

Virtual environment quick start:
# python3 -m venv venv
# source venv/bin/activate
# pip install astropy numpy scipy pybdsf
# deactivate
"""

from __future__ import annotations

import argparse
import glob
import logging
import os
import re
import shlex
import sys
import tempfile
import traceback
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


def _infer_freq_mhz_from_name(name: str) -> Optional[float]:
    """
    Extract a frequency (MHz) from filename, if possible.
    Patterns:
      - '1299~1336MHz'  -> mean of range
      - '1430MHz'
      - tokens like '_1420_' (assumed MHz)
    """
    n = os.path.basename(name)

    # Range like 1299~1336MHz
    m = re.search(r"(\d+(?:\.\d+)?)\s*~\s*(\d+(?:\.\d+)?)\s*MHz", n, re.I)
    if m:
        a, b = float(m.group(1)), float(m.group(2))
        return 0.5 * (a + b)

    # Single like 1430MHz
    m = re.search(r"(\d+(?:\.\d+)?)\s*MHz", n, re.I)
    if m:
        return float(m.group(1))

    # Common tokens: _1420_ etc. (assume MHz)
    m = re.search(r"[_-](\d{3,4})(?:\.\d+)?[_-]", n)
    if m:
        v = float(m.group(1))
        if 100.0 <= v <= 6000.0:
            return v

    return None


def _find_header_freq_hz(hdr) -> Optional[float]:
    """
    Find frequency (Hz) in common header keys.
    """
    for k in ["RESTFRQ", "RESTFRQ1", "FREQ", "CFREQ", "CRVAL3"]:
        if k in hdr:
            try:
                return float(hdr[k])
            except Exception:
                pass

    # FREQ as primary axis (rare for 2D images)
    try:
        if str(hdr.get("CTYPE1", "")).upper().startswith("FREQ") and "CRVAL1" in hdr:
            return float(hdr["CRVAL1"])
    except Exception:
        pass

    return None


def _inject_freq_to_header(hdr, freq_hz: float) -> None:
    """
    Inject a minimal, standards-ish set of cards so old PyBDSF can see frequency.
    Does NOT change NAXIS (keeps image 2D), but declares a 3rd WCS axis via WCSAXES.
    """
    hdr["RESTFRQ"] = float(freq_hz)
    try:
        hdr.comments["RESTFRQ"] = "Injected by pybdsf_srcfind.py (Hz)"
    except Exception:
        pass
    hdr["FREQ"] = float(freq_hz)
    hdr["CFREQ"] = float(freq_hz)

    try:
        wcsaxes = int(hdr.get("WCSAXES", 2))
    except Exception:
        wcsaxes = 2
    if wcsaxes < 3:
        hdr["WCSAXES"] = 3

    hdr["CTYPE3"] = "FREQ"
    hdr["CUNIT3"] = "Hz"
    hdr["CRVAL3"] = float(freq_hz)
    hdr["CRPIX3"] = 1.0
    hdr["CDELT3"] = 1.0
    hdr["SPECSYS"] = hdr.get("SPECSYS", "TOPOCENT")


def is_fits_path(p: Path) -> bool:
    return p.suffix.lower() in {".fits", ".fit", ".fts"}


def expand_images(patterns: Optional[Iterable[str]]) -> List[str]:
    """
    Accept:
      - directory paths (recurse for FITS)
      - glob patterns (e.g., 'subband_images/*highband*.fits')
    Return ordered, de-duplicated absolute paths.
    """
    if not patterns:
        return []

    found: List[str] = []
    for item in patterns:
        p = Path(item).expanduser()
        if p.exists() and p.is_dir():
            for fp in sorted(p.rglob("*")):
                if fp.is_file() and is_fits_path(fp):
                    found.append(str(fp.resolve()))
        else:
            for fp in map(Path, glob.glob(str(p))):
                if fp.is_file() and is_fits_path(fp):
                    found.append(str(fp.resolve()))

    seen, ordered = set(), []
    for ap in found:
        if ap not in seen:
            seen.add(ap)
            ordered.append(ap)
    return ordered


def _parse_frequency_token(token: str) -> float:
    """Return frequency in Hz parsed from token."""
    token = token.strip()
    if not token:
        raise ValueError("Empty frequency token")

    pattern = re.compile(r"^\s*([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)\s*([A-Za-z]*)\s*$")
    m = pattern.match(token)
    if not m:
        raise ValueError("Could not parse frequency token {!r}".format(token))

    value = float(m.group(1))
    suffix = m.group(2).lower()

    if suffix in ("", "mhz"):
        return value * 1.0e6
    if suffix == "hz":
        return value
    if suffix == "ghz":
        return value * 1.0e9
    if suffix == "khz":
        return value * 1.0e3

    raise ValueError("Unsupported frequency unit {!r}".format(suffix))


def parse_filelist(filelist: Optional[str]) -> List[Tuple[str, Optional[float]]]:
    """
    Read a text file where each non-empty line provides:
        <image_path> <frequency>
    Frequency tokens default to MHz but unit suffixes (MHz/Hz/GHz/kHz) are accepted.
    Lines starting with '#' are ignored. Paths are expanded and resolved.
    """
    if not filelist:
        return []

    entries: List[Tuple[str, Optional[float]]] = []
    path_obj = Path(filelist).expanduser()
    if not path_obj.exists():
        raise FileNotFoundError("File list not found: {}".format(filelist))

    with path_obj.open("r") as handle:
        for idx, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue

            tokens = shlex.split(line.replace(",", " "))
            if len(tokens) < 2:
                raise ValueError(
                    "Expected at least two columns (path and frequency) on line {} of {}".format(
                        idx, filelist
                    )
                )

            img_path = Path(tokens[0]).expanduser().resolve()
            freq_hz: Optional[float] = None
            try:
                freq_hz = _parse_frequency_token(tokens[1])
            except ValueError:
                raise ValueError(
                    "Could not parse frequency on line {} of {}: {!r}".format(
                        idx, filelist, tokens[1]
                    )
                )

            entries.append((str(img_path), freq_hz))

    return entries


def build_processing_list(
    images: Optional[Iterable[str]],
    filelist: Optional[str],
) -> Tuple[List[str], Dict[str, float]]:
    """
    Combine CLI images and filelist entries. Returns:
      - ordered list of absolute image paths (deduplicated, preserving filelist order first)
      - mapping of image path -> frequency (Hz) from the filelist
    """
    expanded_from_cli = expand_images(images)
    file_entries = parse_filelist(filelist)

    path_to_freq: Dict[str, float] = {}
    ordered_paths: List[str] = []

    for path, freq in file_entries:
        if path not in path_to_freq:
            ordered_paths.append(path)
        path_to_freq[path] = freq  # later occurrences override

    for path in expanded_from_cli:
        if path not in ordered_paths:
            ordered_paths.append(path)

    return ordered_paths, path_to_freq


def _setup_logger(log_path: str) -> logging.Logger:
    logger = logging.getLogger("pybdsf_srcfind.{}".format(os.path.basename(log_path)))
    logger.setLevel(logging.INFO)
    logger.propagate = False

    formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")

    file_handler = logging.FileHandler(log_path)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)

    logger.info("Log file initialised at %s", datetime.utcnow().isoformat() + "Z")
    return logger


def _close_logger(logger: logging.Logger) -> None:
    handlers = list(logger.handlers)
    for handler in handlers:
        handler.flush()
        handler.close()
        logger.removeHandler(handler)


def find_sources(
    fits_in: str,
    *,
    base: Optional[str] = None,
    isl_thresh: float = 3.0,
    pix_thresh: Optional[float] = None,
    snrthreshold: float = 10.0,
    fluxthreshold: float = 1.0e-6,
    pblimit: float = 0.05,
    freq_hz: Optional[float] = None,
    logger: Optional[logging.Logger] = None,
) -> Tuple[str, str]:
    """
    Run PyBDSF on a single FITS image and export GAUL catalog, model and RMS images.
    Returns (outdir, log_path).
    """
    import os

    try:
        from astropy.io import fits
        from astropy.wcs import WCS
        import bdsf
    except Exception as e:
        if logger:
            logger.error("Missing dependencies: %s", e)
        else:
            print("ERROR: Missing dependencies:", e, file=sys.stderr)
        raise

    if not os.path.isfile(fits_in):
        raise IOError("FITS not found: {0}".format(fits_in))

    imagedir = os.path.dirname(os.path.abspath(fits_in))
    fitsfile = os.path.basename(fits_in)
    stem = os.path.splitext(fitsfile)[0]

    if base is None:
        base = stem[:10]
    if pix_thresh is None:
        pix_thresh = 1.5 * float(isl_thresh)

    out0 = os.path.join(imagedir, "pybdsf.results")
    outdir = os.path.join(out0, base)
    os.makedirs(outdir, exist_ok=True)

    timestamp = datetime.utcnow().strftime("%Y%m%d-%H%M%S")
    log_path = os.path.join(outdir, f"{base}_pybdsf_srcfind_{timestamp}.log")

    close_logger = False
    if logger is None:
        logger = _setup_logger(log_path)
        close_logger = True
    else:
        log_path = getattr(logger, "log_path", log_path)
    setattr(logger, "log_path", log_path)

    logger.info("Starting PyBDSF processing for %s", fits_in)
    logger.info("Output directory: %s", outdir)

    hdr_freq_hz = None
    tmp_path = None

    try:
        with fits.open(fits_in) as hdul:
            hdr_freq_hz = _find_header_freq_hz(hdul[0].header)
    except Exception:
        logger.warning("Could not read header frequency, will try other methods.", exc_info=True)

    if freq_hz is None:
        if hdr_freq_hz is not None:
            freq_hz = hdr_freq_hz
        else:
            mhz = _infer_freq_mhz_from_name(fits_in)
            if mhz is not None:
                freq_hz = mhz * 1.0e6

    if freq_hz is None:
        raise RuntimeError(
            "No frequency in header and none provided. "
            "Pass --freq-mhz/--freq-hz or ensure RESTFRQ is in the FITS header."
        )

    logger.info("Using frequency: %.3f MHz", freq_hz / 1.0e6)

    use_file = fitsfile
    run_dir = imagedir
    wrote_temp = False
    try:
        with fits.open(fits_in) as hdul:
            if _find_header_freq_hz(hdul[0].header) is None:
                _inject_freq_to_header(hdul[0].header, freq_hz)

                stem = os.path.splitext(os.path.basename(fits_in))[0]
                safe_stem = re.sub(r"[^A-Za-z0-9._-]+", "_", stem)[:60]
                tmp_dir = outdir if os.access(outdir, os.W_OK) else imagedir
                tmp_fd, tmp_path = tempfile.mkstemp(
                    prefix=safe_stem + ".", suffix=".freqfix.fits", dir=tmp_dir
                )
                os.close(tmp_fd)

                hdul.writeto(tmp_path, overwrite=True)
                use_file = os.path.basename(tmp_path)
                run_dir = os.path.dirname(tmp_path)
                wrote_temp = True
                logger.info("Wrote temp FITS with injected RESTFRQ to %s", tmp_path)
    except Exception:
        logger.warning(
            "Could not write temp FITS copy with frequency keywords; "
            "falling back to direct freq injection.",
            exc_info=True,
        )

    cwd = os.getcwd()
    try:
        os.chdir(run_dir)
        kwargs = dict(
            output_opts=True,
            output_all=True,
            shapelet_do=True,
            solnname=base,
            thresh="hard",
            thresh_isl=float(isl_thresh),
            thresh_pix=float(pix_thresh),
            adaptive_rms_box=True,
        )
        if not wrote_temp:
            kwargs["freq"] = float(freq_hz)

        logger.info("Invoking bdsf.process_image on %s", use_file)
        img = bdsf.process_image(use_file, **kwargs)

        modimage = os.path.join(outdir, base + "-modimage.fits")
        rmsfile = os.path.join(outdir, base + "-rms.fits")
        srcfilename = os.path.join(outdir, base + "-source-cat.txt")

        img.export_image(outfile=modimage, img_format="fits", img_type="gaus_model", clobber=True)
        img.export_image(outfile=rmsfile, img_format="fits", img_type="rms", clobber=True)
        img.write_catalog(outfile=srcfilename, clobber=True, format="ascii", catalog_type="gaul")

        logger.info("Exported model image to %s", modimage)
        logger.info("Exported RMS image to %s", rmsfile)
        logger.info("Exported GAUL catalog to %s", srcfilename)

        try:
            from astropy.io import fits as _fits
            from astropy.wcs import WCS as _WCS

            hdul = _fits.open(os.path.join(run_dir, use_file))
            _ = _WCS(hdul[0].header)
            data = hdul[0].data
            imax = None
            try:
                imax = data.max() if data is not None else None
            except Exception:
                imax = None
            hdul.close()
            logger.info("Sanity check: max pixel value %.3g", imax if imax is not None else float("nan"))
        except Exception:
            logger.warning("Sanity check failed", exc_info=True)

    except Exception as ex:
        logger.error("Error while running PyBDSF: %s", ex)
        logger.debug("Traceback:\n%s", traceback.format_exc())
        raise
    finally:
        os.chdir(cwd)
        if wrote_temp and tmp_path and os.path.exists(tmp_path):
            try:
                os.remove(tmp_path)
                logger.debug("Removed temporary FITS %s", tmp_path)
            except Exception:
                logger.warning("Could not remove temporary FITS %s", tmp_path, exc_info=True)
        if close_logger:
            _close_logger(logger)

    return outdir, log_path


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run PyBDSF locally on FITS images with optional per-image frequency overrides."
    )
    parser.add_argument(
        "--images",
        nargs="*",
        default=None,
        help=(
            "FITS image paths, directories, or glob patterns. "
            "Directories are traversed recursively."
        ),
    )
    parser.add_argument(
        "--filelist",
        default=None,
        help=(
            "Text file with rows '<image_path> <freq>' (MHz by default). "
            "Frequencies may include unit suffixes."
        ),
    )
    parser.add_argument(
        "--isl",
        type=float,
        default=3.0,
        help="PyBDSF island threshold (default 3.0).",
    )
    parser.add_argument(
        "--pix",
        type=float,
        default=None,
        help="PyBDSF pixel threshold (defaults to 1.5 * isl).",
    )
    parser.add_argument(
        "--snr",
        type=float,
        default=10.0,
        help="Kept for parity; not directly used by PyBDSF.",
    )
    parser.add_argument(
        "--flux",
        type=float,
        default=1.0e-6,
        help="Kept for parity with legacy scripts.",
    )
    parser.add_argument(
        "--pblimit",
        type=float,
        default=0.05,
        help="Kept for parity with legacy scripts.",
    )
    parser.add_argument(
        "--freq-mhz",
        type=float,
        default=None,
        help="Fallback frequency in MHz applied to images lacking other overrides.",
    )
    parser.add_argument(
        "--freq-hz",
        type=float,
        default=None,
        help="Fallback frequency in Hz (overrides --freq-mhz).",
    )
    parser.add_argument(
        "--base-prefix",
        default=None,
        help=(
            "Optional prefix to prepend to the first 10 characters of each FITS stem "
            "for the output directory name."
        ),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="List the images that would be processed and exit.",
    )
    return parser.parse_args(argv)


def determine_frequency(
    image_path: str,
    freq_overrides: Dict[str, float],
    args: argparse.Namespace,
) -> Optional[float]:
    """
    Determine the frequency (Hz) for the image, preferring explicit overrides.
    """
    if image_path in freq_overrides:
        return freq_overrides[image_path]
    if args.freq_hz is not None:
        return float(args.freq_hz)
    if args.freq_mhz is not None:
        return float(args.freq_mhz) * 1.0e6
    return None


def compute_base_name(stem: str, prefix: Optional[str]) -> str:
    base = stem[:10]
    if prefix:
        base = f"{prefix}{base}"
    return base


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)

    try:
        ordered_paths, freq_overrides = build_processing_list(args.images, args.filelist)
    except Exception as exc:
        print("ERROR: {}".format(exc), file=sys.stderr)
        return 2

    if not ordered_paths:
        print("No FITS images found; provide --images or --filelist.", file=sys.stderr)
        return 1

    if args.dry_run:
        print("Dry run: the following images would be processed:")
        for path in ordered_paths:
            freq = determine_frequency(path, freq_overrides, args)
            if freq is None:
                freq_str = "auto/header"
            else:
                freq_str = "{:.3f} MHz".format(freq / 1.0e6)
            print("  {} (freq {})".format(path, freq_str))
        return 0

    exit_code = 0
    for image_path in ordered_paths:
        freq_hz = determine_frequency(image_path, freq_overrides, args)
        stem = os.path.splitext(os.path.basename(image_path))[0]
        base = compute_base_name(stem, args.base_prefix)

        try:
            outdir, log_path = find_sources(
                image_path,
                base=base,
                isl_thresh=args.isl,
                pix_thresh=args.pix,
                snrthreshold=args.snr,
                fluxthreshold=args.flux,
                pblimit=args.pblimit,
                freq_hz=freq_hz,
            )
            print("Completed {} -> outputs in {} (log: {})".format(image_path, outdir, log_path))
        except Exception as exc:
            print("FAILED {}: {}".format(image_path, exc), file=sys.stderr)
            exit_code = 3

    return exit_code


if __name__ == "__main__":
    sys.exit(main())
