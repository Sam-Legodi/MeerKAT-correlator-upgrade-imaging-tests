#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function

"""
Imaging helper that runs `tclean` on each MeasurementSet twice (nominal L-band
and UHF-band sub-ranges, with automatic quartile fallbacks) and optionally per
scan. If an output product already exists the script skips imaging and reruns
QA on the available FITS instead.

Outputs (per MS):
  • `<msdir>/images/<msbase>_lowband.*` and `_highband.*` CASA images
  • Matching FITS exports plus `_scan<N>` variants when per-scan imaging is on
  • `<…>_qa.txt` summaries for every FITS file created or detected

How to run (CASA 5/6):
  `casa --nologger --log2term -c tclean_two_bands.py [--scans=1,2,5] <ms1> <ms2> …`
  or set `MSFILE` in the script to image a single dataset without CLI args.
"""

import os
import sys
import math
import datetime
import numpy as np

# Try CASA 6 task imports; fall back to CASA 5 globals if present.
try:
    from casatasks import tclean, exportfits
except Exception:
    pass

# -------- CASA tools (NO msmd used) --------
# Measurement Set tool (for selection row-count sanity checks)
try:
    from casatools import ms as _ms_cls
    _MS_TOOL_AVAILABLE = True
except Exception:
    try:
        from taskinit import mstool as _ms_cls  # CASA 5
        _MS_TOOL_AVAILABLE = True
    except Exception:
        _MS_TOOL_AVAILABLE = False
        _ms_cls = None

# Table tool (for FIELD/SPECTRAL_WINDOW/MAIN reads)
try:
    from casatools import table as _tb_cls
    _TB_TOOL_AVAILABLE = True
except Exception:
    try:
        from taskinit import tbtool as _tb_cls  # CASA 5
        _TB_TOOL_AVAILABLE = True
    except Exception:
        _TB_TOOL_AVAILABLE = False
        _tb_cls = None

# Prefer Astropy if available (CASA 6 bundles it)
try:
    from astropy.io import fits
except Exception:
    fits = None

# Target bands in Hz (nominal requests)
LOW_LO_HZ_REQ  = int(round(8.98e8))
LOW_HI_HZ_REQ  = int(round(1.00e9))
HIGH_LO_HZ_REQ = int(round(1.46e9))
HIGH_HI_HZ_REQ = int(round(1.70e9))

# Default imaging parameters (edit as needed)
MSFILE     = "data/1757723806_sdp_l0.full.SDP_Allcorrected.ms"  # if non-empty, CLI args are ignored
# MSFILE     = "data/reference_obs/1692135074_sdp_l0.full.SDP.Allcorrected.ms"
# MSFILE     = "data/1757791280_sdp_l0.full.SDP_Allcorrected.ms"

FIELDNAME  = 'J2147-8132'
CELL       = '1.2arcsec'
DATACOL    = "DATA"  # "corrected" or "data"
IMSIZE     = 4096
NITER      = 10000
STOKES     = 'I'
GAIN       = 0.1
DECONV     = 'clark'
THRESH     = 1e-6
GRIDDER    = 'wproject'
WPROJ      = 64
SPECMODE   = 'mfs'
WEIGHTING  = 'briggs'
ROBUST     = 0.0
PBLIMIT    = -1

# Removes intermediate tclean products (psf, model, residual, mask, etc) after imaging.
CLEANUPPRODUCTS = False

# Per-scan imaging toggle
IMAGE_SCANS     = True  # leave this True to enable per-scan imaging


def ensure_dir(path):
    if not os.path.isdir(path):
        os.makedirs(path)

def imagename_for(msfile, band_suffix, scan_id=None):
    msdir  = os.path.dirname(os.path.abspath(msfile))
    msbase = os.path.splitext(os.path.basename(msfile))[0]
    outdir = os.path.join(msdir, 'images')
    ensure_dir(outdir)
    if scan_id is None:
        tag = '{}_{}'.format(msbase, band_suffix)
    else:
        tag = '{}_{}_scan{:d}'.format(msbase, band_suffix, int(scan_id))
    return os.path.join(outdir, tag)

def _pixel_scale_arcsec_from_header(h):
    to_arcsec = 3600.0
    try:
        if 'CDELT1' in h and 'CDELT2' in h:
            dx = abs(float(h['CDELT1'])) * to_arcsec
            dy = abs(float(h['CDELT2'])) * to_arcsec
            return dx, dy
        cd11 = float(h.get('CD1_1', np.nan))
        cd12 = float(h.get('CD1_2', np.nan))
        cd21 = float(h.get('CD2_1', np.nan))
        cd22 = float(h.get('CD2_2', np.nan))
        if not np.isnan([cd11, cd12, cd21, cd22]).any():
            dx = math.hypot(cd11, cd12) * to_arcsec
            dy = math.hypot(cd21, cd22) * to_arcsec
            return dx, dy
    except Exception:
        pass
    return float('nan'), float('nan')

def _beam_from_header(h):
    try:
        bmaj = float(h.get('BMAJ', np.nan)) * 3600.0  # deg -> arcsec
        bmin = float(h.get('BMIN', np.nan)) * 3600.0
        bpa  = float(h.get('BPA',  np.nan))
        return bmaj, bmin, bpa
    except Exception:
        return float('nan'), float('nan'), float('nan')

def qa_fits(fits_path):
    """Quick QA on a FITS image. Writes <fits_prefix>_qa.txt."""
    report_path = os.path.splitext(fits_path)[0] + "_qa.txt"
    now = datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S UTC")

    # --- Try modern FITS I/O first (astropy), then legacy (pyfits) ---
    _fits = None
    if 'fits' in globals() and fits is not None:
        _fits = fits
    else:
        try:
            from astropy.io import fits as _fits  # CASA 6 often bundles astropy
        except Exception:
            try:
                import pyfits as _fits  # legacy name on some CASA 5 installs
            except Exception:
                _fits = None

    if _fits is not None:
        try:
            with _fits.open(fits_path, memmap=False) as hdul:
                hdu = None
                for x in hdul:
                    if getattr(x, 'data', None) is not None:
                        hdu = x
                        break
                if hdu is None or hdu.data is None:
                    raise RuntimeError("no image data")
                data = np.asarray(hdu.data, dtype=float)
                hdr  = hdu.header

                arr = np.ravel(data)
                finite = np.isfinite(arr)
                arr_f = arr[finite]
                n_total = arr.size
                n_finite = arr_f.size
                n_nans = n_total - n_finite
                nan_frac = float(n_nans) / float(n_total) if n_total > 0 else float('nan')

                if n_finite == 0:
                    raise RuntimeError("all values are NaN or inf")

                median = float(np.nanmedian(arr_f))
                mad = float(np.nanmedian(np.abs(arr_f - median)))
                robust_rms = 1.4826 * mad
                mean = float(np.nanmean(arr_f))
                std  = float(np.nanstd(arr_f))
                vmin = float(np.nanmin(arr_f))
                vmax = float(np.nanmax(arr_f))
                p5, p50, p95, p99 = [float(x) for x in np.nanpercentile(arr_f, [5, 50, 95, 99])]
                peak = vmax
                dyn_range = abs(peak) / robust_rms if robust_rms > 0 else float('inf')

                dx_as, dy_as = _pixel_scale_arcsec_from_header(hdr)
                bmaj_as, bmin_as, bpa_deg = _beam_from_header(hdr)
                unit = hdr.get('BUNIT', '')

            lines = [
                "QA report",
                "Time: {}".format(now),
                "File: {}".format(fits_path),
                "Shape: {}  dtype: {}".format(data.shape, data.dtype),
                "Unit: {}".format(unit),
                "Pixels: total={} finite={} NaNs={} nan_frac={:.6f}".format(n_total, n_finite, n_nans, nan_frac),
                "Stats: min={:.6g} max={:.6g} mean={:.6g} median={:.6g} std={:.6g}".format(vmin, vmax, mean, median, std),
                "Robust: mad={:.6g} robust_rms={:.6g}".format(mad, robust_rms),
                "Percentiles: p5={:.6g} p50={:.6g} p95={:.6g} p99={:.6g}".format(p5, p50, p95, p99),
                "DynamicRange: peak/robust_rms = {:.6g}".format(dyn_range),
                "PixelScale: dx={:.3f}\" dy={:.3f}\"".format(dx_as, dy_as),
                "Beam: BMAJ={:.3f}\" BMIN={:.3f}\" BPA={:.2f} deg".format(bmaj_as, bmin_as, bpa_deg),
            ]
            with open(report_path, "w") as f:
                f.write("\n".join(lines) + "\n")

            print("QA: wrote {}".format(report_path))
            return report_path

        except Exception as e:
            # Fall through to CASA-based QA on the .image sibling.
            print("QA(FITS) failed on {}: {}. Falling back to CASA image.".format(fits_path, e))

    # --- CASA-based fallback on the sibling .image (no FITS I/O needed) ---
    casa_image = os.path.splitext(fits_path)[0] + ".image"
    if not os.path.exists(casa_image):
        with open(report_path, "w") as f:
            f.write("QA report\nTime: {}\nFITS I/O unavailable and CASA image missing.\nFile: {}\n"
                    .format(now, fits_path))
        print("QA: fallback failed; no CASA image at {}.".format(casa_image))
        return report_path

    # Prefer casatasks.imstat/imhead to avoid tool/version pitfalls
    try:
        from casatasks import imstat, imhead
    except Exception:
        imstat = None
        imhead = None

    stats = {}
    beam_bmaj = beam_bmin = beam_bpa = float('nan')
    dx_as = dy_as = float('nan')
    unit = ""

    try:
        if imstat is not None:
            s = imstat(imagename=casa_image, chans="", stokes="")
            # imstat returns dict with keys like 'min', 'max', 'mean', 'sigma', 'median', 'medabsdevmed', 'npts'
            stats = s or {}
        if imhead is not None:
            h = imhead(imagename=casa_image, mode='get', hdkey='')
            # Pixel scale: use cdelt1/2 (radians) -> arcsec, or CDELT in degrees depending on build
            try:
                cdelt = h.get('cdelt', None) or [h.get('cdelt1', None), h.get('cdelt2', None)]
                if cdelt and len(cdelt) >= 2:
                    # cdelt may be in radians; detect via magnitude
                    c1, c2 = float(cdelt[0]), float(cdelt[1])
                    # if |c| < 1e-3 assume radians, else degrees
                    to_arcsec = 206264.806 if (abs(c1) < 1e-3 and abs(c2) < 1e-3) else 3600.0
                    dx_as, dy_as = abs(c1) * to_arcsec, abs(c2) * to_arcsec
            except Exception:
                pass
            try:
                bmaj = h.get('bmaj', None); bmin = h.get('bmin', None); bpa = h.get('bpa', None)
                if bmaj is not None and bmin is not None:
                    beam_bmaj = float(bmaj) * 3600.0  # deg -> arcsec
                    beam_bmin = float(bmin) * 3600.0
                    beam_bpa  = float(bpa) if bpa is not None else float('nan')
            except Exception:
                pass
            try:
                bunit = h.get('bunit', None)
                if bunit:
                    unit = str(bunit)
            except Exception:
                pass
    except Exception as e:
        print("QA(CASA fallback) read error on {}: {}".format(casa_image, e))

    # Build robust-ish metrics from imstat result
    try:
        arr_min = float(stats.get('min', [float('nan')])[0])
        arr_max = float(stats.get('max', [float('nan')])[0])
        mean    = float(stats.get('mean', [float('nan')])[0])
        std     = float(stats.get('sigma', [float('nan')])[0])
        median  = float(stats.get('median', [float('nan')])[0])
        mad     = float(stats.get('medabsdevmed', [float('nan')])[0])
        robust_rms = 1.4826 * mad if math.isfinite(mad) else std
        dyn_range  = (abs(arr_max) / robust_rms) if (robust_rms and robust_rms > 0) else float('inf')
        n_total    = int(stats.get('npts', 0))
    except Exception:
        arr_min = arr_max = mean = std = median = mad = robust_rms = dyn_range = float('nan')
        n_total = 0

    lines = [
        "QA report (CASA fallback)",
        "Time: {}".format(now),
        "File: {}".format(fits_path),
        "Stats: min={:.6g} max={:.6g} mean={:.6g} median={:.6g} std={:.6g}".format(arr_min, arr_max, mean, median, std),
        "Robust: mad={:.6g} robust_rms={:.6g}".format(mad, robust_rms),
        "DynamicRange: peak/robust_rms = {:.6g}".format(dyn_range),
        "PixelScale: dx={:.3f}\" dy={:.3f}\"".format(dx_as, dy_as),
        "Beam: BMAJ={:.3f}\" BMIN={:.3f}\" BPA={:.2f} deg".format(beam_bmaj, beam_bmin, beam_bpa),
        "Unit: {}".format(unit),
        "Note: FITS I/O not available; statistics computed from CASA image: {}".format(casa_image),
    ]
    with open(report_path, "w") as f:
        f.write("\n".join(lines) + "\n")

    print("QA: wrote (CASA fallback) {}".format(report_path))
    return report_path


# ---------- MS utilities (NO msmd) ----------

def _read_spw_freq_minmax(msfile):
    """
    Read per-SPW min/max frequencies (Hz) from SPECTRAL_WINDOW/CHAN_FREQ using the table tool.
    Returns list of (fmin, fmax) per spw index. On failure, returns [].
    """
    if not _TB_TOOL_AVAILABLE:
        return []
    tb = None
    spw_path = os.path.join(msfile, 'SPECTRAL_WINDOW')
    try:
        tb = _tb_cls()
        tb.open(spw_path)
        chan_freq_col = tb.getcol('CHAN_FREQ')  # (nchan, nspw) or vector per row
        if chan_freq_col is None:
            return []
        arr = np.array(chan_freq_col, dtype=float)
        if arr.ndim == 2:
            arr = np.transpose(arr)  # (nspw, nchan)
        elif arr.ndim == 1:
            arr = arr[None, :]
        spw_minmax = []
        for spw_idx in range(arr.shape[0]):
            freqs = np.asarray(arr[spw_idx], dtype=float)
            if freqs.size == 0 or not np.isfinite(freqs).any():
                spw_minmax.append((None, None))
            else:
                spw_minmax.append((float(np.nanmin(freqs)), float(np.nanmax(freqs))))
        return spw_minmax
    except Exception:
        return []
    finally:
        try:
            tb.close()
        except Exception:
            pass

def _ms_total_freq_bounds(msfile):
    """Return (fmin, fmax) across all SPWs via table reads. None if unavailable."""
    spw_minmax = _read_spw_freq_minmax(msfile)
    vals = [(a, b) for (a, b) in spw_minmax if (a is not None and b is not None)]
    if not vals:
        return (None, None)
    fmins = [a for a, _ in vals]
    fmaxs = [b for _, b in vals]
    return (float(min(fmins)), float(max(fmaxs)))

def _band_exists(msfile, lo_hz, hi_hz):
    """True if any SPW overlaps [lo_hz, hi_hz] using table-based SPW bounds."""
    if lo_hz >= hi_hz:
        return False
    spw_minmax = _read_spw_freq_minmax(msfile)
    found = False
    for fmin, fmax in spw_minmax:
        if fmin is None or fmax is None:
            continue
        lo = max(lo_hz, fmin)
        hi = min(hi_hz, fmax)
        if lo < hi:
            found = True
            break
    return found if spw_minmax else True

def _quartile_bands(msfile):
    """Lower 25% and upper 25% of total frequency span (table-based)."""
    fmin, fmax = _ms_total_freq_bounds(msfile)
    if fmin is None or fmax is None or not np.isfinite([fmin, fmax]).all():
        return (LOW_LO_HZ_REQ, LOW_HI_HZ_REQ), (HIGH_LO_HZ_REQ, HIGH_HI_HZ_REQ)
    span = fmax - fmin
    if span <= 0:
        return (LOW_LO_HZ_REQ, LOW_HI_HZ_REQ), (HIGH_LO_HI_HZ_REQ, HIGH_HI_HZ_REQ)
    q = 0.25 * span
    lowband  = (int(round(fmin)),     int(round(fmin + q)))
    highband = (int(round(fmax - q)), int(round(fmax)))
    return lowband, highband

def decide_bands(msfile):
    """Choose effective low/high bands. Fallback to quartiles if either requested band is missing."""
    req_low  = (LOW_LO_HZ_REQ,  LOW_HI_HZ_REQ)
    req_high = (HIGH_LO_HZ_REQ, HIGH_HI_HZ_REQ)
    low_ok   = _band_exists(msfile, req_low[0],  req_low[1])
    high_ok  = _band_exists(msfile, req_high[0], req_high[1])
    if low_ok and high_ok:
        return req_low, req_high
    qb_low, qb_high = _quartile_bands(msfile)
    print("Band fallback in effect for {}: using quartiles.".format(msfile))
    print("  lowband:  {} Hz ~ {} Hz".format(qb_low[0],  qb_low[1]))
    print("  highband: {} Hz ~ {} Hz".format(qb_high[0], qb_high[1]))
    return qb_low, qb_high

def _spw_overlap_selector(msfile, lo_hz, hi_hz):
    """
    Build an SPW selection that intersects [lo_hz, hi_hz] with actual SPW coverages
    using table-based SPW min/max. If we cannot read SPWs, fall back to '*:lo~hiHz'.
    """
    spw_minmax = _read_spw_freq_minmax(msfile)
    if not spw_minmax:
        return "*:{0}~{1}Hz".format(int(lo_hz), int(hi_hz))
    chunks = []
    for spw, (fmin, fmax) in enumerate(spw_minmax):
        if fmin is None or fmax is None:
            continue
        lo = max(float(lo_hz), fmin)
        hi = min(float(hi_hz), fmax)
        if lo < hi:
            chunks.append("{0}:{1}Hz~{2}Hz".format(int(spw), int(lo), int(hi)))
    return ",".join(chunks) if chunks else "*:{0}~{1}Hz".format(int(lo_hz), int(hi_hz))

# ---------- Discover scans for a FIELD (NO msmd) ----------

def _field_ids_for_name(msfile, fieldname):
    """
    Return a list of FIELD_IDs whose NAME equals fieldname.
    """
    if not _TB_TOOL_AVAILABLE:
        return []
    tb = None
    try:
        tb = _tb_cls()
        tb.open(os.path.join(msfile, 'FIELD'))
        names = tb.getcol('NAME')  # vector of strings
        ids = []
        for idx, nm in enumerate(names):
            try:
                nm_str = nm.decode() if hasattr(nm, 'decode') else str(nm)
            except Exception:
                nm_str = str(nm)
            if nm_str == fieldname:
                ids.append(int(idx))
        return ids
    except Exception:
        return []
    finally:
        try:
            tb.close()
        except Exception:
            pass

def discover_scans_for_field(msfile, fieldname, chunk_rows=500000):
    """
    Efficiently discover unique SCAN_NUMBERs for a given field name.

    Strategy (most robust first):
      1) Try msmetadata (msmd): fieldsforname -> scansforfield
      2) Try ms tool: msselect({'field': <name>}) -> getscansummary()
      3) Fallback to table reads: FIELD -> MAIN (chunked), match FIELD_ID then collect SCAN_NUMBER

    Returns a sorted list of int scan IDs (may be empty if none found).
    """
    # -------------- Attempt 1: msmetadata (if available) --------------
    try:
        msmd = _get_msmd_tool()  # must exist in the surrounding module
    except Exception:
        msmd = None

    if msmd is not None:
        try:
            msmd.open(msfile)
            # fieldsforname(name) usually returns an array of field IDs; some CASA builds use a list.
            try:
                fids = msmd.fieldsforname(fieldname)
            except Exception:
                fids = []
            if fids is None:
                fids = []
            # Normalize to ints
            fids = [int(x) for x in (list(fids) if hasattr(fids, '__iter__') else [fids])]
            scans = set()
            for fid in fids:
                try:
                    s = msmd.scansforfield(int(fid))
                    if s is None:
                        continue
                    for sc in (list(s) if hasattr(s, '__iter__') else [s]):
                        scans.add(int(sc))
                except Exception:
                    # Some CASA builds may not have scansforfield; ignore and continue
                    pass
            try:
                msmd.close()
            except Exception:
                pass
            if scans:
                out = sorted(scans)
                print("Discovered {} scan(s) for field='{}' in {} via msmd: {}".format(len(out), fieldname, msfile, out))
                return out
        except Exception:
            try:
                msmd.close()
            except Exception:
                pass
            # fall through to next method

    # -------------- Attempt 2: ms tool selection + getscansummary --------------
    if '_MS_TOOL_AVAILABLE' in globals() and _MS_TOOL_AVAILABLE:
        try:
            ms = _ms_cls()
            ms.open(msfile)
            # select by *name* is supported by ms.msselect
            ms.msselect({'field': str(fieldname)})
            summ = ms.getscansummary()  # dict keyed by scan number (strings or ints)
            ms.close()
            scans = sorted(int(k) for k in summ.keys())
            if scans:
                print("Discovered {} scan(s) for field='{}' in {} via ms tool: {}".format(len(scans), fieldname, msfile, scans))
                return scans
        except Exception:
            try:
                ms.close()
            except Exception:
                pass
            # fall through

    # -------------- Attempt 3: table fallback (FIELD -> MAIN, chunked) --------------
    if not ('_TB_TOOL_AVAILABLE' in globals() and _TB_TOOL_AVAILABLE):
        print("Table tool unavailable; cannot discover scans.")
        return []

    def _field_ids_for_name_local(msfile_local, fname):
        tb = None
        out_ids = []
        try:
            tb = _tb_cls()
            tb.open(os.path.join(msfile_local, 'FIELD'))
            names = tb.getcol('NAME')
            for idx, nm in enumerate(names):
                try:
                    nm_str = nm.decode() if hasattr(nm, 'decode') else str(nm)
                except Exception:
                    nm_str = str(nm)
                if nm_str == fname:
                    out_ids.append(int(idx))
        except Exception:
            out_ids = []
        finally:
            try:
                tb.close()
            except Exception:
                pass
        return out_ids

    fids = _field_ids_for_name_local(msfile, fieldname)
    if not fids:
        print("No FIELD_ID found with NAME='{}' in {}.".format(fieldname, msfile))
        return []

    fid_set = set(int(x) for x in fids)

    tb = None
    scans_collected = set()
    try:
        tb = _tb_cls()
        tb.open(msfile)  # MAIN
        nrows = int(tb.nrows())
        start = 0

        # Helper that avoids np.isin dependency on older NumPy:
        def _mask_in_set(arr, sset):
            # arr may be numpy or list-like; return a Python list of booleans
            return [int(v) in sset for v in arr]

        while start < nrows:
            nget = min(int(chunk_rows), nrows - start)
            field_ids = tb.getcol('FIELD_ID', startrow=start, nrow=nget)
            scan_nums = tb.getcol('SCAN_NUMBER', startrow=start, nrow=nget)
            if field_ids is None or scan_nums is None:
                break
            m = _mask_in_set(field_ids, fid_set)
            # Collect scans where mask is True
            for flag, sc in zip(m, scan_nums):
                if flag:
                    try:
                        scans_collected.add(int(sc))
                    except Exception:
                        pass
            start += nget
    except Exception as e:
        print("Scan discovery failed for {}: {}".format(msfile, e))
        return []
    finally:
        try:
            tb.close()
        except Exception:
            pass

    scans = sorted(scans_collected)
    if scans:
        print("Discovered {} scan(s) for field='{}' in {} (table fallback): {}".format(len(scans), fieldname, msfile, scans))
    else:
        print("No scans found for field='{}' in {}.".format(fieldname, msfile))
    return scans

# ---------- Selection row-count sanity check (optional) ----------

def _nrows_for_selection(msfile, field="", scan="", spw=""):
    """
    Return number of rows matching a selection, using the 'ms' tool; None if unavailable/error.
    """
    if not _MS_TOOL_AVAILABLE:
        return None
    ms = None
    try:
        ms = _ms_cls()
        ms.open(msfile)
        sel = {}
        if field:
            sel['field'] = str(field)
        if scan != "":
            sel['scan'] = str(scan)
        if spw != "":
            sel['spw'] = str(spw)
        ms.msselect(sel)
        n = int(ms.nrow())
        return n
    except Exception:
        return None
    finally:
        try:
            ms.close()
        except Exception:
            pass

# ---------- Imaging ----------

def run_one(msfile, lo_hz, hi_hz, suffix):
    """Band-level imaging for the full dataset selection."""
    spw_sel = _spw_overlap_selector(msfile, lo_hz, hi_hz)
    imgname = imagename_for(msfile, suffix)
    casa_image = imgname + '.image'
    fits_out   = imgname + '.fits'

    if os.path.exists(casa_image) and os.path.exists(fits_out):
        print("Found existing products for {}. Running QA only.".format(imgname))
        qa_fits(fits_out)
        return

    if os.path.exists(casa_image) and not os.path.exists(fits_out):
        print("Found existing CASA image. Exporting FITS -> {}".format(fits_out))
        exportfits(imagename=casa_image, fitsimage=fits_out, overwrite=True, dropdeg=False, stokeslast=True)
        qa_fits(fits_out)
        return

    # Optional sanity check
    nrows = _nrows_for_selection(msfile, field=FIELDNAME, scan="", spw=spw_sel)
    if nrows == 0:
        print("[band:{}] Warning: selection has 0 rows (spw='{}'). Proceeding anyway to let tclean report."
              .format(suffix, spw_sel))

    try:
        tclean(vis=msfile,
               imagename=imgname,
               field=FIELDNAME,
               datacolumn=DATACOL,
               cell=CELL,
               imsize=IMSIZE,
               niter=NITER,
               stokes=STOKES,
               gain=GAIN,
               deconvolver=DECONV,
               threshold=THRESH,
               gridder=GRIDDER,
               wprojplanes=WPROJ,
               specmode=SPECMODE,
               weighting=WEIGHTING,
               spw=spw_sel,
               robust=ROBUST,
               pblimit=PBLIMIT)
    except Exception as e:
        print("Band tclean failed for {} ({}–{} Hz): {}".format(msfile, lo_hz, hi_hz, e))

    if os.path.exists(casa_image):
        exportfits(imagename=casa_image, fitsimage=fits_out, overwrite=True, dropdeg=False, stokeslast=True)
        qa_fits(fits_out)
    else:
        print("No CASA image produced for band {} — skipping FITS/QA.".format(imgname))

def run_one_scan_fullband(msfile, suffix, scan_id):
    """
    Per-scan imaging for a single scan ID using the ENTIRE band (spw="").
    """
    imgname = imagename_for(msfile, suffix, scan_id=scan_id)
    casa_image = imgname + '.image'
    fits_out   = imgname + '.fits'

    if os.path.exists(casa_image) and os.path.exists(fits_out):
        print("Found existing products for scan {}. QA only.".format(scan_id))
        qa_fits(fits_out)
        return

    if os.path.exists(casa_image) and not os.path.exists(fits_out):
        print("Found existing CASA image for scan {}. Exporting FITS.".format(scan_id))
        exportfits(imagename=casa_image, fitsimage=fits_out, overwrite=True, dropdeg=False, stokeslast=True)
        qa_fits(fits_out)
        return

    # Sanity check: make sure field+scan actually select rows with full band
    nrows = _nrows_for_selection(msfile, field=FIELDNAME, scan=str(int(scan_id)), spw="")
    if nrows is not None and nrows == 0:
        print("[scan:{}:{}] 0 rows for field='{}', scan='{}'. Skipping this scan."
              .format(suffix, scan_id, FIELDNAME, scan_id))
        return

    try:
        tclean(vis=msfile,
               imagename=imgname,
               field=FIELDNAME,
               scan=str(int(scan_id)),
               datacolumn=DATACOL,
               cell=CELL,
               imsize=IMSIZE,
               niter=NITER,
               stokes=STOKES,
               gain=GAIN,
               deconvolver=DECONV,
               threshold=THRESH,
               gridder=GRIDDER,
               wprojplanes=WPROJ,
               specmode=SPECMODE,
               weighting=WEIGHTING,
               spw="",  # ENTIRE BAND for per-scan images
               robust=ROBUST,
               pblimit=PBLIMIT)
    except Exception as e:
        print("Per-scan tclean failed for scan {} (full band): {}".format(scan_id, e))

    if os.path.exists(casa_image):
        exportfits(imagename=casa_image, fitsimage=fits_out, overwrite=True, dropdeg=False, stokeslast=True)
        qa_fits(fits_out)
    else:
        print("No CASA image produced for scan {} — skipping FITS/QA.".format(scan_id))

def run_scans_fullband(msfile, scan_ids):
    """
    Image each discovered scan ID but keep spw='' (full band).
    """
    if not IMAGE_SCANS:
        print("Per-scan imaging disabled by IMAGE_SCANS toggle.")
        return
    if not scan_ids:
        print("Per-scan imaging skipped: no scan IDs discovered for field='{}'.".format(FIELDNAME))
        return
    for sc in scan_ids:
        try:
            sc_int = int(sc)
        except Exception:
            print("Skip invalid scan id '{}'".format(sc))
            continue
        # Use full band but keep separate output name suffixes to mirror your layout
        run_one_scan_fullband(msfile, suffix='msf', scan_id=sc_int)

# ---------- Driver ----------

def resolve_ms_list(argv):
    msfile_override = (MSFILE is not None) and isinstance(MSFILE, str) and MSFILE.strip() != ""
    if msfile_override:
        ms = MSFILE.strip()
        if not os.path.exists(ms):
            print("Error: MSFILE override not found -> {}".format(ms))
            sys.exit(2)
        print("MSFILE override in use. Ignoring CLI args. Target: {}".format(ms))
        return [ms]
    if len(argv) < 1:
        print("Error: provide one or more .ms paths, or set MSFILE.")
        sys.exit(2)
    return argv

def preflight_all(targets):
    """
    If all band products exist for all targets, run QA on those FITS files only and exit.
    (We still proceed to per-scan imaging because scans are auto-discovered now.)
    """
    needed = []
    for ms in targets:
        base_low  = imagename_for(ms, 'lowband')
        base_high = imagename_for(ms, 'highband')
        needed.append((base_low  + '.image', base_low  + '.fits'))
        needed.append((base_high + '.image', base_high + '.fits'))
    all_exist = all(os.path.exists(img) and os.path.exists(fits_) for img, fits_ in needed)
    if all_exist:
        print("All band products exist. Running QA only for bands; continuing to per-scan imaging.")
        for _, fits_ in needed:
            qa_fits(fits_)

def main(argv):
    targets = resolve_ms_list(argv)
    preflight_all(targets)

    for ms in targets:
        if not os.path.isdir(ms):
            print("Skip: not found -> {}".format(ms))
            continue

        # Decide effective bands for this MS (table-based, no msmd) for band-level images
        (low_lo, low_hi), (high_lo, high_hi) = decide_bands(ms)

        # Band-level images (with band SPW selection)
        run_one(ms, low_lo,  low_hi,  'lowband')
        run_one(ms, high_lo, high_hi, 'highband')

        # Auto-discover scan IDs for this field and image each scan with full band
        if IMAGE_SCANS:
            if "1757723806" in ms:
               scan_ids = [3, 5, 7, 9] 
            elif "1757723807" in ms:
                scan_ids = [3,5,7]
            elif "1692135074" in ms:
                scan_ids = [3,5]
            else:
                scan_ids = discover_scans_for_field(ms, FIELDNAME)
        run_scans_fullband(ms, scan_ids)

if __name__ == '__main__':
    # CASA passes script name in sys.argv[0] and additional args afterward
    main(sys.argv[1:])

    if CLEANUPPRODUCTS:
        import subprocess
        cmd = (
            r"find . "
            r"-path '*/*images*/*' "
            r"\( -type d "
            r"\( -name '*.pb' "
            r"-o -name '*.psf' "
            r"-o -name '*.residual' "
            r"-o -name '*.model' "
            r"-o -name '*.mask' "
            r"-o -name '*.sumwt' "
            r"-o -name '*.weight' "
            r"-o -name '*.flux' "
            r"-o -name '*.image.pbcor' "
            r"-o -name '*.tt0' "
            r"-o -name '*.tt1' "
            r"-o -name '*.tt2' \) "
            r"-o -type f "
            r"! -name '*.image' "
            r"! -name '*.fits' \) "
            r"-exec rm -rf {} +"
        )
        subprocess.check_call(cmd, shell=True)
        print("Cleaned up CASA image products, leaving only .image and .fits files.")
