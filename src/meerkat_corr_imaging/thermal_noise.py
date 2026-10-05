"""Natural-weighting Stokes-I radiometer noise for the existing flux analysis.

SEFDs are per antenna, from ESDKB-Sensitivity calculators-051026-112758.pdf,
pages 5–6. No imaging-weight, taper, confusion or calibration correction.
"""
from __future__ import annotations

from contextlib import contextmanager
from pathlib import Path
import math
import re

import numpy as np

SEFD_JY = {"L": 425.0, "UHF": 550.0, "S0": 365.0, "S1": 364.0,
           "S2": 365.0, "S3": 366.0, "S4": 369.0}
REFERENCE = "ESDKB-Sensitivity calculators-051026-112758.pdf, pp. 5–6"
FORMULA = "1e6 * SEFD_Jy / sqrt(2 * effective_bandwidth_hz * on_source_integration_s * N * (N-1))"
ASSUMPTIONS = (
    "Natural weighting; identical band-mean per-antenna SEFD; independent parallel hands; "
    "positive weights are eligibility gates, not sensitivity weights. No Briggs, taper, "
    "confusion, calibration, dynamic-range or primary-beam correction."
)


def radiometer_rms_ujy(antenna_count, effective_bandwidth_hz, on_source_integration_s, sefd_jy):
    """SEFD / sqrt(2 * bandwidth * time * N * (N-1)), returned in µJy/beam."""
    n = float(antenna_count)
    if isinstance(antenna_count, bool) or not math.isfinite(n) or n != int(n) or n < 2:
        raise ValueError("antenna_count must be an integer >= 2")
    bandwidth, seconds, sefd = map(float, (effective_bandwidth_hz, on_source_integration_s, sefd_jy))
    if not all(math.isfinite(v) and v > 0 for v in (bandwidth, seconds, sefd)):
        raise ValueError("bandwidth, integration and SEFD must be finite and positive")
    return 1e6 * sefd / math.sqrt(2 * bandwidth * seconds * n * (n - 1))


def band_from_label(label):
    """Only explicit S0–S4 labels identify overlapping S sub-bands reliably."""
    match = re.search(r"S[0-4](?!\d)", str(label).upper())
    return match.group() if match else None


def _band(value, frequencies=()):
    if value is not None:
        value = str(value).upper().replace("-BAND", "").strip()
        if value not in SEFD_JY:
            raise ValueError("band must be L, UHF or S0–S4; generic S is ambiguous")
        return value
    frequency = np.asarray(frequencies, dtype=float)
    candidates = [b for b, lo, hi in (("UHF", .544e9, 1.088e9), ("L", .856e9, 1.712e9))
                  if len(frequency) and frequency.min() >= lo-1 and frequency.max() <= hi+1]
    if len(candidates) != 1:
        raise ValueError("Observing band is ambiguous; specify band (especially S0–S4)")
    return candidates[0]


@contextmanager
def _table(factory, path):
    tb = factory(str(path), readonly=True, ack=False)
    try:
        yield tb
    finally:
        tb.close()


def _union_length(intervals):
    end, length = -math.inf, 0.0
    for lo, hi in sorted(intervals):
        length += max(0.0, hi - max(lo, end))
        end = max(end, hi)
    return length


def _cell_array(rows, column, start, count, shape, *, default=np.nan):
    """Read regular chunks efficiently; undefined cells remain unusable."""
    try:
        data = np.asarray(rows.getcol(column, startrow=start, nrow=count))
        if data.shape == (count,) + shape:
            return data
    except (RuntimeError, KeyError):
        pass
    data = np.full((count,) + shape, default)
    for i in range(count):
        if rows.iscelldefined(column, start+i):
            cell = np.asarray(rows.getcell(column, start+i))
            if cell.shape != shape:
                raise ValueError("Malformed " + column + " cell")
            data[i] = cell
    return data


def ms_sensitivity(spec, table_factory=None):
    """Derive effective bandwidth from usable baseline/parallel-hand exposure.

    With C = sum(valid_cell * abs(CHAN_WIDTH) * EXPOSURE), define
    bandwidth_eff = C / [N(N-1) * t_on_source]. This reduces to the actual
    unflagged bandwidth for complete dual-hand data, and accounts for partial
    flagging, missing baselines and varying antenna participation without
    treating the union of intermittently usable channels as fully exposed.
    """
    if table_factory is None:
        from casacore.tables import table as table_factory
    path = Path(spec["ms"]).expanduser().resolve()
    if not path.is_dir():
        raise FileNotFoundError(path)
    with _table(table_factory, path / "FIELD") as tb:
        names = list(map(str, tb.getcol("NAME")))
        field = spec.get("field")
        if field is None:
            raise ValueError("field is required; do not include calibrator time in target sensitivity")
        ids = [int(field)] if str(field).isdigit() else [i for i, name in enumerate(names) if name == field]
        if len(ids) != 1 or not 0 <= ids[0] < len(names):
            raise ValueError("field must identify exactly one MS field")
        field_id = ids[0]
    windows, full_frequencies, intervals = {}, [], []
    with _table(table_factory, path / "SPECTRAL_WINDOW") as tb:
        for spw in range(tb.nrows()):
            freq = np.asarray(tb.getcell("CHAN_FREQ", spw), dtype=float)
            widths = np.abs(np.asarray(tb.getcell("CHAN_WIDTH", spw), dtype=float))
            if freq.ndim != 1 or widths.shape != freq.shape or not np.isfinite(freq).all() or not np.all(np.isfinite(widths) & (widths > 0)):
                raise ValueError("Invalid channel frequencies/widths")
            full_frequencies.extend(freq)
            selected = np.ones(len(freq), dtype=bool)
            if spec.get("channels") not in (None, "all"):
                selected[:] = False
                indices = spec["channels"].get(str(spw), spec["channels"].get(spw, []))
                if any(isinstance(i, bool) or int(i) != i or not 0 <= int(i) < len(freq) for i in indices):
                    raise ValueError("Invalid channel index")
                selected[list(map(int, indices))] = True
            if spec.get("frequency_range_hz") is not None:
                lo, hi = map(float, spec["frequency_range_hz"])
                if not math.isfinite(lo+hi) or not lo < hi:
                    raise ValueError("frequency_range_hz must be finite and increasing")
                selected &= (freq >= lo) & (freq <= hi)
            windows[spw] = (freq, widths, selected)
    band = _band(spec.get("band"), full_frequencies)
    with _table(table_factory, path / "POLARIZATION") as tb:
        polarizations = [list(map(int, tb.getcell("CORR_TYPE", i))) for i in range(tb.nrows())]
    with _table(table_factory, path / "DATA_DESCRIPTION") as tb:
        ddids = list(zip(tb.getcol("SPECTRAL_WINDOW_ID"), tb.getcol("POLARIZATION_ID")))
    scans = spec.get("scans")
    if scans == "all":
        scans = None
    if isinstance(scans, str):
        # CASA's richer scan syntax is intentionally not guessed.
        scans = [int(s) for s in scans.split(",") if s.strip()]
    scans = set(map(int, scans)) if scans else None
    datacolumn = {"corrected": "CORRECTED_DATA", "data": "DATA"}.get(
        str(spec.get("datacolumn", "data")).lower(), str(spec.get("datacolumn", "DATA")).upper())
    participants, dumps, seen, active_spws, observed_scans = set(), {}, set(), set(), set()
    cell_exposure = flagged_exposure = selected_exposure = 0.0
    valid_channels = {}
    with _table(table_factory, path) as main:
        for ddid, (spw, pol) in enumerate(ddids):
            freq, widths, selected = windows[int(spw)]
            if not selected.any():
                continue
            query = f"FIELD_ID == {field_id} && DATA_DESC_ID == {ddid} && ANTENNA1 != ANTENNA2"
            if scans:
                query += " && SCAN_NUMBER IN [" + ",".join(map(str, sorted(scans))) + "]"
            rows = main.query(query)
            try:
                if not rows.nrows():
                    continue
                active_spws.add(int(spw))
                corr = polarizations[int(pol)]
                pair = [9, 12] if 9 in corr and 12 in corr else [5, 8]
                if not all(c in corr for c in pair):
                    raise ValueError("Stokes I requires both XX/YY or RR/LL")
                hands = [corr.index(c) for c in pair]
                chunk_size = max(1, min(1024, 32_000_000 // (len(freq)*len(corr)*40)))
                support = valid_channels.setdefault(int(spw), np.zeros(len(freq), dtype=bool))
                for start in range(0, rows.nrows(), chunk_size):
                    count = min(chunk_size, rows.nrows()-start)
                    scalar = {c: np.asarray(rows.getcol(c, startrow=start, nrow=count))
                              for c in ("ANTENNA1", "ANTENNA2", "TIME", "INTERVAL", "EXPOSURE", "FLAG_ROW", "SCAN_NUMBER")}
                    observed_scans.update(map(int, scalar["SCAN_NUMBER"]))
                    exposure = scalar["EXPOSURE"]
                    if not np.all(np.isfinite(exposure) & (exposure > 0) & (exposure <= scalar["INTERVAL"] * (1+1e-8))) or not np.all(np.isfinite(scalar["INTERVAL"]) & (scalar["INTERVAL"] > 0)) or not np.isfinite(scalar["TIME"]).all():
                        raise ValueError("Invalid MS exposure, interval or time")
                    shape = (len(freq), len(corr))
                    flags = _cell_array(rows, "FLAG", start, count, shape, default=True).astype(bool)[:, :, hands]
                    flags |= scalar["FLAG_ROW"][:, None, None]
                    data = _cell_array(rows, datacolumn, start, count, shape, default=complex(np.nan, np.nan))[:, :, hands]
                    weight = _cell_array(rows, "WEIGHT", start, count, (len(corr),))[:, None, hands]
                    weight = np.broadcast_to(weight, flags.shape).copy()
                    if "WEIGHT_SPECTRUM" in rows.colnames():
                        # Undefined spectrum falls back to WEIGHT; defined invalid spectrum never does.
                        spectrum = _cell_array(rows, "WEIGHT_SPECTRUM", start, count, shape)[:, :, hands]
                        for i in range(count):
                            if rows.iscelldefined("WEIGHT_SPECTRUM", start+i):
                                weight[i] = spectrum[i]
                    good = ~flags & np.isfinite(data) & np.isfinite(weight) & (weight > 0) & selected[None, :, None]
                    cell_exposure += float(np.sum(good * widths[None, :, None] * exposure[:, None, None]))
                    selected_exposure += float(2 * widths[selected].sum() * exposure.sum())
                    flagged_exposure += float(np.sum(flags[:, selected] * widths[None, selected, None] * exposure[:, None, None]))
                    support |= good.any(axis=(0, 2))
                    for i in range(count):
                        a, b = int(scalar["ANTENNA1"][i]), int(scalar["ANTENNA2"][i])
                        time, dt = float(scalar["TIME"][i]), float(scalar["INTERVAL"][i])
                        key = (int(spw), min(a,b), max(a,b), time)
                        if key in seen:
                            raise ValueError("Duplicate baseline/SPW/time samples; exposure would be double-counted")
                        seen.add(key)
                        dumps[(time, dt)] = (time-dt/2, time+dt/2)
                        if good[i].any():
                            participants.update((a, b))
            finally:
                rows.close()
    for spw in active_spws:
        freq, widths, selected = windows[spw]
        intervals.extend((float(f-w/2), float(f+w/2)) for f,w in zip(freq[selected], widths[selected]))
    nominal = sum(windows[s][1][windows[s][2]].sum() for s in active_spws)
    if nominal and not math.isclose(_union_length(intervals), nominal, rel_tol=1e-8):
        raise ValueError("Overlapping selected spectral channels; independent bandwidth cannot be assumed")
    n, seconds = len(participants), _union_length(dumps.values())
    if n < 2 or seconds <= 0 or cell_exposure <= 0:
        raise ValueError("No usable on-source dual-polarization cross-correlation data")
    bandwidth = cell_exposure / (n*(n-1)*seconds)
    return dict(antenna_count=n, antenna_ids=sorted(participants), effective_bandwidth_hz=bandwidth,
                on_source_integration_s=seconds, band=band, ms=str(path), field=names[field_id],
                field_id=field_id, scans=sorted(observed_scans), datacolumn=datacolumn,
                nominal_selected_bandwidth_hz=float(nominal),
                usable_channel_bandwidth_hz=float(sum(windows[s][1][mask].sum() for s,mask in valid_channels.items())),
                channel_selection={str(s): np.flatnonzero(windows[s][2]).tolist() for s in active_spws},
                flag_fraction=flagged_exposure/selected_exposure,
                usable_fraction_of_present_samples=cell_exposure/selected_exposure,
                parallel_hand_exposure_hz_s=cell_exposure,
                input_method="MS flags, finite data and positive weights; exposure-equivalent bandwidth",
                approximation="A single effective bandwidth represents variable flagging/array participation. "
                    "t_int is the union of selected on-source dump intervals (not baseline-summed or elapsed scan span). "
                    "MS EXPOSURE accounts for integration loss within each dump.")


def sensitivity(spec, table_factory=None):
    """Never invent sensitivity from source catalogues or nominal bandwidth."""
    if spec.get("manifest"):
        import json
        manifest = json.loads(Path(spec["manifest"]).read_text())
        if manifest.get("processing_status") != "complete":
            raise ValueError("Imaging manifest must be complete")
        contract = manifest["contract"]
        params = contract["tclean_parameters"]
        spec = {**spec, "ms": params["vis"], "field": contract["field_id"],
                "scans": contract["scan_ids"], "datacolumn": params["datacolumn"],
                "channels": {}}
        spec.pop("frequency_range_hz", None)
        for channel in contract["requested_channels"]:
            spec["channels"].setdefault(str(channel["spw"]), []).append(channel["channel"])
    if spec.get("ms"):
        result = ms_sensitivity(spec, table_factory)
    elif spec.get("visibility_results"):
        if spec.get("allow_approximate") is not True:
            raise ValueError("Cached visibility approximation requires allow_approximate: true")
        result = summary_sensitivity(spec)
    else:
        required = ("antenna_count", "effective_bandwidth_hz", "on_source_integration_s", "band", "provenance")
        if any(k not in spec for k in required) or not str(spec["provenance"]).strip():
            raise ValueError("MS/complete imaging manifest unavailable; explicit N, usable bandwidth, on-source time, band and provenance required")
        result = {k: spec[k] for k in required}
        result.update(band=_band(spec["band"]), input_method="explicit metadata",
                      approximation=spec.get("approximation", "User-supplied usable bandwidth and time; not independently verified from MS"))
    result.update(sefd_jy=SEFD_JY[result["band"]], sefd_reference=REFERENCE, formula=FORMULA,
                  sefd_assumption="Band-mean SEFD; S1 inferred from overlapping sub-bands in the reference",
                  assumptions=ASSUMPTIONS)
    result["theoretical_rms_ujy_beam"] = radiometer_rms_ujy(
        result["antenna_count"], result["effective_bandwidth_hz"], result["on_source_integration_s"], result["sefd_jy"])
    result["status"] = "available"
    return result


def summary_sensitivity(spec):
    """Approximation from existing visibility QA CSVs when the MS is unavailable.

    Channel-width/dump-duration metadata and weight validity were not retained
    by that analyser. Infer uniform widths/cadence, preserve actual flag losses,
    and clearly label these assumptions rather than using nominal bandwidth.
    """
    import pandas as pd
    path = Path(spec['visibility_results']).expanduser().resolve()
    rows = pd.read_csv(path / 'perrow_amp_stats.csv', usecols=[
        'TIME', 'ANT1', 'ANT2', 'CLASS', 'SCAN', 'SPW', 'POL', 'FLAG_FRAC'])
    channels = pd.read_csv(path / 'rfi_free_channel_mask.csv')
    rows = rows[rows.CLASS == 'cross']
    channels = channels[channels.CLASS == 'cross']
    if spec.get('scans') or spec.get('channels'):
        raise ValueError('Cached channel flag fractions cannot recover scan/channel-index selections; use MS or explicit metadata')
    hands = set(rows.POL.unique())
    pair = ['XX', 'YY'] if {'XX', 'YY'} <= hands else ['RR', 'LL']
    if not set(pair) <= hands:
        raise ValueError('Cached Stokes I needs both parallel hands')
    rows = rows[rows.POL.isin(pair)]
    channels = channels[channels.POL.isin(pair)]
    eligible = rows[np.isfinite(rows.FLAG_FRAC) & (rows.FLAG_FRAC < 1)]
    antennas = sorted(set(map(int, eligible.ANT1)) | set(map(int, eligible.ANT2)))
    n = len(antennas)
    # No elapsed calibrator gaps; infer duration from within-scan dump spacings.
    steps = np.concatenate([np.diff(np.sort(group.TIME.unique()))
                            for _, group in rows.groupby('SCAN')]) if len(rows) else np.array([])
    steps = steps[np.isfinite(steps) & (steps > 0)]
    if n < 2 or not len(steps):
        raise ValueError('Cached rows lack usable antennas or within-scan dump cadence')
    dump = float(np.median(steps))
    seconds = len(rows.TIME.unique()) * dump
    band = _band(spec.get('band'), channels.FREQ_HZ.unique())
    exposure, flagged, present = 0.0, 0.0, 0.0
    widths = {}
    selected_channels = {}
    for spw, group in channels.groupby('SPW'):
        grid = group[['CHAN', 'FREQ_HZ']].drop_duplicates().sort_values('CHAN')
        if grid.CHAN.duplicated().any() or len(grid) < 2:
            raise ValueError('Cached frequency grid does not identify channel widths')
        spacing = np.abs(np.diff(grid.FREQ_HZ) / np.diff(grid.CHAN))
        width = float(np.median(spacing))
        if not width > 0 or not np.allclose(spacing, width, rtol=1e-6):
            raise ValueError('Nonuniform cached channel widths cannot be recovered')
        widths[str(int(spw))] = width
        if spec.get('frequency_range_hz'):
            lo, hi = map(float, spec['frequency_range_hz'])
            if not math.isfinite(lo+hi) or not lo < hi:
                raise ValueError('frequency_range_hz must be finite and increasing')
            group = group[(group.FREQ_HZ >= lo) & (group.FREQ_HZ <= hi)]
        selected_channels[str(int(spw))] = sorted(map(int, group.CHAN.unique()))
        for hand in pair:
            selected_hand = group[group.POL == hand]
            if selected_hand.CHAN.duplicated().any() or set(selected_hand.CHAN) != set(group.CHAN):
                raise ValueError('Cached channels must have one entry for each parallel hand')
            fractions = selected_hand.FLAG_FRAC.to_numpy(float)
            if not np.all(np.isfinite(fractions) & (fractions >= 0) & (fractions <= 1)):
                raise ValueError('Invalid cached channel flag fractions')
            count = len(rows[(rows.SPW == spw) & (rows.POL == hand)])
            exposure += float((1-fractions).sum()*count*width*dump)
            flagged += float(fractions.sum()*count*width*dump)
            present += float(len(fractions)*count*width*dump)
    bandwidth = exposure/(n*(n-1)*seconds)
    if bandwidth <= 0:
        raise ValueError('No usable cached bandwidth in selected frequency range')
    return dict(antenna_count=n, antenna_ids=antennas, effective_bandwidth_hz=bandwidth,
                on_source_integration_s=seconds, band=band, visibility_results=str(path),
                channel_widths_hz=widths, inferred_dump_seconds=dump, channel_selection=selected_channels,
                flag_fraction=flagged/present, parallel_hand_exposure_hz_s=exposure,
                input_method='approximate from existing visibility QA CSVs',
                approximation='Uniform channel width inferred from CHAN/FREQ_HZ; uniform dump duration '
                    'inferred from median within-scan cadence; EXPOSURE/INTERVAL and positive-weight/finite-data '
                    'validity unavailable. Cached input must describe this field and all selected scans. '
                    'RFI_FREE spectral-diagnostic mask is not an imaging mask and is not applied. '
                    'Extra imaging cuts require an MS/complete image manifest or explicit metadata.')


def thermal_noise_results(config=None):
    """Per-product/per-catalogue-side diagnostics, preserving catalogue-only flux runs."""
    config = config or {}
    unknown = set(config) - {"low", "high", "mfs"}
    if unknown:
        raise ValueError("Unknown thermal_noise products: " + ", ".join(sorted(unknown)))
    results = {}
    for product in ("low", "high", "mfs"):
        results[product] = {}
        for role in ("reference", "test"):
            spec = config.get(product, {}).get(role, {})
            try:
                results[product][role] = sensitivity(spec)
            except (ValueError, OSError, ImportError, RuntimeError, KeyError) as exc:
                results[product][role] = dict(status="unavailable", theoretical_rms_ujy_beam=None,
                                             reason=str(exc), inputs=spec, assumptions=ASSUMPTIONS)
    return results
