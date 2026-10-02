"""CASA 6 paired frequency planning and provenance; no CASA import at import time.

MAIN is streamed in bounded, fixed-DDID chunks. Inputs are opened read-only.
The planner retains whole native channels; overlap never implies resampling.
"""
import copy
import bisect
import datetime
import glob
import hashlib
import inspect
import json
import math
import os
import re
from contextlib import contextmanager

import numpy as np

SCHEMA = 1
FREQ_FRAMES = ['REST', 'LSRK', 'LSRD', 'BARY', 'GEO', 'TOPO', 'GALACTO', 'LGROUP', 'CMB']
CORRELATIONS = {5: 'RR', 8: 'LL', 9: 'XX', 12: 'YY'}
POLICY = {
    'denominator': 'present cross-correlation rows x selected parallel hands, per channel; including flagged/zero-weight/absent cells',
    'flag_bit_policy': 'CASA boolean FLAG OR FLAG_ROW; all stored correlations must be unflagged for strict Stokes I; raw BITFLAG/BITFLAG_ROW are reported but not interpreted or applied by tclean',
    'weights': 'WEIGHT_SPECTRUM when defined and nonempty, otherwise broadcast WEIGHT; no fallback for zero/negative/nonfinite spectrum weights',
    'absent_data': 'undefined FLAG/DATA/weight cells count as unusable in present-row denominator; SPWs without rows have denominator zero and are ineligible; unscheduled rows are not invented',
    'stokes_i': 'both XX/YY or RR/LL must have finite DATA and positive finite weights; jointly valid hands contribute together',
    'effective_frequency': 'sum of native channel centres weighted by summed positive input weights on jointly valid Stokes-I opportunities; not Briggs/gridding weights',
}
DEFAULTS = dict(enabled=False, occupancy_threshold=0.8, min_common_coverage=0.9,
                min_channel_overlap=0.9, chunk_rows=128, max_chunk_bytes=67108864,
                datacolumn='data', dry_run=False, output_dir='', fields=[], tclean={})
IMAGING_DEFAULTS = dict(cell='1.2arcsec', imsize=4096, niter=10000, gain=0.1,
                        deconvolver='clark', threshold='1uJy', gridder='wproject',
                        wprojplanes=64, weighting='briggs', robust=0.0, pblimit=-1.0)
# Selections and side effects are controlled by the experiment, not task overrides.
RESERVED = {'vis', 'imagename', 'field', 'scan', 'spw', 'datacolumn', 'stokes',
            'specmode', 'antenna', 'selectdata', 'restart', 'savemodel', 'parallel',
            'timerange', 'uvrange', 'observation', 'intent', 'outlierfile', 'startmodel',
            'phasecenter', 'interactive', 'restoration', 'calcpsf', 'calcres'}


class SelectionPlanningError(ValueError):
    def __init__(self, message, diagnostics):
        super().__init__(message)
        self.diagnostics = diagnostics


def validate_config(raw):
    if not isinstance(raw, dict):
        raise ValueError('paired_astrometry must be a mapping')
    unknown = set(raw) - set(DEFAULTS)
    if unknown:
        raise ValueError('Unknown paired_astrometry options: ' + ', '.join(sorted(unknown)))
    cfg = copy.deepcopy(DEFAULTS)
    cfg.update(copy.deepcopy(raw))
    for name in ('enabled', 'dry_run'):
        if not isinstance(cfg[name], bool):
            raise ValueError(name + ' must be boolean')
    for name in ('occupancy_threshold', 'min_common_coverage', 'min_channel_overlap'):
        value = float(cfg[name])
        if not math.isfinite(value) or not 0 < value <= 1:
            raise ValueError(name + ' must be in (0, 1]')
        cfg[name] = value
    for name in ('chunk_rows', 'max_chunk_bytes'):
        if isinstance(cfg[name], bool) or int(cfg[name]) != cfg[name] or cfg[name] <= 0:
            raise ValueError(name + ' must be a positive integer')
        cfg[name] = int(cfg[name])
    if str(cfg['datacolumn']).lower() not in ('data', 'corrected'):
        raise ValueError('datacolumn must be data or corrected')
    cfg['datacolumn'] = str(cfg['datacolumn']).lower()
    if not isinstance(cfg['tclean'], dict) or RESERVED.intersection(cfg['tclean']):
        raise ValueError('tclean overrides cannot change controlled selections/side effects: ' + ', '.join(sorted(RESERVED)))
    if cfg['tclean'].get('deconvolver', 'clark') == 'mtmfs':
        raise ValueError('Use single-term MFS deconvolution (e.g. clark); mtmfs is not supported by this product contract')
    if cfg['enabled']:
        if not cfg['output_dir'] or not isinstance(cfg['fields'], list) or not cfg['fields']:
            raise ValueError('paired_astrometry needs output_dir and fields')
        kinds = []
        names = []
        for field in cfg['fields']:
            if not isinstance(field, dict) or set(field) - {'name', 'kind', 'reference_ms', 'test_ms'}:
                raise ValueError('Each field needs name, kind, reference_ms, test_ms')
            if any(not isinstance(field.get(key), str) or not field[key].strip()
                   for key in ('name', 'kind', 'reference_ms', 'test_ms')):
                raise ValueError('Each field needs nonempty name, kind, reference_ms, test_ms')
            if field['kind'] not in ('gain_calibrator', 'target'):
                raise ValueError('field kind must be gain_calibrator or target')
            if os.path.realpath(field['reference_ms']) == os.path.realpath(field['test_ms']):
                raise ValueError('Reference and test must be distinct MSs')
            kinds.append(field['kind'])
            names.append(field['name'])
        if sorted(kinds) != ['gain_calibrator', 'target'] or len(set(names)) != len(names):
            raise ValueError('Specify one distinct gain_calibrator and one target field')
    return cfg


def write_json(path, value):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    temp = path + '.tmp'
    with open(temp, 'w') as handle:
        json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write('\n')
    os.replace(temp, path)


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


@contextmanager
def table_open(factory, path):
    tb = factory()
    try:
        tb.open(path, nomodify=True)
        yield tb
    finally:
        tb.close()


def ms_identity(path):
    """File-stat identity supplements a digest of selected visibility content.

    Locks are excluded because read-only CASA opens may change them. No MS files
    are hashed in one allocation. Never put the output directory inside the MS.
    """
    path = os.path.realpath(os.path.expanduser(path))
    if not os.path.isdir(path):
        raise IOError('MeasurementSet not found: ' + path)
    entries = []
    for root, dirs, files in os.walk(path):
        dirs.sort()
        for filename in sorted(files):
            if filename == 'table.lock':
                continue
            full = os.path.join(root, filename)
            stat = os.stat(full)
            entries.append([os.path.relpath(full, path), stat.st_size, stat.st_mtime_ns])
    return dict(path=path, file_stat_sha256=digest(entries), file_count=len(entries))


def metadata(path, fieldname, factory):
    with table_open(factory, os.path.join(path, 'FIELD')) as tb:
        names = [str(x) for x in tb.getcol('NAME')]
        ids = [i for i, name in enumerate(names) if name == fieldname]
        if len(ids) != 1:
            raise ValueError('{}: expected exactly one FIELD NAME={!r}; found {}'.format(path, fieldname, ids))
        fid = ids[0]
        if int(tb.getcell('NUM_POLY', fid)) != 0:
            raise ValueError('Moving/polynomial phase centres are unsupported')
        if 'EPHEMERIS_ID' in tb.colnames() and int(tb.getcell('EPHEMERIS_ID', fid)) >= 0:
            raise ValueError('Ephemeris phase centres are unsupported')
        direction = np.asarray(tb.getcell('PHASE_DIR', fid), dtype=float).reshape(2, -1)[:, 0]
        info = tb.getcolkeyword('PHASE_DIR', 'MEASINFO')
        frame = info.get('Ref')
        if 'VarRefCol' in info:
            code = int(tb.getcell(info['VarRefCol'], fid))
            codes = list(info['TabRefCodes'])
            frame = list(info['TabRefTypes'])[codes.index(code)]
        if not frame or not np.isfinite(direction).all():
            raise ValueError('Missing/nonfinite PHASE_DIR reference metadata')
        phase = dict(ra_rad=float(direction[0]), dec_rad=float(direction[1]), frame=str(frame))
    channels = {}
    with table_open(factory, os.path.join(path, 'SPECTRAL_WINDOW')) as tb:
        for spw in range(tb.nrows()):
            freqs = np.asarray(tb.getcell('CHAN_FREQ', spw), dtype=float)
            widths = np.asarray(tb.getcell('CHAN_WIDTH', spw), dtype=float)
            if freqs.ndim != 1 or freqs.shape != widths.shape:
                raise ValueError('Invalid SPW frequency/width dimensions')
            code = int(tb.getcell('MEAS_FREQ_REF', spw))
            if not 0 <= code < len(FREQ_FRAMES):
                raise ValueError('Unknown frequency frame code: ' + str(code))
            channels[spw] = []
            for chan, (freq, width) in enumerate(zip(freqs, widths)):
                valid = bool(np.isfinite(freq) and np.isfinite(width) and abs(width) > 0)
                channels[spw].append(dict(spw=spw, channel=chan,
                    centre_hz=float(freq) if np.isfinite(freq) else None,
                    width_hz=float(abs(width)) if np.isfinite(width) else None,
                    interval_hz=[float(freq-abs(width)/2), float(freq+abs(width)/2)] if valid else None,
                    frequency_frame=FREQ_FRAMES[code], valid_frequency=valid))
    with table_open(factory, os.path.join(path, 'POLARIZATION')) as tb:
        pols = [list(map(int, np.asarray(tb.getcell('CORR_TYPE', i)).ravel())) for i in range(tb.nrows())]
    with table_open(factory, os.path.join(path, 'DATA_DESCRIPTION')) as tb:
        ddids = [(int(tb.getcell('SPECTRAL_WINDOW_ID', i)), int(tb.getcell('POLARIZATION_ID', i)))
                 for i in range(tb.nrows())]
    return dict(field_id=fid, field_name=fieldname, phase_centre=phase,
                channels=channels, pols=pols, ddids=ddids)


STAT_KEYS = ('denominator', 'flagged', 'any_correlation_flagged', 'positive_finite_weights', 'zero_weights',
             'negative_weights', 'nonfinite_weights', 'nonfinite_data', 'absent_samples',
             'absent_flag_samples', 'absent_weight_samples', 'joint_usable', 'weight_sum')


def empty_stats(nchan):
    return {key: np.zeros(nchan, dtype=np.float64 if key == 'weight_sum' else np.int64) for key in STAT_KEYS}


def eligibility(flags, flag_rows, weights, data, selected, absent_rows=None):
    """Arrays are CASA (correlation, channel, row). Strict I opportunities."""
    flags = np.asarray(flags, dtype=bool)
    weights = np.asarray(weights, dtype=float)
    data = np.asarray(data)
    if flags.shape != weights.shape or flags.shape != data.shape or flags.ndim != 3:
        raise ValueError('FLAG, weight and DATA shapes must agree (corr, channel, row)')
    bad = flags | np.asarray(flag_rows, dtype=bool)[None, None, :]
    w = weights[selected]
    d = data[selected]
    finite_w = np.isfinite(w)
    positive = finite_w & (w > 0)
    finite_d = np.isfinite(d)
    joint = ~bad.any(axis=0) & positive.all(axis=0) & finite_d.all(axis=0)
    if absent_rows is not None:
        joint &= ~np.asarray(absent_rows, dtype=bool)[None, :]
    ncorr = len(selected)
    stats = empty_stats(flags.shape[1])
    stats['denominator'][:] = flags.shape[2] * ncorr
    stats['flagged'] = bad[selected].sum(axis=(0, 2))
    stats['any_correlation_flagged'] = bad.any(axis=0).sum(axis=1) * ncorr
    stats['positive_finite_weights'] = positive.sum(axis=(0, 2))
    stats['zero_weights'] = (finite_w & (w == 0)).sum(axis=(0, 2))
    stats['negative_weights'] = (finite_w & (w < 0)).sum(axis=(0, 2))
    stats['nonfinite_weights'] = (~finite_w).sum(axis=(0, 2))
    stats['nonfinite_data'] = (~finite_d).sum(axis=(0, 2))
    stats['joint_usable'] = joint.sum(axis=1) * ncorr
    stats['weight_sum'] = np.where(joint[None, :, :], w, 0).sum(axis=(0, 2))
    return stats


def add_stats(dest, source):
    for key in STAT_KEYS:
        dest[key] += source[key]


def _hash_array(parts, label, arr):
    """Column/row-order hashing makes provenance independent of chunk size."""
    arr = np.asarray(arr)
    kind = arr.dtype.kind
    dtype = {'f': '<f8', 'c': '<c16', 'i': '<i8', 'u': '<u8', 'b': '?'}[kind]
    arr = arr.astype(dtype, copy=False)
    if label not in parts:
        parts[label] = hashlib.sha256()
        parts[label].update(str((arr.dtype.str, arr.shape[:-1])).encode())
    parts[label].update(np.ascontiguousarray(np.moveaxis(arr, -1, 0)).tobytes())


def _cell_defined(tb, column, row):
    return column in tb.colnames() and tb.iscelldefined(column, row)


def _read_chunk(tb, start, count, ncorr, nchan, datacol):
    """Fast getcol path; undefined spectrum cells fall back individually.

    Missing required cells yield an explicit absence mask. A defined malformed
    cell fails rather than silently using a different weight convention.
    """
    shape = (ncorr, nchan, count)
    try:
        flags = np.asarray(tb.getcol('FLAG', startrow=start, nrow=count), dtype=bool)
        data = np.asarray(tb.getcol(datacol, startrow=start, nrow=count))
        row_weights = np.asarray(tb.getcol('WEIGHT', startrow=start, nrow=count), dtype=float)
        if flags.shape != shape or data.shape != shape or row_weights.shape != (ncorr, count):
            raise ValueError('Unexpected fixed-DDID cell shape')
        weights = np.broadcast_to(row_weights[:, None, :], shape).copy()
        sources = ['WEIGHT'] * count
        absent = np.zeros(count, dtype=bool)
        missing_flags = np.zeros(count, dtype=bool)
        missing_weights = np.zeros(count, dtype=bool)
    except RuntimeError:
        flags = np.zeros(shape, dtype=bool)
        data = np.full(shape, np.nan, dtype=complex)
        weights = np.full(shape, np.nan, dtype=float)
        sources = ['absent'] * count
        absent = np.ones(count, dtype=bool)
        missing_flags = np.ones(count, dtype=bool)
        missing_weights = np.ones(count, dtype=bool)
        for i in range(count):
            row = start + i
            for col in ('FLAG', datacol, 'WEIGHT'):
                if _cell_defined(tb, col, row):
                    value = np.asarray(tb.getcell(col, row))
                    expected = (ncorr,) if col == 'WEIGHT' else (ncorr, nchan)
                    if value.shape != expected:
                        raise ValueError('Malformed required cell: ' + col)
                    if col == 'FLAG':
                        flags[:, :, i], missing_flags[i] = value, False
                    elif col == 'WEIGHT':
                        weights[:, :, i], missing_weights[i] = value[:, None], False
                    else:
                        data[:, :, i] = value
            if all(_cell_defined(tb, c, row) for c in ('FLAG', datacol, 'WEIGHT')):
                sources[i], absent[i] = 'WEIGHT', False
    if 'WEIGHT_SPECTRUM' in tb.colnames():
        try:
            spectrum = np.asarray(tb.getcol('WEIGHT_SPECTRUM', startrow=start, nrow=count), dtype=float)
        except RuntimeError:
            spectrum = None
        if spectrum is not None and spectrum.size:
            if spectrum.shape != shape:
                raise ValueError('Malformed WEIGHT_SPECTRUM')
            weights = spectrum
            missing_weights[:] = False
            sources = ['absent' if missing else 'WEIGHT_SPECTRUM' for missing in absent]
        else:
            for i in range(count):
                row = start + i
                if _cell_defined(tb, 'WEIGHT_SPECTRUM', row):
                    w = np.asarray(tb.getcell('WEIGHT_SPECTRUM', row), dtype=float)
                    if w.size:
                        if w.shape != (ncorr, nchan):
                            raise ValueError('Malformed WEIGHT_SPECTRUM')
                        weights[:, :, i], missing_weights[i] = w, False
                        sources[i] = 'WEIGHT_SPECTRUM' if not absent[i] else 'absent'
    return flags, data, weights, absent, sources, missing_flags, missing_weights


def survey_ms(path, fieldname, cfg, factory):
    path = os.path.realpath(os.path.expanduser(path))
    identity = ms_identity(path)
    meta = metadata(path, fieldname, factory)
    aggregate = {spw: empty_stats(len(chans)) for spw, chans in meta['channels'].items()}
    by_scan = {}
    hasher = hashlib.sha256()
    hasher.update(digest(meta).encode())
    hash_parts = {}
    weights_used = {'WEIGHT': 0, 'WEIGHT_SPECTRUM': 0, 'absent': 0}
    raw_bits = {}
    correlations = []
    datacol = 'DATA' if cfg['datacolumn'] == 'data' else 'CORRECTED_DATA'
    with table_open(factory, path) as main:
        if datacol not in main.colnames():
            raise ValueError(path + ': missing ' + datacol)
        for ddid, (spw, pol) in enumerate(meta['ddids']):
            corr = meta['pols'][pol]
            pair = [9, 12] if 9 in corr and 12 in corr else [5, 8]
            if not all(c in corr for c in pair):
                # Only reject unsupported polarization if this DDID has selected rows.
                selected = None
            else:
                selected = [corr.index(c) for c in pair]
            rows = main.query('FIELD_ID == {} && DATA_DESC_ID == {} && ANTENNA1 != ANTENNA2'.format(meta['field_id'], ddid))
            try:
                total = rows.nrows()
                if total == 0:
                    continue
                if selected is None:
                    raise ValueError('Stokes I needs both XX/YY or RR/LL in DDID ' + str(ddid))
                correlations.append(dict(ddid=ddid, spw=spw, stored=corr,
                    selected=[CORRELATIONS[c] for c in pair], indices=selected))
                ncorr, nchan = len(corr), len(meta['channels'][spw])
                # Account for copies of DATA, weights, flags and temporary masks.
                bytes_per_row = max(1, ncorr * nchan * 96)
                count_max = max(1, min(cfg['chunk_rows'], cfg['max_chunk_bytes'] // bytes_per_row))
                for start in range(0, total, count_max):
                    count = min(count_max, total-start)
                    scalars = {c: np.asarray(rows.getcol(c, startrow=start, nrow=count))
                               for c in ('SCAN_NUMBER', 'FLAG_ROW', 'TIME', 'UVW', 'ANTENNA1', 'ANTENNA2')}
                    flags, data, weights, absent, sources, missing_flags, missing_weights = _read_chunk(rows, start, count, ncorr, nchan, datacol)
                    for source in sources:
                        weights_used[source] += 1
                    for c, arr in list(scalars.items()) + [('FLAG', flags), (datacol, data), ('weights', weights), ('absent', absent)]:
                        _hash_array(hash_parts, str(ddid) + ':' + c, arr)
                    for column in ('BITFLAG', 'BITFLAG_ROW'):
                        if column in rows.colnames():
                            values = np.asarray(rows.getcol(column, startrow=start, nrow=count))
                            unique, counts = np.unique(values, return_counts=True)
                            hist = raw_bits.setdefault(column, {})
                            for value, number in zip(unique, counts):
                                hist[str(int(value))] = hist.get(str(int(value)), 0) + int(number)
                            _hash_array(hash_parts, str(ddid) + ':' + column, values)
                    scans = scalars['SCAN_NUMBER']
                    for scan in np.unique(scans):
                        mask = scans == scan
                        stats = eligibility(flags[:, :, mask], scalars['FLAG_ROW'][mask], weights[:, :, mask], data[:, :, mask], selected, absent[mask])
                        stats['absent_samples'][:] = int(np.count_nonzero(absent[mask])) * len(selected)
                        stats['absent_flag_samples'][:] = int(np.count_nonzero(missing_flags[mask])) * len(selected)
                        stats['absent_weight_samples'][:] = int(np.count_nonzero(missing_weights[mask])) * len(selected)
                        stats['nonfinite_weights'] -= stats['absent_weight_samples']
                        scan_stats = by_scan.setdefault(int(scan), {s: empty_stats(len(c)) for s, c in meta['channels'].items()})
                        add_stats(scan_stats[spw], stats)
                        add_stats(aggregate[spw], stats)
            finally:
                rows.close()
    if not by_scan:
        raise ValueError('{}: no cross-correlation rows for field {}'.format(path, fieldname))
    for label, part in sorted(hash_parts.items()):
        hasher.update(label.encode())
        hasher.update(part.digest())
    identity['selected_content_sha256'] = hasher.hexdigest()
    if identity['file_stat_sha256'] != ms_identity(path)['file_stat_sha256']:
        raise RuntimeError('MS changed during read-only survey: ' + path)
    return dict(identity=identity, metadata=meta, aggregate=aggregate, by_scan=by_scan,
                weight_sources=weights_used, raw_flag_bits=raw_bits, correlations=correlations)


def channel_records(survey, scan=None, threshold=0.8):
    stats = survey['aggregate'] if scan is None else survey['by_scan'][scan]
    out = []
    for spw, chans in survey['metadata']['channels'].items():
        for chan in chans:
            entry = dict(chan)
            i = chan['channel']
            entry.update({key: float(stats[spw][key][i]) if key == 'weight_sum' else int(stats[spw][key][i]) for key in STAT_KEYS})
            denominator = entry['denominator']
            entry['usable_fraction'] = entry['joint_usable'] / denominator if denominator else 0.0
            entry['flagged_fraction'] = entry['flagged'] / denominator if denominator else None
            entry['eligible'] = bool(entry['valid_frequency'] and entry['joint_usable'] > 0 and entry['usable_fraction'] >= threshold)
            out.append(entry)
    return sorted(out, key=lambda c: (c['centre_hz'] if c['centre_hz'] is not None else float('inf'), c['spw'], c['channel']))


def quartiles(channels):
    """Half-open slices: floor(N/4), floor(3N/8):floor(5N/8), N-floor(N/4).

    Integer arithmetic avoids banker rounding and float ambiguity. N<4 fails.
    """
    ordered = sorted(channels, key=lambda c: (c['centre_hz'], c['spw'], c['channel']))
    n = len(ordered)
    if n < 4:
        raise ValueError('At least four usable channels are required for quartiles')
    q = n // 4
    return dict(low=ordered[:q], middle=ordered[(3*n)//8:(5*n)//8], high=ordered[n-q:])


def intervals(channels):
    return union([c['interval_hz'] for c in channels if c['interval_hz'] is not None])


def union(ranges):
    result = []
    for lo, hi in sorted(ranges):
        if hi <= lo:
            continue
        if result and lo <= result[-1][1]:
            result[-1][1] = max(hi, result[-1][1])
        else:
            result.append([float(lo), float(hi)])
    return result


def intersection(a, b):
    a, b = union(a), union(b)
    out = []
    i = j = 0
    while i < len(a) and j < len(b):
        lo, hi = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if hi > lo:
            out.append([lo, hi])
        if a[i][1] <= b[j][1]:
            i += 1
        else:
            j += 1
    return union(out)


def bandwidth(ranges):
    return sum(hi-lo for lo, hi in ranges)


def overlap_filter(channels, support, limit):
    """O(N log M) interval coverage; never scan a 32k-channel grid per channel."""
    ranges = union(support)
    starts = [r[0] for r in ranges]
    prefix = [0.0]
    for lo, hi in ranges:
        prefix.append(prefix[-1] + hi-lo)

    def area(value):
        i = bisect.bisect_right(starts, value) - 1
        if i < 0:
            return 0.0
        return prefix[i] + min(value-ranges[i][0], ranges[i][1]-ranges[i][0])

    return [c for c in channels if
            (area(c['interval_hz'][1])-area(c['interval_hz'][0])) / c['width_hz'] + 1e-12 >= limit]


def difference(a, b):
    out = []
    b = union(b)
    j = 0
    for lo, hi in union(a):
        cursor = lo
        while j < len(b) and b[j][1] <= lo:
            j += 1
        k = j
        while k < len(b) and b[k][0] < hi:
            blo, bhi = b[k]
            if blo > cursor:
                out.append([cursor, min(blo, hi)])
            cursor = max(cursor, min(bhi, hi))
            k += 1
        if cursor < hi:
            out.append([cursor, hi])
    return out


def channel_key(c):
    return c['spw'], c['channel']


def common_selection(reference, test, cfg):
    """Reference-defined band; symmetric whole-channel overlap fixed point.

    Remove channels with < min_channel_overlap of their own width supported by
    the other MS. Then require common union bandwidth >= min_common_coverage of
    original reference band AND each retained MS support. Never clip channels.
    """
    if not reference or not test:
        raise ValueError('Inadequate common frequency overlap: empty eligible channel selection')
    frames = {c['frequency_frame'] for c in reference + test}
    if len(frames) != 1:
        raise ValueError('Frequency frames differ; transform to a documented common frame before comparison: ' + str(sorted(frames)))
    initial = intervals(reference)
    ref, other = list(reference), list(test)
    limit = cfg['min_channel_overlap']
    while True:
        a, b = intervals(ref), intervals(other)
        new_ref = overlap_filter(ref, b, limit)
        new_other = overlap_filter(other, a, limit)
        if len(new_ref) == len(ref) and len(new_other) == len(other):
            break
        ref, other = new_ref, new_other
        if not ref or not other:
            raise ValueError('Inadequate common frequency overlap after whole-channel overlap filtering')
    a, b = intervals(ref), intervals(other)
    common = intersection(a, b)
    common_bw = bandwidth(common)
    coverages = dict(requested_reference=common_bw/bandwidth(initial),
                     reference=common_bw/bandwidth(a), test=common_bw/bandwidth(b))
    if min(coverages.values()) + 1e-12 < cfg['min_common_coverage']:
        raise ValueError('Inadequate common frequency coverage: {} (required {})'.format(coverages, cfg['min_common_coverage']))
    detail = dict(rule='whole native channels; symmetric interval-overlap fixed point; no regridding',
        min_channel_overlap=limit, min_common_coverage=cfg['min_common_coverage'],
        frequency_frame=next(iter(frames)), common_intervals_hz=common,
        requested_reference_intervals_hz=initial, coverage_fractions=coverages,
        reference_channel_ids=[list(channel_key(c)) for c in ref],
        test_channel_ids=[list(channel_key(c)) for c in other],
        reference_intervals_hz=a, test_intervals_hz=b,
        unmatched_reference_edges_hz=difference(a, b), unmatched_test_edges_hz=difference(b, a),
        excluded_reference_intervals_hz=difference(initial, a),
        excluded_test_intervals_hz=difference(intervals(test), b),
        exact_frequency_support=(a == b), common_bandwidth_hz=common_bw)
    return ref, other, detail


def spw_selector(channels):
    groups = {}
    for c in channels:
        groups.setdefault(c['spw'], set()).add(c['channel'])
    if not groups:
        raise ValueError('Empty SPW channel selection must never mean all channels')
    chunks = []
    for spw, ids in sorted(groups.items()):
        runs = []
        start = previous = None
        for chan in sorted(ids):
            if previous is not None and chan != previous + 1:
                runs.append(str(start) if start == previous else '{}~{}'.format(start, previous))
                start = chan
            elif start is None:
                start = chan
            previous = chan
        runs.append(str(start) if start == previous else '{}~{}'.format(start, previous))
        chunks.append('{}:{}'.format(spw, ';'.join(runs)))
    return ','.join(chunks)


def selection_summary(channels):
    supported = [c for c in channels if c['joint_usable'] > 0 and c['valid_frequency']]
    ranges = intervals(supported)
    sums = {key: sum(c[key] for c in channels) for key in STAT_KEYS}
    denom = sums['denominator']
    weight = sums['weight_sum']
    return dict(channel_count=len(channels), contributing_channel_count=len(supported),
        contributed_channel_ids=[[c['spw'], c['channel']] for c in supported],
        actual_frequency_intervals_hz=ranges, actual_frequency_span_hz=[ranges[0][0], ranges[-1][1]] if ranges else None,
        actual_bandwidth_hz=bandwidth(ranges),
        weighted_effective_frequency_hz=sum(c['centre_hz'] * c['weight_sum'] for c in supported)/weight if weight > 0 else None,
        usable_fraction=sums['joint_usable']/denom if denom else 0.0,
        flagged_fraction=sums['flagged']/denom if denom else None,
        flag_and_weight_summary=sums)


def task_defaults(task):
    defaults = {}
    for name, param in inspect.signature(task).parameters.items():
        if param.default is not inspect.Parameter.empty:
            defaults[name] = copy.deepcopy(param.default)
    if not {'vis', 'imagename', 'stokes', 'specmode'}.issubset(defaults):
        raise RuntimeError('Cannot introspect complete CASA tclean defaults; use a CASA 6 task with its public signature')
    return defaults


def slug(name):
    return re.sub('[^A-Za-z0-9._-]', '_', name) + '_' + hashlib.sha256(name.encode()).hexdigest()[:8]


def imaging_dependencies(params):
    """Record mask/voltage-pattern inputs, including contents at fixed paths."""
    result = {}
    for key in ('mask', 'vptable'):
        values = params.get(key, '')
        values = values if isinstance(values, (list, tuple)) else [values]
        records = []
        for value in values:
            if not isinstance(value, str) or not value or not os.path.exists(os.path.expanduser(value)):
                continue  # inline region expression or task-resolved value
            path = os.path.realpath(os.path.expanduser(value))
            if os.path.isdir(path):
                records.append(ms_identity(path))
            else:
                hasher = hashlib.sha256()
                with open(path, 'rb') as handle:
                    for block in iter(lambda: handle.read(1024*1024), b''):
                        hasher.update(block)
                records.append(dict(path=path, content_sha256=hasher.hexdigest()))
        if records:
            result[key] = records
    return result


def build_plan(cfg, factory, task, version):
    cfg = validate_config(cfg)
    if not cfg['enabled']:
        raise ValueError('paired_astrometry.enabled must be true')
    if not str(version).startswith('6.'):
        raise RuntimeError('Paired astrometry requires CASA 6; found ' + str(version))
    defaults = task_defaults(task)
    unknown = set(cfg['tclean']) - set(defaults)
    if unknown:
        raise ValueError('Unsupported tclean parameters in this CASA version: ' + ', '.join(sorted(unknown)))
    output = os.path.realpath(os.path.expanduser(cfg['output_dir']))
    surveys = {}
    diagnostics = []
    products = []
    field_plans = []
    dependency_cache = None
    # Survey everything before imaging: all fields/MSs must pass the overlap contract.
    for field in cfg['fields']:
        pair = {}
        for role in ('reference', 'test'):
            path = os.path.realpath(os.path.expanduser(field[role + '_ms']))
            if os.path.commonpath([output, path]) == path:
                raise ValueError('output_dir must not be inside an input MS')
            cache_key = (path, field['name'])
            if cache_key not in surveys:
                print('Surveying {} field={} (read-only)'.format(path, field['name']), flush=True)
                surveys[cache_key] = survey_ms(path, field['name'], cfg, factory)
            pair[role] = surveys[cache_key]
        native = {role: channel_records(s, threshold=cfg['occupancy_threshold']) for role, s in pair.items()}
        usable = {role: [c for c in channels if c['eligible']] for role, channels in native.items()}
        native_bands = {role: quartiles(channels) if len(channels) >= 4 else {} for role, channels in usable.items()}
        for role, s in pair.items():
            diagnostics.append(dict(ms_identity=s['identity'], role=role, field_name=field['name'], field_kind=field['kind'],
                field_id=s['metadata']['field_id'], phase_centre=s['metadata']['phase_centre'], scan_ids=sorted(s['by_scan']),
                occupancy_threshold=cfg['occupancy_threshold'], occupancy_policy=POLICY,
                selected_correlations=s['correlations'], weight_sources=s['weight_sources'], raw_flag_bits=s['raw_flag_bits'],
                channels=native[role], native_usable_summary=selection_summary(usable[role]),
                native_quartile_summaries={band: selection_summary(chans) for band, chans in native_bands[role].items()},
                native_quartile_channel_ids={band: [list(channel_key(c)) for c in chans] for band, chans in native_bands[role].items()}))
        field_plans.append((field, pair, usable, native_bands))
    for field, pair, usable, native_bands in field_plans:
        if any(not bands for bands in native_bands.values()):
            raise SelectionPlanningError('Field {}: fewer than four native usable channels; inspect native diagnostics'.format(field['name']), diagnostics)
        bands = dict(native_bands['reference'])
        bands['fullband'] = usable['reference']
        phase = pair['reference']['metadata']['phase_centre']
        phase_string = '{} {:.17g}rad {:.17g}rad'.format(phase['frame'], phase['ra_rad'], phase['dec_rad'])
        for band, ref_channels in bands.items():
            try:
                ref, other, common = common_selection(ref_channels, usable['test'], cfg)
            except ValueError as error:
                raise SelectionPlanningError('Field {} product {}: {}'.format(field['name'], band, error), diagnostics) from error
            for role, selected in (('reference', ref), ('test', other)):
                s = pair[role]
                scans = sorted(s['by_scan'])
                contexts = [[scan] for scan in scans] if band == 'fullband' else [scans]
                for scan_ids in contexts:
                    actual = channel_records(s, scan=scan_ids[0] if band == 'fullband' else None, threshold=cfg['occupancy_threshold'])
                    lookup = {channel_key(c): c for c in actual}
                    requested = [lookup[channel_key(c)] for c in selected]
                    summary = selection_summary(requested)
                    actual_common = intersection(summary['actual_frequency_intervals_hz'], common['common_intervals_hz'])
                    coverage = bandwidth(actual_common) / common['common_bandwidth_hz']
                    basename = '{}_{}_{}{}'.format(role, slug(field['name']), band,
                                  '_scan{}'.format(scan_ids[0]) if band == 'fullband' else '')
                    name = os.path.join(output, 'images', basename)
                    params = copy.deepcopy(defaults)
                    params.update(IMAGING_DEFAULTS)
                    params.update(cfg['tclean'])
                    params.update(vis=s['identity']['path'], imagename=name, field=str(s['metadata']['field_id']),
                        scan=','.join(map(str, scan_ids)), spw=spw_selector(selected),
                        datacolumn=cfg['datacolumn'], stokes='I', specmode='mfs', antenna='*&*',
                        phasecenter=phase_string, selectdata=True, restart=False, savemodel='none',
                        parallel=False, interactive=False, restoration=True, calcpsf=True, calcres=True)
                    if dependency_cache is None:
                        dependency_cache = imaging_dependencies(params)
                    contract = dict(schema_version=SCHEMA, ms_identity=s['identity'], ms_role=role,
                        paired_inputs={r: dict(ms_identity=sv['identity'],
                            field_id=sv['metadata']['field_id'], field_name=field['name']) for r, sv in pair.items()},
                        field_id=s['metadata']['field_id'], field_name=field['name'], field_kind=field['kind'], scan_ids=scan_ids,
                        product=band, requested_channels=selected,
                        selected_channels_with_actual_data=requested, selection=summary,
                        occupancy_threshold=cfg['occupancy_threshold'], occupancy_policy=POLICY,
                        selected_correlations=s['correlations'], weight_sources=s['weight_sources'], raw_flag_bits=s['raw_flag_bits'],
                        input_phase_centre=s['metadata']['phase_centre'], imaging_phase_centre=phase,
                        common_frequency_selection=common, scan_actual_common_coverage_fraction=coverage,
                        frequency_comparability=('aggregate interval-overlap contract' if band != 'fullband' else
                            'common requested fullband; actual scan support is independent; scan IDs do not establish temporal pairing'),
                        per_scan_support_incomplete=(band == 'fullband' and coverage < 1.0 - 1e-12),
                        actual_support_comparability=('not verified across independent scans' if band == 'fullband' else
                            ('identical aggregate native interval support' if common['exact_frequency_support'] else 'overlap contract only; native support differs')),
                        imaging_input_dependencies=dependency_cache,
                        casa_version=version, tclean_parameters=params, data_column=('DATA' if cfg['datacolumn'] == 'data' else 'CORRECTED_DATA'),
                        output_paths=dict(casa_image=name+'.image', fits=name+'.fits', qa=name+'_qa.txt', manifest=name+'.manifest.json'))
                    products.append(dict(contract=contract, fingerprint=digest(contract), processing_status='planned',
                        fits_wcs_centre=None, restoring_beam=None, contribution_status='predicted_from_MS_flags_weights_finite_DATA; gridding contribution unmeasured'))
    return dict(schema_version=SCHEMA, casa_version=version, config=cfg, native_diagnostics=diagnostics, images=products)


def output_state(product):
    """Only a completed manifest with identical contract permits reuse."""
    contract = product['contract']
    paths = contract['output_paths']
    prefix = contract['tclean_parameters']['imagename']
    existing = [p for p in glob.glob(prefix + '.*') if not p.endswith('.manifest.json') and not p.endswith('.tmp')]
    if not existing and not os.path.exists(paths['manifest']):
        return 'new'
    try:
        with open(paths['manifest']) as handle:
            previous = json.load(handle)
        if (previous['fingerprint'] == product['fingerprint'] and
                digest(previous['contract']) == product['fingerprint'] and
                previous['processing_status'] == 'complete' and
                os.path.exists(paths['casa_image']) and os.path.isfile(paths['fits'])):
            return 'reuse'
    except (OSError, ValueError, KeyError, TypeError):
        pass
    raise RuntimeError('Stale or incomplete output at {}. Inputs/parameters/provenance differ; choose a new output_dir or explicitly remove the old product family.'.format(prefix))


def fits_metadata(path):
    from astropy.io import fits
    from astropy.wcs import WCS
    with fits.open(path, memmap=True) as hdus:
        h = hdus[0].header
        wcs = WCS(h).celestial
        x, y = (h['NAXIS1']-1)/2.0, (h['NAXIS2']-1)/2.0
        ra, dec = wcs.all_pix2world([[x, y]], 0)[0]
        if not np.isfinite([ra, dec]).all():
            raise RuntimeError('Nonfinite FITS WCS centre')
        centre = dict(ra_deg=float(ra), dec_deg=float(dec),
            frame=h.get('RADESYS', h.get('RADECSYS', 'unspecified')), equinox=h.get('EQUINOX'),
            definition='geometric image centre, zero-based pixel ((NAXIS1-1)/2,(NAXIS2-1)/2)',
            reference_pixel=[float(h['CRPIX1']), float(h['CRPIX2'])],
            reference_world_deg=[float(h['CRVAL1']), float(h['CRVAL2'])],
            ctype=[h['CTYPE1'], h['CTYPE2']])
        beam = dict(major_arcsec=float(h['BMAJ'])*3600, minor_arcsec=float(h['BMIN'])*3600, pa_deg=float(h['BPA']))
        if not np.isfinite(list(beam.values())).all() or beam['major_arcsec'] <= 0 or beam['minor_arcsec'] <= 0:
            raise RuntimeError('Missing/invalid restoring beam in FITS')
    return centre, beam


def execute_plan(plan, task, export, qa=None):
    """Preflight every output before any expensive tclean; never auto-delete."""
    cfg = plan['config']
    output = os.path.realpath(os.path.expanduser(cfg['output_dir']))
    # A dry run uses a separate plan file and never overwrites completed image manifests.
    write_json(os.path.join(output, 'selection_plan.json'), plan)
    write_json(os.path.join(output, 'native_channel_diagnostics.json'), plan['native_diagnostics'])
    for p in plan['images']:
        c = p['contract']
        s = c['selection']
        print('{} {} {} scans={} channels={} effective={:.6f} MHz coverage={}'.format(
            c['ms_role'], c['field_name'], c['product'], c['scan_ids'], s['channel_count'],
            (s['weighted_effective_frequency_hz'] or 0)/1e6, c['common_frequency_selection']['coverage_fractions']), flush=True)
    if cfg['dry_run']:
        print('Selection-only run: no imaging or FITS export performed.', flush=True)
        return
    states = [output_state(p) for p in plan['images']]
    failures = []
    for p, state in zip(plan['images'], states):
        c = p['contract']
        paths = c['output_paths']
        if state == 'reuse':
            print('Reusing verified product: ' + paths['fits'], flush=True)
            with open(paths['manifest']) as handle:
                p.update(json.load(handle))
            p['reused'] = True
            continue
        p['processing_status'] = 'running'
        p['started_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
        write_json(paths['manifest'], p)
        try:
            if c['selection']['contributing_channel_count'] == 0:
                raise RuntimeError('No unflagged positive-weight data in requested scan selection')
            task(**c['tclean_parameters'])
            if not os.path.isdir(paths['casa_image']):
                raise RuntimeError('tclean did not produce the expected .image')
            export(imagename=paths['casa_image'], fitsimage=paths['fits'], overwrite=False, dropdeg=False, stokeslast=True)
            p['fits_wcs_centre'], p['restoring_beam'] = fits_metadata(paths['fits'])
            if qa:
                qa(paths['fits'])
            p['processing_status'] = 'complete'
        except Exception as error:
            p['processing_status'] = 'failed'
            p['error'] = '{}: {}'.format(type(error).__name__, error)
            failures.append(p['error'])
        finally:
            p['finished_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
            write_json(paths['manifest'], p)
    write_json(os.path.join(output, 'run_manifest.json'), plan)
    if failures:
        raise RuntimeError('Paired imaging failed: ' + '; '.join(failures))


def run(cfg, factory=None, task=None, export=None, qa=None):
    import casatasks
    from casatools import table
    version = casatasks.version_string()
    try:
        plan = build_plan(cfg, factory or table, task or casatasks.tclean, version)
    except SelectionPlanningError as error:
        output = os.path.realpath(os.path.expanduser(cfg['output_dir']))
        write_json(os.path.join(output, 'native_channel_diagnostics.json'), error.diagnostics)
        write_json(os.path.join(output, 'planning_failure.json'), dict(processing_status='failed',
            phase='selection', error=str(error), casa_version=version, config=cfg,
            finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
        raise
    execute_plan(plan, task or casatasks.tclean, export or casatasks.exportfits, qa)
    return plan
