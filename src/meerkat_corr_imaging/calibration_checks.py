"""Task and output checks compatible with CASA 5's Python 2.7."""
from __future__ import print_function
import os
import time
import numpy as np


def table_tool():
    try:
        from casatools import table
    except ImportError:
        from taskinit import tbtool as table
    return table()


def inspect_solutions(path):
    tool = table_tool()
    try:
        tool.open(path)
        if not tool.nrows():
            raise RuntimeError('Empty calibration table: ' + path)
        column = 'CPARAM' if 'CPARAM' in tool.colnames() else 'FPARAM'
        total = usable = 0
        coverage = {}
        for row in range(tool.nrows()):
            values = np.asarray(tool.getcell(column, row))
            flags = np.asarray(tool.getcell('FLAG', row), dtype=bool)
            good = ~flags & np.isfinite(values)
            key = '%s/%s/%s' % tuple(tool.getcell(c, row) for c in
                                    ('ANTENNA1', 'SPECTRAL_WINDOW_ID', 'FIELD_ID'))
            coverage[key] = coverage.get(key, 0) + int(good.sum())
            total += values.size
            usable += int(good.sum())
        if not usable:
            raise RuntimeError('No finite unflagged solutions: ' + path)
        return {'path': path, 'usable_fraction': float(usable)/total, 'coverage': coverage}
    finally:
        tool.done()


def inspect_corrected(vis, fields):
    tool = table_tool()
    try:
        tool.open(vis)
        if 'CORRECTED_DATA' not in tool.colnames():
            raise RuntimeError('CORRECTED_DATA missing')
        counts = {}
        # Bounded sampling per field, rather than loading an entire MS cube.
        for field in fields.split(','):
            # CASA 5 cannot reliably query an unnamed temporary query table.
            # Apply both predicates directly to the open MeasurementSet.
            rows = tool.query('FIELD_ID==%s AND ANTENNA1!=ANTENNA2' % field)
            try:
                good_count = 0
                for row in range(0, rows.nrows(), max(1, rows.nrows() // 256)):
                    values = np.asarray(rows.getcell('CORRECTED_DATA', row))
                    flags = np.asarray(rows.getcell('FLAG', row), dtype=bool)
                    if rows.getcell('FLAG_ROW', row):
                        continue
                    if np.any(~flags & ~np.isfinite(values)):
                        raise RuntimeError('Nonfinite unflagged corrected data in field ' + field)
                    good_count += int((~flags).sum())
                if not good_count:
                    raise RuntimeError('No usable corrected samples in field ' + field)
                counts[field] = good_count
            finally:
                rows.close()
        return counts
    finally:
        tool.done()


def quality_residual(vis, flux_field, applied_fields):
    names = table_tool()
    try:
        names.open(os.path.join(vis, 'FIELD'))
        field_names = list(names.getcol('NAME'))
    finally:
        names.done()
    ids = []
    for token in str(flux_field).split(','):
        ids.extend([int(token)] if token.isdigit() else
                   [i for i, name in enumerate(field_names) if name == token])
    if not ids:
        raise RuntimeError('Quality check flux field not found')
    if not set(ids).issubset(set(map(int, applied_fields.split(',')))):
        raise RuntimeError('Quality flux field excluded from calibration application')
    tool = table_tool()
    rows = None
    try:
        tool.open(vis)
        if 'MODEL_DATA' not in tool.colnames():
            raise RuntimeError('Quality check requires MODEL_DATA')
        rows = tool.query('FIELD_ID IN [%s] AND ANTENNA1!=ANTENNA2' % ','.join(map(str, ids)))
        residuals = []
        for row in range(0, rows.nrows(), max(1, rows.nrows() // 256)):
            if rows.getcell('FLAG_ROW', row):
                continue
            data = np.asarray(rows.getcell('CORRECTED_DATA', row))
            model = np.asarray(rows.getcell('MODEL_DATA', row))
            flags = np.asarray(rows.getcell('FLAG', row), dtype=bool)
            good = ~flags & np.isfinite(data) & np.isfinite(model) & (np.abs(model) > 0)
            if good.any():
                residuals.extend((np.abs(data[good]-model[good])/np.abs(model[good])).tolist())
        if not residuals:
            raise RuntimeError('No usable modeled calibrator samples for quality check')
        return float(np.median(residuals))
    finally:
        if rows is not None:
            rows.close()
        tool.done()


class CheckedTask(object):
    """Keep CASA 5 CLI metadata (defaults, parameters, etc.) accessible."""
    def __init__(self, original, checked):
        self._original = original
        self._checked = checked

    def __call__(self, **kwargs):
        return self._checked(**kwargs)

    def __getattr__(self, name):
        return getattr(self._original, name)


def install_checks(namespace, result):
    quality = os.environ.get('MCI_QUALITY_ENABLED') == '1'
    threshold = float(os.environ.get('MCI_QUALITY_MIN_SOLUTION_FRACTION', '0.95'))
    result['quality'] = 'pending' if quality else 'not_run'
    result['solution_tables'] = []
    def wrap(name, task):
        def checked(**kwargs):
            print('[CAL TASK] Starting ' + name)
            import sys
            sys.stdout.flush()
            logger = namespace.get('casalog')
            if logger is None:
                try:
                    from taskinit import casalog as logger
                except ImportError:
                    from casatasks import casalog as logger
            logfile = logger.logfile() if logger is not None else None
            offset = os.path.getsize(logfile) if logfile and os.path.exists(logfile) else 0
            started = time.time()
            if name == 'applycal':
                # Flag samples without applicable solutions so unchanged DATA
                # cannot masquerade as successfully corrected unflagged samples.
                kwargs['applymode'] = 'calonly' # 'calflagstrict'
            if name == 'setjy' and quality:
                kwargs['usescratch'] = True
            value = task(**kwargs)
            if value is False:
                raise RuntimeError(name + ' returned False')
            if logfile:
                with open(logfile) as handle:
                    handle.seek(offset)
                    lines = handle.readlines()
                fatal = [line for line in lines if ('SEVERE' in line or 'ERROR' in line)
                         and not any(text in line for text in ('MeasTable::dUTC', 'Leap second', 'TAI_UTC'))]
                if fatal:
                    raise RuntimeError(name + ' logged errors: ' + ''.join(fatal[-5:]))
            output = kwargs.get('fluxtable') or kwargs.get('caltable')
            if output:
                result['solution_tables'].append(inspect_solutions(output))
            if name == 'applycal':
                result['corrected_samples_by_field'] = inspect_corrected(
                    kwargs['vis'], kwargs['field'])
                result['application_statistics'] = [line.strip() for line in lines
                    if any(term in line.lower() for term in ('calibration apply', 'flagged', 'jones'))]
                result['calibration_applied'] = True
                print('Calibration applied successfully (task/output checks; sampled data validation).')
                sys.stdout.flush()
                if quality:
                    fractions = [entry['usable_fraction'] for entry in result['solution_tables']]
                    result['quality'] = 'failed'
                    residual = quality_residual(kwargs['vis'], os.environ['MCI_CAL_FLUX_FIELD'], kwargs['field'])
                    limit = float(os.environ.get('MCI_QUALITY_MAX_RESIDUAL', '0.1'))
                    passed = bool(fractions) and min(fractions) >= threshold and residual <= limit
                    result['quality_metrics'] = {'minimum_solution_fraction': min(fractions) if fractions else 0,
                                                 'required_fraction': threshold, 'median_fractional_residual': residual, 'maximum_residual': limit}
                    result['quality'] = 'passed' if passed else 'failed'
                    reasons = []
                    if not fractions:
                        reasons.append('no solution tables available')
                    for entry in result['solution_tables']:
                        if entry['usable_fraction'] < threshold:
                            reasons.append('%s usable solution fraction %.6g < %.6g' %
                                           (entry.get('path', 'unnamed table'), entry['usable_fraction'], threshold))
                    if not np.isfinite(residual) or residual > limit:
                        reasons.append('median fractional complex residual %.6g exceeds limit %.6g or is nonfinite' % (residual, limit))
                    result['quality_failure_reasons'] = reasons
                    print('[CAL QUALITY] minimum solution fraction=%s (required >= %s); median residual=%s (required <= %s)' %
                          (min(fractions) if fractions else 0, threshold, residual, limit))
                    sys.stdout.flush()
                    if not passed:
                        raise RuntimeError('Calibration quality failed: ' + '; '.join(reasons))
                    print('Calibration quality passed (configured solution coverage and sampled calibrator residual checks).')
            result['tasks'].append({'name': name, 'seconds': time.time()-started, 'status': 'passed'})
            print('[CAL TASK] Finished ' + name)
            sys.stdout.flush()
            return value
        return checked
    for name in ('delmod', 'setjy', 'gaincal', 'bandpass', 'fluxscale', 'polcal', 'applycal', 'split'):
        task = namespace.get(name)
        if task is None:
            try:
                import casatasks
                task = getattr(casatasks, name)
            except ImportError:
                raise RuntimeError('Missing CASA task: ' + name)
        namespace[name] = CheckedTask(task, wrap(name, task))
