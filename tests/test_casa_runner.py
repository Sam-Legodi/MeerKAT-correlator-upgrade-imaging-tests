import json
import os
import sys
import pytest
from meerkat_corr_imaging.casa_runner import supervise, run_calibration


def test_no_newline_and_clean_exit(capsys):
    supervise([sys.executable, '-c', "print('prompt', end='')"], dict(os.environ), timeout=3)
    assert 'prompt' in capsys.readouterr().out


def test_timeout():
    with pytest.raises(TimeoutError):
        supervise([sys.executable, '-c', 'import time; time.sleep(30)'], dict(os.environ), timeout=.1)


def test_missing_result(tmp_path):
    with pytest.raises(RuntimeError, match='without a calibration result'):
        run_calibration([sys.executable, '-c', 'pass'], dict(os.environ, MCI_CAL_MSFILE='x'), tmp_path, 3, 1)


def test_matching_result(tmp_path):
    script = "import os,json; json.dump(dict(run_id=os.environ['MCI_RUN_ID'],ms='x',calibration_applied=True,quality='not_run'),open(os.environ['MCI_RESULT_PATH'],'w'))"
    result = run_calibration([sys.executable, '-c', script], dict(os.environ, MCI_CAL_MSFILE='x'), tmp_path, 3, 1)
    assert result['calibration_applied']


def test_completion_hang(tmp_path):
    path = tmp_path/'result.json'
    script = "import time; open(%r,'w').write('{}'); time.sleep(30)" % str(path)
    with pytest.raises(TimeoutError):
        supervise([sys.executable, '-c', script], dict(os.environ, MCI_RESULT_PATH=str(path)), 5, .1)


def test_failed_task_stops(monkeypatch):
    from meerkat_corr_imaging.calibration_checks import install_checks
    names = ('delmod','setjy','gaincal','bandpass','fluxscale','polcal','applycal','split')
    ns = {name: lambda **kw: None for name in names}
    from types import SimpleNamespace
    ns['casalog'] = SimpleNamespace(logfile=lambda: None)
    ns['gaincal'] = lambda **kw: False
    result = {'tasks': []}
    install_checks(ns, result)
    with pytest.raises(RuntimeError, match='returned False'):
        ns['gaincal']()
    assert not result['tasks']


def test_log_errors_and_warning_exception(tmp_path):
    from types import SimpleNamespace
    from meerkat_corr_imaging.calibration_checks import install_checks
    log = tmp_path/'casa.log'
    log.write_text('')
    names = ('delmod','setjy','gaincal','bandpass','fluxscale','polcal','applycal','split')
    ns = {name: lambda **kw: None for name in names}
    ns['casalog'] = SimpleNamespace(logfile=lambda: str(log))
    def warning(**kw):
        with log.open('a') as f:
            f.write('SEVERE MeasTable::dUTC Leap second out of date\n')
    ns['setjy'] = warning
    def failure(**kw):
        with log.open('a') as f:
            f.write('SEVERE gaincal solve failed\n')
    ns['gaincal'] = failure
    result = {'tasks': []}
    install_checks(ns, result)
    ns['setjy']()
    with pytest.raises(RuntimeError, match='logged errors'):
        ns['gaincal']()
    assert [t['name'] for t in result['tasks']] == ['setjy']


def test_application_and_quality_are_separate(monkeypatch, tmp_path):
    from types import SimpleNamespace
    import meerkat_corr_imaging.calibration_checks as checks
    log = tmp_path/'casa.log'
    log.write_text('')
    monkeypatch.setenv('MCI_QUALITY_ENABLED', '1')
    monkeypatch.setenv('MCI_CAL_FLUX_FIELD', 'flux')
    monkeypatch.setattr(checks, 'inspect_corrected', lambda *a: {'0': 100})
    monkeypatch.setattr(checks, 'quality_residual', lambda *a: 0.5)
    names = ('delmod','setjy','gaincal','bandpass','fluxscale','polcal','applycal','split')
    ns = {name: lambda **kw: None for name in names}
    ns['casalog'] = SimpleNamespace(logfile=lambda: str(log))
    result = {'tasks': []}
    checks.install_checks(ns, result)
    result['solution_tables'] = [{'usable_fraction': 1.0}]
    with pytest.raises(RuntimeError, match='quality failed'):
        ns['applycal'](vis='x', field='0')
    assert result['calibration_applied'] is True
    assert result['quality'] == 'failed'


@pytest.mark.parametrize('flags,empty', [([True, True], False), ([], True)])
def test_invalid_solution_tables(monkeypatch, flags, empty):
    import numpy as np
    import meerkat_corr_imaging.calibration_checks as checks
    class Table:
        closed = False
        def open(self, path): pass
        def nrows(self): return 0 if empty else 1
        def colnames(self): return ['CPARAM']
        def getcell(self, name, row):
            return {'CPARAM': np.array([1., 1.]), 'FLAG': np.array(flags),
                    'ANTENNA1': 0, 'SPECTRAL_WINDOW_ID': 0, 'FIELD_ID': 0}[name]
        def done(self): self.closed = True
    table = Table()
    monkeypatch.setattr(checks, 'table_tool', lambda: table)
    with pytest.raises(RuntimeError):
        checks.inspect_solutions('x')
    assert table.closed


def test_stage_timeout():
    with pytest.raises(TimeoutError, match='stage'):
        supervise([sys.executable, '-u', '-c', "import time; print('[CAL TASK] Starting gaincal'); time.sleep(30)"],
                  dict(os.environ, MCI_STAGE_TIMEOUT_SECONDS='0.1'), timeout=5)
