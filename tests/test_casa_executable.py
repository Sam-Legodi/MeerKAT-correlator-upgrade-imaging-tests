"""Configured CASA selection applies to calibration and both imaging modes."""
import pytest

from meerkat_corr_imaging.config import _dict_to_dataclass
from meerkat_corr_imaging.steps import step3_calibrate_image as step


@pytest.mark.parametrize('value', [None, '', '  ', False, 6])
def test_invalid_casa_executable(value):
    with pytest.raises(ValueError, match='casa.executable'):
        _dict_to_dataclass(dict(project_name='test', casa=dict(executable=value)))


@pytest.mark.parametrize('configured,override,expected', [
    (None, None, 'casa'),
    ('/opt/casa6/bin/casa', None, '/opt/casa6/bin/casa'),
    ('/opt/casa6/bin/casa', '/other/casa', '/other/casa'),
])
def test_executable_used_for_legacy_calibration_and_imaging(tmp_path, monkeypatch, configured, override, expected):
    ms = tmp_path/'input.ms'; ms.mkdir()
    casa = {} if configured is None else dict(executable=configured)
    cfg = _dict_to_dataclass(dict(project_name='test', casa=casa,
        reference=dict(ms_paths=[str(ms)]), extra=dict(force_calibrate=True)))
    if override is None:
        monkeypatch.delenv('CASA', raising=False)
    else:
        monkeypatch.setenv('CASA', override)
    calls = []
    monkeypatch.setattr(step, 'run_calibration', lambda cmd, *a, **kw: calls.append(cmd))
    monkeypatch.setattr(step, '_run_cmd', lambda cmd, **kw: calls.append(cmd))
    step.run(cfg)
    assert len(calls) == 2
    assert all(c[0] == expected for c in calls)


def test_executable_used_for_target_only_pair(tmp_path, monkeypatch):
    paths = [tmp_path/'reference.ms', tmp_path/'test.ms']
    for p in paths:
        p.mkdir()
    cfg = _dict_to_dataclass(dict(project_name='test', casa=dict(
        executable='/opt/casa6/bin/casa', paired_astrometry=dict(enabled=True,
            experiment_mode='target_only', output_dir=str(tmp_path/'output'),
            fields=[dict(name='J2147-8132', kind='target',
                reference_ms=str(paths[0]), test_ms=str(paths[1]))]))))
    monkeypatch.delenv('CASA', raising=False)
    calls = []
    monkeypatch.setattr(step, '_run_cmd', lambda cmd, **kw: calls.append(cmd))
    step.run(cfg)
    assert len(calls) == 1 and calls[0][0] == '/opt/casa6/bin/casa'
