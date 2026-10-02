"""Selection/provenance tests use fixed-shape CASA table/task mocks, no CASA."""
import copy
import json
from pathlib import Path

import numpy as np
import pytest

from meerkat_corr_imaging import paired_astrometry as pa
from meerkat_corr_imaging.config import _dict_to_dataclass
from meerkat_corr_imaging.steps import step3_calibrate_image


def channel(chan, freq=None, width=1.0, spw=0, **updates):
    freq = float(chan) if freq is None else freq
    result = dict(spw=spw, channel=chan, centre_hz=freq, width_hz=width,
                  interval_hz=[freq-width/2, freq+width/2], frequency_frame='TOPO', valid_frequency=True,
                  eligible=True, usable_fraction=1.0, flagged_fraction=0.0)
    result.update({k: 2.0 if k == 'weight_sum' else (2 if k in ('denominator', 'joint_usable', 'positive_finite_weights') else 0) for k in pa.STAT_KEYS})
    result.update(updates)
    return result


@pytest.mark.parametrize('n,low,middle,high', [
    (4, [0], [1], [3]), (5, [0], [1, 2], [4]),
    (7, [0], [2, 3], [6]), (8, [0, 1], [3, 4], [6, 7]),
    (9, [0, 1], [3, 4], [7, 8]),
    (12, [0, 1, 2], [4, 5, 6], [9, 10, 11]),
    (16, list(range(4)), list(range(6, 10)), list(range(12, 16))),
])
def test_deterministic_quarter_boundaries(n, low, middle, high):
    bands = pa.quartiles([channel(i) for i in reversed(range(n))])
    assert [[c['channel'] for c in bands[b]] for b in ('low', 'middle', 'high')] == [low, middle, high]


def test_quartiles_count_channels_and_preserve_irregular_gaps():
    chans = [channel(i, f) for i, f in enumerate([20, 21, 24, 28, 80, 81, 90, 100])]
    assert [c['centre_hz'] for c in pa.quartiles(chans)['middle']] == [28, 80]
    assert pa.intervals(pa.quartiles(chans)['middle']) == [[27.5, 28.5], [79.5, 80.5]]
    assert pa.spw_selector([channel(i) for i in (0, 1, 4, 8, 9)] + [channel(2, spw=1)]) == '0:0~1;4;8~9,1:2'
    with pytest.raises(ValueError, match='four'):
        pa.quartiles([channel(0)])
    with pytest.raises(ValueError, match='Empty'):
        pa.spw_selector([])


def test_eligibility_flags_weights_and_strict_I():
    flags = np.zeros((4, 7, 5), bool)
    weights = np.ones(flags.shape)
    data = np.ones(flags.shape, complex)
    flags[:, 0, :] = True # wholly flagged
    weights[:, 1, :] = 0 # wholly zero weight
    weights[:, 2, :] = np.nan
    weights[:, 3, :] = -1
    flags[0, 4, 0] = True # one parallel hand: both opportunities excluded
    flags[1, 5, 0] = True # cross hand: strict I conservative flag policy
    data[3, 6, 0] = np.nan
    stats = pa.eligibility(flags, [False]*5, weights, data, [0, 3])
    assert stats['denominator'].tolist() == [10]*7
    assert stats['joint_usable'].tolist() == [0, 0, 0, 0, 8, 8, 8]
    assert stats['weight_sum'].tolist() == [0, 0, 0, 0, 8, 8, 8]
    assert stats['flagged'][4] == 1
    assert stats['flagged'][5] == 0 # separate crosshand restriction
    flagged_rows = pa.eligibility(flags, [True]*5, weights, data, [0, 3])
    assert not flagged_rows['joint_usable'].any()
    with pytest.raises(ValueError, match='shapes'):
        pa.eligibility(flags, [False]*5, weights[:, :2], data, [0, 3])


def config(tmp_path):
    fields = []
    for name, kind in [('gain', 'gain_calibrator'), ('target', 'target')]:
        paths = {}
        for role in ('reference', 'test'):
            ms = tmp_path / (role + '_' + name + '.ms')
            ms.mkdir()
            (ms / 'table.dat').write_text('metadata')
            paths[role + '_ms'] = str(ms)
        fields.append(dict(name=name, kind=kind, **paths))
    return pa.validate_config(dict(enabled=True, output_dir=str(tmp_path / 'output'), fields=fields))


def test_paired_exact_and_different_widths():
    cfg = pa.validate_config({})
    ref = [channel(i, i+0.5) for i in range(8)]
    other = [channel(i, 2*i+1, width=2, spw=2) for i in range(4)]
    a, b, detail = pa.common_selection(ref, other, cfg)
    assert len(a) == 8 and len(b) == 4
    assert detail['exact_frequency_support']
    assert detail['coverage_fractions'] == dict(requested_reference=1, reference=1, test=1)
    # Different centres leave measurable unmatched native edges.
    shifted = [channel(i, i+0.55, spw=3) for i in range(8)]
    a, b, detail = pa.common_selection(ref, shifted, cfg)
    assert not detail['exact_frequency_support']
    assert detail['unmatched_reference_edges_hz'] == [[0, pytest.approx(.05)]]
    assert detail['unmatched_test_edges_hz'] == [[8, pytest.approx(8.05)]]
    assert detail['coverage_fractions']['reference'] == pytest.approx(7.95/8)


def test_common_intersection_excludes_gaps_and_inadequate_overlap():
    cfg = pa.validate_config({})
    ref = [channel(i) for i in range(16)]
    test = [channel(i, spw=1) for i in range(16) if i != 4]
    a, b, detail = pa.common_selection(ref, test, cfg)
    assert 4 not in [c['channel'] for c in a]
    assert detail['coverage_fractions']['requested_reference'] == 15/16
    assert [3.5, 4.5] not in detail['common_intervals_hz']
    with pytest.raises(ValueError, match='coverage'):
        pa.common_selection(ref[2:6], test[:4], cfg)
    with pytest.raises(ValueError, match='overlap'):
        pa.common_selection(ref, [channel(0, 200)], cfg)
    with pytest.raises(ValueError, match='frames'):
        pa.common_selection(ref, [channel(0, frequency_frame='BARY')], cfg)
    # A wide test channel straddling a reference flag gap cannot quietly fill it.
    with pytest.raises(ValueError, match='overlap'):
        pa.common_selection([channel(0, .5), channel(2, 2.5)], [channel(0, 1.5, 3)], cfg)


class MockTable:
    def __init__(self, db, reads):
        self.db, self.reads, self.rows, self.path = db, reads, None, None
    def open(self, path, nomodify=True):
        assert nomodify is True
        self.path = path
        self.rows = self.db[path]
    def close(self):
        self.reads.append(('close', self.path))
    def colnames(self):
        return list(self.rows[0]) if self.rows else []
    def nrows(self):
        return len(self.rows)
    def iscelldefined(self, col, row):
        return self.rows[row].get(col) is not None
    def getcell(self, col, row):
        result = self.rows[row][col]
        if result is None:
            raise RuntimeError('undefined')
        return result
    def getcol(self, col, startrow=0, nrow=-1):
        nrow = len(self.rows)-startrow if nrow < 0 else nrow
        if self.path.endswith('.ms'):
            self.reads.append((col, startrow, nrow))
        values = [self.getcell(col, i) for i in range(startrow, startrow+nrow)]
        return np.stack(values, axis=-1)
    def getcolkeyword(self, col, key):
        return {'Ref': 'J2000'}
    def query(self, expression):
        import re
        fid, ddid = [int(x) for x in re.findall(r'== (\d+)', expression)]
        result = MockTable(self.db, self.reads)
        result.path = self.path
        result.rows = [r for r in self.rows if r['FIELD_ID'] == fid and r['DATA_DESC_ID'] == ddid and r['ANTENNA1'] != r['ANTENNA2']]
        return result


def fixture_tables(cfg, nchans=(16,), scans=(2, 8)):
    db = {}
    reads = []
    for field in cfg['fields']:
        for role in ('reference', 'test'):
            path = field[role + '_ms']
            db[path + '/FIELD'] = [dict(NAME=field['name'], NUM_POLY=0, PHASE_DIR=np.array([[1.2], [-.8]])),
                                   dict(NAME='unrelated', NUM_POLY=0, PHASE_DIR=np.array([[2.0], [-.7]]))]
            db[path + '/SPECTRAL_WINDOW'] = [dict(CHAN_FREQ=1e9+100*np.arange(n)+j*10000,
                CHAN_WIDTH=np.full(n, 100), MEAS_FREQ_REF=5) for j, n in enumerate(nchans)]
            db[path + '/POLARIZATION'] = [dict(CORR_TYPE=np.array([9, 12]))]
            db[path + '/DATA_DESCRIPTION'] = [dict(SPECTRAL_WINDOW_ID=j, POLARIZATION_ID=0) for j in range(len(nchans))]
            rows = []
            own_scans = scans if role == 'reference' else tuple(s+10 for s in scans)
            for ddid, nchan in enumerate(nchans):
                for scan in own_scans:
                    for baseline in range(3):
                        rows.append(dict(FIELD_ID=0, DATA_DESC_ID=ddid, SCAN_NUMBER=scan, ANTENNA1=baseline,
                            ANTENNA2=baseline+1, FLAG=np.zeros((2, nchan), bool), FLAG_ROW=False,
                            DATA=np.ones((2, nchan), complex), WEIGHT=np.array([1., 2.]),
                            WEIGHT_SPECTRUM=np.ones((2, nchan)), TIME=1e9+scan, UVW=np.array([20., 30., 40.])))
            # Ignore wrong-field and autocorrelation rows even if their scan IDs differ.
            wrong = copy.deepcopy(rows[0]); wrong['FIELD_ID'] = 1; wrong['SCAN_NUMBER'] = 999
            auto = copy.deepcopy(rows[0]); auto['ANTENNA2'] = auto['ANTENNA1']; auto['SCAN_NUMBER'] = 998
            db[path] = rows + [wrong, auto]
    return db, reads, lambda: MockTable(db, reads)


def mock_task(vis='', imagename='', stokes='I', specmode='mfs',
              cell='1arcsec', imsize=100, niter=0, gain=.1, deconvolver='hogbom',
              threshold=0, gridder='standard', wprojplanes=1, weighting='natural', robust=.5,
              pblimit=.2, field='', scan='', spw='', datacolumn='corrected', antenna='',
              phasecenter='', selectdata=True, restart=True, savemodel='none', parallel=False,
              interactive=False, restoration=True, calcpsf=True, calcres=True, restoringbeam='',
              pbcor=False):
    pass


def test_survey_multiple_variable_spws_weights_absence_and_chunking(tmp_path):
    cfg = config(tmp_path)
    cfg['chunk_rows'] = 2
    db, reads, factory = fixture_tables(cfg, nchans=(8, 16))
    path = cfg['fields'][0]['reference_ms']
    for row in db[path][:-2]:
        row['WEIGHT_SPECTRUM'][:, 1] = 0 # spectrum overrides positive WEIGHT
        row['FLAG'][:, 3] = True
    db[path][0]['WEIGHT_SPECTRUM'] = None # explicit undefined -> row weight
    db[path][1]['DATA'] = None # absent required cell -> unusable, denominator retained
    survey = pa.survey_ms(path, 'gain', cfg, factory)
    assert sorted(survey['by_scan']) == [2, 8]
    assert survey['weight_sources'] == dict(WEIGHT=1, WEIGHT_SPECTRUM=10, absent=1)
    chans = pa.channel_records(survey)
    assert len(chans) == 24
    assert not next(c for c in chans if c['spw'] == 0 and c['channel'] == 1)['eligible']
    assert not next(c for c in chans if c['channel'] == 3)['eligible']
    absent = next(c for c in chans if c['spw'] == 0 and c['channel'] == 0)
    assert absent['denominator'] == 12
    assert absent['absent_samples'] == 2
    assert absent['usable_fraction'] == pytest.approx(10/12)
    assert all(entry[2] <= 2 for entry in reads if entry[0] != 'close')


def test_plan_separate_field_and_independent_scan_matrix_and_manifest(tmp_path):
    cfg = config(tmp_path)
    db, reads, factory = fixture_tables(cfg)
    # Full-band request retains a channel whose first test scan has no valid data.
    path = cfg['fields'][1]['test_ms']
    for row in db[path]:
        if row['SCAN_NUMBER'] == 12:
            row['FLAG'][:, 5] = True
    # Lower aggregate threshold so it remains eligible across the pair.
    cfg['occupancy_threshold'] = .5
    plan = pa.build_plan(cfg, factory, mock_task, '6.6.1')
    assert len(plan['images']) == 20 # 2 roles x 2 fields x (2 scans + 3 bands)
    assert len(plan['native_diagnostics']) == 4
    for product in plan['images']:
        c = product['contract']
        params = c['tclean_parameters']
        assert c['data_column'] == 'DATA' and params['datacolumn'] == 'data'
        assert params['stokes'] == 'I' and params['specmode'] == 'mfs'
        assert params['field'] == '0' and params['antenna'] == '*&*'
        assert params['savemodel'] == 'none' and params['restart'] is False
        assert c['input_phase_centre']['frame'] == 'J2000'
        assert c['imaging_phase_centre']['frame'] == 'J2000'
        assert len(c['requested_channels']) > 0
        assert c['ms_identity']['selected_content_sha256']
        assert c['paired_inputs']['reference']['ms_identity']['path'].endswith(c['field_name'] + '.ms')
        assert c['common_frequency_selection']['reference_channel_ids']
        assert c['common_frequency_selection']['test_channel_ids']
        assert c['selected_correlations'][0]['selected'] == ['XX', 'YY']
        assert c['occupancy_policy']['absent_data']
        assert c['common_frequency_selection']['coverage_fractions']['reference'] == 1
        assert c['field_kind'] in ('gain_calibrator', 'target')
        expected = [2, 8] if c['ms_role'] == 'reference' else [12, 18]
        assert c['scan_ids'] == expected if c['product'] != 'fullband' else c['scan_ids'][0] in expected
    p = next(p for p in plan['images'] if p['contract']['field_name'] == 'target' and p['contract']['ms_role'] == 'test' and p['contract']['scan_ids'] == [12])
    assert p['contract']['per_scan_support_incomplete']
    assert p['contract']['selection']['channel_count'] == 16
    assert p['contract']['selection']['contributing_channel_count'] == 15
    assert [0, 5] not in p['contract']['selection']['contributed_channel_ids']
    assert len({p['contract']['output_paths']['manifest'] for p in plan['images']}) == 20
    # All parameters, including unmodified task defaults, are captured.
    assert plan['images'][0]['contract']['tclean_parameters']['pbcor'] is False


def test_manifest_execution_reuse_and_stale_detection(tmp_path, monkeypatch):
    cfg = config(tmp_path)
    _, _, factory = fixture_tables(cfg)
    plan = pa.build_plan(cfg, factory, mock_task, '6.6.1')
    calls = []
    def task(**params):
        calls.append(params)
        Path(params['imagename'] + '.image').mkdir()
    def export(**params):
        Path(params['fitsimage']).write_text('mock fits')
    monkeypatch.setattr(pa, 'fits_metadata', lambda path: ({'ra_deg': 1.2, 'dec_deg': -.8}, {'major_arcsec': 8}))
    pa.execute_plan(plan, task, export)
    assert len(calls) == 20
    for p in plan['images']:
        saved = json.loads(Path(p['contract']['output_paths']['manifest']).read_text())
        assert saved['processing_status'] == 'complete'
        assert saved['fits_wcs_centre'] and saved['restoring_beam']
        assert saved['contribution_status']
        assert pa.output_state(p) == 'reuse'
    again = pa.build_plan(cfg, factory, mock_task, '6.6.1')
    pa.execute_plan(again, task, export)
    assert len(calls) == 20
    assert all(p['processing_status'] == 'complete' and p['reused'] for p in again['images'])
    cfg['tclean']['robust'] = .1
    stale = pa.build_plan(cfg, factory, mock_task, '6.6.1')
    with pytest.raises(RuntimeError, match='Stale'):
        pa.execute_plan(stale, task, export)
    assert len(calls) == 20 # preflight catches everything before tclean
    changed = copy.deepcopy(plan['images'][0]); changed['contract']['ms_identity']['selected_content_sha256'] = 'changed'
    changed['fingerprint'] = pa.digest(changed['contract'])
    with pytest.raises(RuntimeError, match='Stale'):
        pa.output_state(changed)
    Path(plan['images'][0]['contract']['output_paths']['manifest']).unlink()
    with pytest.raises(RuntimeError, match='Stale'):
        pa.output_state(plan['images'][0])


def test_dry_run_and_task_failure(tmp_path, monkeypatch):
    cfg = config(tmp_path)
    _, _, factory = fixture_tables(cfg)
    cfg['dry_run'] = True
    plan = pa.build_plan(cfg, factory, mock_task, '6.6.1')
    def forbidden(**kwargs):
        pytest.fail('dry run must not image')
    pa.execute_plan(plan, forbidden, forbidden)
    assert Path(cfg['output_dir'], 'selection_plan.json').is_file()
    assert not list(Path(cfg['output_dir']).glob('images/*.manifest.json'))
    cfg['dry_run'] = False
    plan = pa.build_plan(cfg, factory, mock_task, '6.6.1')
    def failed(**params):
        raise RuntimeError('CASA failure')
    with pytest.raises(RuntimeError, match='Paired imaging failed'):
        pa.execute_plan(plan, failed, forbidden)
    assert all(p['processing_status'] == 'failed' for p in plan['images'])
    assert 'CASA failure' in json.loads(Path(plan['images'][0]['contract']['output_paths']['manifest']).read_text())['error']


def test_pipeline_pair_launch_and_legacy_default(tmp_path, monkeypatch):
    pair = config(tmp_path)
    cfg = _dict_to_dataclass(dict(project_name='test', casa=dict(paired_astrometry=pair)))
    calls = []
    monkeypatch.setenv('MCI_TCLEAN_MSFILE', 'stale.ms')
    monkeypatch.setattr(step3_calibrate_image, '_run_cmd', lambda cmd, **kwargs: calls.append((cmd, kwargs)))
    step3_calibrate_image.run(cfg)
    assert len(calls) == 1
    command, kw = calls[0]
    assert command[-1].endswith('casa_image_batch.py')
    assert 'MCI_TCLEAN_MSFILE' not in kw['env']
    assert len(kw['inputs']) == 4
    assert json.loads(kw['env']['MCI_PAIRED_ASTROMETRY_JSON'])['fields'] == pair['fields']
    cfg.extra['force_calibrate'] = True
    with pytest.raises(ValueError, match='force_calibrate'):
        step3_calibrate_image.run(cfg)
    assert _dict_to_dataclass(dict(project_name='legacy')).casa.paired_astrometry == {}


def test_fits_wcs_centre_and_beam(tmp_path):
    from astropy.io import fits
    from astropy.wcs import WCS
    w = WCS(naxis=2)
    w.wcs.crpix = [3., 3.]
    w.wcs.cdelt = [-1/3600, 1/3600]
    w.wcs.crval = [180., -45.]
    w.wcs.ctype = ['RA---SIN', 'DEC--SIN']
    h = w.to_header(); h['BMAJ'] = 8/3600; h['BMIN'] = 7/3600; h['BPA'] = 12
    path = tmp_path/'test.fits'
    fits.writeto(path, np.zeros((5, 5)), h)
    centre, beam = pa.fits_metadata(str(path))
    assert centre['ra_deg'] == pytest.approx(180.)
    assert centre['dec_deg'] == pytest.approx(-45.)
    assert beam == pytest.approx(dict(major_arcsec=8, minor_arcsec=7, pa_deg=12))


@pytest.mark.parametrize('update', [dict(occupancy_threshold=0), dict(min_common_coverage=1.1),
    dict(min_channel_overlap=float('nan')), dict(chunk_rows=-1), dict(unknown=True),
    dict(tclean={'stokes': 'pseudoI'}), dict(tclean={'deconvolver': 'mtmfs'})])
def test_config_rejects_invalid_options(update):
    with pytest.raises(ValueError):
        pa.validate_config(update)


def test_supported_casa_and_sample_configs(tmp_path):
    import yaml
    cfg = config(tmp_path)
    _, _, factory = fixture_tables(cfg)
    with pytest.raises(RuntimeError, match='CASA 6'):
        pa.build_plan(cfg, factory, mock_task, '5.8.0')
    root = Path(__file__).resolve().parents[1]
    for name in ('l_band.yaml', 's4.yaml'):
        sample = _dict_to_dataclass(yaml.safe_load((root/'configs/paired_astrometry'/name).read_text()))
        pair = sample.casa.paired_astrometry
        assert pair['enabled']
        assert [f['name'] for f in pair['fields']] == ['J1619-8418', 'J2147-8132']
        assert all('PLACEHOLDER' in f['reference_ms'] for f in pair['fields'])
        assert pair['datacolumn'] == 'data' and pair['min_common_coverage'] == .9


def test_empty_spw_denominator_frequency_metadata_and_bits(tmp_path):
    cfg = config(tmp_path)
    db, _, factory = fixture_tables(cfg, nchans=(8, 4))
    path = cfg['fields'][0]['reference_ms']
    db[path] = [row for row in db[path] if row['DATA_DESC_ID'] == 0]
    db[path+'/SPECTRAL_WINDOW'][0]['CHAN_WIDTH'][2] = 0
    for row in db[path]:
        row['BITFLAG'] = np.full((2, 8), 16, dtype=int)
        row['BITFLAG_ROW'] = 0
    survey = pa.survey_ms(path, 'gain', cfg, factory)
    chans = pa.channel_records(survey)
    empty = [c for c in chans if c['spw'] == 1]
    assert all(c['denominator'] == 0 and not c['eligible'] for c in empty)
    assert not next(c for c in chans if c['spw'] == 0 and c['channel'] == 2)['eligible']
    assert survey['raw_flag_bits']['BITFLAG'] == {'16': 96}
    assert all(c['eligible'] for c in chans if c['spw'] == 0 and c['channel'] != 2)
    # Raw bits are not silently substituted for CASA flags.
    assert all(c['flagged'] == 0 for c in chans)


def test_absent_flags_do_not_invent_flag_bits(tmp_path):
    cfg = config(tmp_path)
    db, _, factory = fixture_tables(cfg)
    path = cfg['fields'][0]['reference_ms']
    db[path][0]['FLAG'] = None
    survey = pa.survey_ms(path, 'gain', cfg, factory)
    records = pa.channel_records(survey)
    assert all(c['absent_flag_samples'] == 2 and c['absent_samples'] == 2 for c in records)
    assert all(c['flagged'] == 0 for c in records)
    assert all(c['joint_usable'] == 10 for c in records)


def test_overlap_uses_reference_quartiles_with_different_channel_grids(tmp_path):
    cfg = config(tmp_path)
    db, _, factory = fixture_tables(cfg)
    # Reference 16 x 100 Hz; test 8 x 200 Hz over identical frequency support.
    for field in cfg['fields']:
        path = field['test_ms']
        spw = db[path+'/SPECTRAL_WINDOW'][0]
        spw['CHAN_FREQ'] = 1e9+50+200*np.arange(8)
        spw['CHAN_WIDTH'] = np.full(8, 200)
        for row in db[path]:
            row['FLAG'] = np.zeros((2, 8), bool)
            row['DATA'] = np.ones((2, 8), complex)
            row['WEIGHT_SPECTRUM'] = np.ones((2, 8))
    plan = pa.build_plan(cfg, factory, mock_task, '6.6.1')
    for p in plan['images']:
        c = p['contract']
        assert c['common_frequency_selection']['exact_frequency_support']
        expected = 16 if c['product'] == 'fullband' else 4
        assert c['selection']['channel_count'] == (expected if c['ms_role'] == 'reference' else expected//2)
    middle = [p['contract']['tclean_parameters']['spw'] for p in plan['images'] if p['contract']['product'] == 'middle']
    assert middle == ['0:6~9', '0:3~4']*2


def test_selection_failure_preserves_all_native_diagnostics(tmp_path, monkeypatch):
    import sys
    from types import SimpleNamespace
    cfg = config(tmp_path)
    db, _, factory = fixture_tables(cfg)
    db[cfg['fields'][0]['test_ms']+'/SPECTRAL_WINDOW'][0]['CHAN_FREQ'] += 1e8
    monkeypatch.setitem(sys.modules, 'casatasks', SimpleNamespace(version_string=lambda: '6.6.1', tclean=mock_task, exportfits=None))
    monkeypatch.setitem(sys.modules, 'casatools', SimpleNamespace(table=factory))
    with pytest.raises(pa.SelectionPlanningError, match='gain product low'):
        pa.run(cfg)
    diagnostics = json.loads(Path(cfg['output_dir'], 'native_channel_diagnostics.json').read_text())
    assert len(diagnostics) == 4
    failure = json.loads(Path(cfg['output_dir'], 'planning_failure.json').read_text())
    assert failure['processing_status'] == 'failed' and failure['phase'] == 'selection'
    assert not Path(cfg['output_dir'], 'images').exists()


@pytest.mark.parametrize('target_only', [False, True])
def test_compatibility_entrypoint_dispatches_paired_config(tmp_path, monkeypatch, target_only):
    import os
    import sys
    import meerkat_corr_imaging
    package = Path(meerkat_corr_imaging.__file__).parent
    script = package.parent.parent/'scripts/tclean_two_bands.py'
    cfg = target_only_config(tmp_path) if target_only else config(tmp_path)
    cfg['dry_run'] = True
    _, _, factory = fixture_tables(cfg)
    from types import SimpleNamespace
    monkeypatch.setitem(sys.modules, 'casatasks', SimpleNamespace(version_string=lambda: '6.6.1', tclean=mock_task, exportfits=lambda **kw: pytest.fail('export')))
    monkeypatch.setitem(sys.modules, 'casatools', SimpleNamespace(table=factory))
    monkeypatch.setenv('MCI_CASA_SCRIPT_DIR', str(package))
    monkeypatch.setenv('MCI_PAIRED_ASTROMETRY_JSON', json.dumps(cfg))
    monkeypatch.setenv('MCI_TCLEAN_MSFILE', 'ignored_legacy.ms')
    monkeypatch.setattr(sys, 'argv', [str(script)])
    exec(compile(script.read_text(), str(script), 'exec'), {'__name__': '__main__', '__file__': str(script)})
    assert len(json.loads(Path(cfg['output_dir'], 'selection_plan.json').read_text())['images']) == (10 if target_only else 20)


def test_interval_operations_randomized_against_unit_grid():
    rng = np.random.default_rng(19)
    for _ in range(50):
        a = [list(map(int, sorted(rng.choice(30, 2, replace=False)))) for _ in range(8)]
        b = [list(map(int, sorted(rng.choice(30, 2, replace=False)))) for _ in range(8)]
        def cells(ranges):
            return {i for lo, hi in ranges for i in range(int(lo), int(hi))}
        assert cells(pa.intersection(a, b)) == cells(a) & cells(b)
        assert cells(pa.difference(a, b)) == cells(a) - cells(b)
        cs = [channel(i, i+.5) for i in range(30)]
        kept = pa.overlap_filter(cs, b, .9)
        assert {c['channel'] for c in kept} == cells(b)


def test_content_identity_independent_of_chunk_size_and_tracks_visibility_changes(tmp_path):
    cfg = config(tmp_path)
    db, _, factory = fixture_tables(cfg)
    path = cfg['fields'][0]['reference_ms']
    cfg['chunk_rows'] = 1
    first = pa.survey_ms(path, 'gain', cfg, factory)
    cfg['chunk_rows'] = 3
    second = pa.survey_ms(path, 'gain', cfg, factory)
    assert first['identity'] == second['identity']
    db[path][0]['DATA'][0, 0] = 20+1j
    third = pa.survey_ms(path, 'gain', cfg, factory)
    assert first['identity']['selected_content_sha256'] != third['identity']['selected_content_sha256']


def test_external_mask_content_changes_fingerprint(tmp_path):
    cfg = config(tmp_path)
    _, _, factory = fixture_tables(cfg)
    mask = tmp_path/'mask.crtf'; mask.write_text('region one')
    # Extend the mocked task signature with CASA's mask default.
    import inspect
    def task(**params):
        pass
    sig = inspect.signature(mock_task)
    task.__signature__ = sig.replace(parameters=list(sig.parameters.values()) +
        [inspect.Parameter('mask', inspect.Parameter.POSITIONAL_OR_KEYWORD, default='')])
    cfg['tclean']['mask'] = str(mask)
    first = pa.build_plan(cfg, factory, task, '6.6.1')
    mask.write_text('region two')
    second = pa.build_plan(cfg, factory, task, '6.6.1')
    assert first['images'][0]['fingerprint'] != second['images'][0]['fingerprint']


def target_only_config(tmp_path):
    cfg = config(tmp_path)
    cfg['experiment_mode'] = 'target_only'
    cfg['fields'] = [f for f in cfg['fields'] if f['kind'] == 'target']
    return pa.validate_config(cfg)


def test_experiment_mode_default_and_explicit_target_only(tmp_path):
    cfg = config(tmp_path)
    assert cfg['experiment_mode'] == 'calibrator_and_target'
    raw = copy.deepcopy(cfg); raw.pop('experiment_mode')
    assert pa.validate_config(raw)['experiment_mode'] == 'calibrator_and_target'
    raw['fields'] = [raw['fields'][1]]
    with pytest.raises(ValueError, match='explicit experiment_mode'):
        pa.validate_config(raw)
    raw['experiment_mode'] = 'target_only'
    assert pa.validate_config(raw)['fields'][0]['kind'] == 'target'
    loaded = _dict_to_dataclass(dict(project_name='target', casa=dict(paired_astrometry=raw)))
    assert loaded.casa.paired_astrometry['experiment_mode'] == 'target_only'


@pytest.mark.parametrize('mode,kinds', [
    ('automatic', ['target']), (None, ['target']),
    ('target_only', ['gain_calibrator']),
    ('target_only', ['target', 'target']), ('target_only', []),
    ('calibrator_and_target', ['target']), ('calibrator_and_target', ['gain_calibrator']),
])
def test_reject_conflicting_experiment_modes(tmp_path, mode, kinds, monkeypatch):
    cfg = config(tmp_path)
    originals = {f['kind']: f for f in cfg['fields']}
    cfg['experiment_mode'] = mode
    cfg['fields'] = [copy.deepcopy(originals[kind]) for kind in kinds]
    # Invalid mode/field combinations are rejected before resolving MS paths.
    monkeypatch.setattr(pa.os.path, 'realpath', lambda *args, **kw: pytest.fail('Conflicting config accessed an MS path'))
    with pytest.raises(ValueError):
        pa.validate_config(cfg)


def test_target_only_requires_both_paths_and_keeps_weight_policy(tmp_path):
    cfg = target_only_config(tmp_path)
    invalid = copy.deepcopy(cfg); invalid['fields'][0].pop('reference_ms')
    with pytest.raises(ValueError, match='reference_ms'):
        pa.validate_config(invalid)
    db, _, factory = fixture_tables(cfg)
    path = cfg['fields'][0]['reference_ms']
    for row in db[path]:
        row['WEIGHT'][:] = 0
        row['WEIGHT_SPECTRUM'][:] = 0
    with pytest.raises(pa.SelectionPlanningError, match='fewer than four'):
        pa.build_plan(cfg, factory, mock_task, '6.6.5')
    # Even with positive WEIGHT, a defined zero spectrum is still unusable.
    for row in db[path]:
        row['WEIGHT'][:] = 1
    with pytest.raises(pa.SelectionPlanningError, match='fewer than four'):
        pa.build_plan(cfg, factory, mock_task, '6.6.5')


def test_target_only_matrix_frequency_coverage_and_provenance(tmp_path):
    cfg = target_only_config(tmp_path)
    db, _, factory = fixture_tables(cfg)
    # Preserve a real irregular gap in both roles, without filling it in.
    for field in cfg['fields']:
        for role in ('reference', 'test'):
            for row in db[field[role+'_ms']]:
                row['FLAG'][:, 7] = True
    # Different scan count and IDs in test: remove its second scan.
    path = cfg['fields'][0]['test_ms']
    db[path] = [r for r in db[path] if r['SCAN_NUMBER'] != 18]
    plan = pa.build_plan(cfg, factory, mock_task, '6.6.5')
    assert len(plan['images']) == 6+2+1
    assert len(plan['native_diagnostics']) == 2
    assert plan['experiment_mode'] == 'target_only'
    assert plan['calibrator_imaging_qa_status'] == 'not_performed_target_only'
    assert cfg['occupancy_threshold'] == .8 and cfg['min_common_coverage'] == .9
    assert 'phase transfer' in plan['interpretation_limitations']
    for p in plan['images']:
        c = p['contract']
        assert c['field_kind'] == 'target'
        assert c['experiment_mode'] == 'target_only'
        assert c['calibrator_imaging_qa_status'] == 'not_performed_target_only'
        assert 'calibrator-image check' in c['interpretation_limitations']
        assert c['data_column'] == 'DATA'
        assert c['tclean_parameters']['stokes'] == 'I' and c['tclean_parameters']['specmode'] == 'mfs'
        assert c['common_frequency_selection']['coverage_fractions'] == dict(requested_reference=1, reference=1, test=1)
        assert [0, 7] not in c['selection']['contributed_channel_ids']
        if c['product'] == 'fullband':
            assert c['scan_ids'][0] in ([2, 8] if c['ms_role'] == 'reference' else [12])
            assert 'independent scans' in c['actual_support_comparability']
        else:
            assert c['scan_ids'] == ([2, 8] if c['ms_role'] == 'reference' else [12])
    for role in ('reference', 'test'):
        assert {p['contract']['product'] for p in plan['images'] if p['contract']['ms_role'] == role} == {'low', 'middle', 'high', 'fullband'}
    # Keep the same quartile helper and integer rounding as combined mode.
    native = pa.quartiles([c for c in plan['native_diagnostics'][0]['channels'] if c['eligible']])
    for band in ('low', 'middle', 'high'):
        c = next(p['contract'] for p in plan['images'] if p['contract']['ms_role'] == 'reference' and p['contract']['product'] == band)
        assert c['common_frequency_selection']['reference_channel_ids'] == [list(pa.channel_key(ch)) for ch in native[band]]


def test_target_only_pipeline_never_accesses_calibrator_paths(tmp_path, monkeypatch):
    import builtins
    from meerkat_corr_imaging.config import Target
    cfg = target_only_config(tmp_path)
    cfg['dry_run'] = True
    _, _, factory = fixture_tables(cfg) # no calibrator tables exist in the mock
    forbidden_paths = [str(tmp_path / (role+'_gain.ms')) for role in ('reference', 'test')]
    cfg['fields'].insert(0, dict(name='gain', kind='gain_calibrator',
        reference_ms=forbidden_paths[0], test_ms=forbidden_paths[1]))
    master = _dict_to_dataclass(dict(project_name='target', casa=dict(paired_astrometry=cfg)))
    # Legacy MS lists and inherited overrides must not enter the paired matrix.
    master.reference = Target(name='legacy-cal', ms_paths=[forbidden_paths[0]])
    master.tests = [Target(name='legacy-test-cal', ms_paths=[forbidden_paths[1]])]
    monkeypatch.setenv('MCI_TCLEAN_MSFILE', forbidden_paths[0])
    accessed = []
    def guard(path):
        path = str(path)
        assert not any(path == p or path.startswith(p + '/') for p in forbidden_paths), path
        accessed.append(path)
    original_open, original_stat, original_realpath = builtins.open, pa.os.stat, pa.os.path.realpath
    def checked_open(path, *args, **kw):
        guard(path); return original_open(path, *args, **kw)
    def checked_stat(path, *args, **kw):
        guard(path); return original_stat(path, *args, **kw)
    def checked_realpath(path, *args, **kw):
        guard(path); return original_realpath(path, *args, **kw)
    monkeypatch.setattr(builtins, 'open', checked_open)
    monkeypatch.setattr(pa.os, 'stat', checked_stat)
    monkeypatch.setattr(pa.os.path, 'realpath', checked_realpath)
    def forbidden(**kw):
        pytest.fail('Selection-only target run must not calibrate/image/export')
    monkeypatch.setattr(step3_calibrate_image, 'run_calibration', forbidden)
    calls = []
    def launch(command, env=None, inputs=()):
        calls.append(command)
        assert len(inputs) == 2 and 'MCI_TCLEAN_MSFILE' not in env
        pair = json.loads(env['MCI_PAIRED_ASTROMETRY_JSON'])
        plan = pa.build_plan(pair, factory, mock_task, '6.6.5')
        pa.execute_plan(plan, forbidden, forbidden)
    monkeypatch.setattr(step3_calibrate_image, '_run_cmd', launch)
    step3_calibrate_image.run(master)
    assert len(calls) == 1
    assert accessed and all('gain.ms' not in p for p in accessed)


def test_default_mode_does_not_fall_back_when_calibrator_has_zero_weights(tmp_path):
    cfg = config(tmp_path)
    db, _, factory = fixture_tables(cfg)
    for role in ('reference', 'test'):
        for row in db[cfg['fields'][0][role+'_ms']]:
            row['WEIGHT'][:] = 0; row['WEIGHT_SPECTRUM'][:] = 0
    with pytest.raises(pa.SelectionPlanningError, match='Field gain') as error:
        pa.build_plan(cfg, factory, mock_task, '6.6.5')
    assert len(error.value.diagnostics) == 4
    assert all(d['experiment_mode'] == 'calibrator_and_target' for d in error.value.diagnostics)


def test_mode_changes_reject_existing_target_outputs(tmp_path, monkeypatch):
    cfg = config(tmp_path)
    _, _, factory = fixture_tables(cfg)
    combined = pa.build_plan(cfg, factory, mock_task, '6.6.5')
    def task(**params):
        Path(params['imagename']+'.image').mkdir()
    def export(**params):
        Path(params['fitsimage']).write_text('mock fits')
    monkeypatch.setattr(pa, 'fits_metadata', lambda path: ({'ra_deg': 1}, {'major_arcsec': 8}))
    pa.execute_plan(combined, task, export)
    cfg['experiment_mode'] = 'target_only'; cfg['fields'] = [cfg['fields'][1]]
    target = pa.build_plan(cfg, factory, mock_task, '6.6.5')
    with pytest.raises(RuntimeError, match='Stale'):
        pa.execute_plan(target, lambda **kw: pytest.fail('stale output imaged'), export)
    # Complete target-only manifests carry mode/omitted-QA information too.
    cfg['output_dir'] = str(tmp_path/'target-output')
    target = pa.build_plan(cfg, factory, mock_task, '6.6.5')
    pa.execute_plan(target, task, export)
    for p in target['images']:
        saved = json.loads(Path(p['contract']['output_paths']['manifest']).read_text())
        assert saved['processing_status'] == 'complete'
        assert saved['contract']['experiment_mode'] == 'target_only'
        assert saved['contract']['calibrator_imaging_qa_status'] == 'not_performed_target_only'
        assert saved['fits_wcs_centre'] and saved['restoring_beam']
    cfg['experiment_mode'] = 'calibrator_and_target'; cfg['fields'] = combined['config']['fields']
    reverse = pa.build_plan(cfg, factory, mock_task, '6.6.5')
    p = next(p for p in reverse['images'] if p['contract']['field_kind'] == 'target')
    with pytest.raises(RuntimeError, match='Stale'):
        pa.output_state(p)


def test_target_only_samples_separate_filename_tags_from_field_names():
    import yaml
    root = Path(__file__).resolve().parents[1]
    for name in ('l_band_target_only.yaml', 's4_target_only.yaml'):
        master = _dict_to_dataclass(yaml.safe_load((root/'configs/paired_astrometry'/name).read_text()))
        cfg = master.casa.paired_astrometry
        assert cfg['experiment_mode'] == 'target_only' and cfg['dry_run']
        assert cfg['datacolumn'] == 'data' and len(cfg['fields']) == 1
        field = cfg['fields'][0]
        assert field['name'] == 'J2147-8132' and field['kind'] == 'target'
        for role in ('reference', 'test'):
            path = field[role+'_ms']
            assert '_J2147_SDPflags+cal.' in path and 'J2147-8132' not in path and 'PLACEHOLDER' in path
        assert master.extra['force_calibrate'] is False
        assert '_target_only' in cfg['output_dir']


@pytest.mark.parametrize('ignored', [
    dict(kind='gain_calibrator'),
    dict(kind='gain_calibrator', name='gain', reference_ms='/missing/gain.ms', test_ms='/also-missing/gain.ms'),
    dict(kind='gain_calibrator', name='', reference_ms=None, test_ms='', extra_metadata='ignored'),
])
def test_target_only_ignores_calibrator_entries_before_path_validation(tmp_path, monkeypatch, ignored):
    cfg = target_only_config(tmp_path)
    target = copy.deepcopy(cfg['fields'][0])
    cfg['fields'].insert(0, ignored)
    original = copy.deepcopy(cfg)
    realpath = pa.os.path.realpath
    allowed = {target['reference_ms'], target['test_ms']}
    def guarded(path, *args, **kw):
        assert str(path) in allowed, 'Ignored calibrator path was resolved'
        return realpath(path, *args, **kw)
    monkeypatch.setattr(pa.os.path, 'realpath', guarded)
    normalized = pa.validate_config(cfg)
    assert normalized['fields'] == [target]
    assert pa.validate_config(normalized) == normalized
    assert cfg == original # caller's original entries remain intact


def test_default_mode_still_validates_calibrator_entries(tmp_path):
    cfg = config(tmp_path)
    cfg['fields'][0].pop('reference_ms')
    with pytest.raises(ValueError, match='reference_ms'):
        pa.validate_config(cfg)
