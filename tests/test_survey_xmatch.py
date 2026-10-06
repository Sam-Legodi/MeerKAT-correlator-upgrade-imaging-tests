"""Synthetic fixtures test ingestion/associations; no downloaded catalogues in git."""
from copy import deepcopy
import json
from pathlib import Path

import astropy.units as u
from astropy.io import fits
from astropy.table import MaskedColumn, Table
import numpy as np
import pytest

from meerkat_corr_imaging.audit import StepInputsFailed, run_step_with_audit
from meerkat_corr_imaging.config import Config, PathsCfg, Target, load_config
from meerkat_corr_imaging.image_pipeline import build_image_pipeline_plan, wire_image_pipeline
from meerkat_corr_imaging.steps import step6_xmatch
from meerkat_corr_imaging.survey_xmatch import (
    CatalogueSpec, SurveyJob, execute_survey_job, frequency_provenance, resolve_survey_jobs,
)
from meerkat_corr_imaging.xmatch_pybdsf import (
    cross_match_catalogues, parse_angle, prepare_catalog,
)


def _fixtures(tmp_path, *, xml=False):
    # One MeerKAT source, two Gaussian components of a common RACS parent.
    a = Table({'RA': [10.], 'DEC': [-20.], 'Source_id': [7]})
    b = Table({'RAJ2000': [10., 10.0001, 50.], 'DEJ2000': [-20., -20., 0.],
               'GID': ['a', 'b', 'c'], 'ID': ['parent', 'parent', 'other'],
               'Ftot': [20., 20., 40.]})
    b['RAJ2000'].unit = b['DEJ2000'].unit = u.deg
    b['Ftot'].unit = u.mJy
    p1, p2 = tmp_path / 'a.fits', tmp_path / ('b.xml' if xml else 'b.fits')
    a.write(p1)
    if xml:
        b.write(p2, format='votable')
    else:
        # Actual ASCII TableHDU, not a binary table disguised by its extension.
        columns = [fits.Column(name=c, format='A8' if b[c].dtype.kind == 'U' else 'E20.10',
                               unit=str(b[c].unit) if b[c].unit else None, array=b[c])
                   for c in b.colnames]
        fits.HDUList([fits.PrimaryHDU(), fits.TableHDU.from_columns(columns)]).writeto(p2)
    return {'input1': {'profile': 'pybdsf_source', 'path': str(p1)},
            'input2': {'profile': 'racs_low_dr1_gausscut', 'path': str(p2)},
            'output': str(tmp_path / 'matches.fits'), 'radius': '2 arcsec'}


@pytest.mark.parametrize('xml', [False, True], ids=['synthetic-FITS-ASCII', 'synthetic-VOTable'])
def test_survey_ingestion_components_and_provenance(tmp_path, xml):
    raw = _fixtures(tmp_path, xml=xml)
    report = execute_survey_job(raw)
    output = tmp_path / 'Sky-CrossMatches' / 'matches.fits'
    t = Table.read(output)
    assert len(t) == 2
    assert list(t['GID_2']) == ['a', 'b']
    assert list(t['ID_2']) == ['parent', 'parent']
    assert list(t['Ftot_2']) == [20., 20.]  # no summation/deduplication
    assert t['Ftot_2'].unit == u.mJy
    assert t['sep_arcsec'].unit == u.arcsec
    assert list(t['input_row_1']) == [0, 0]
    assert list(t['input_row_2']) == [0, 1]
    assert list(t['candidate_count_1']) == [2, 2]
    assert all(t['ambiguous'])
    assert 'total_flux_ratio' not in t.colnames  # survey analysis excluded
    assert report['catalogues'][0] == {'input': 1, 'valid': 1, 'invalid': 0, 'matched': 1,
                                      'unmatched': 0, 'ambiguous': 1, 'grain': 'source', 'distinct_source_ids': 1,
                                      'valid_source_ids': 1, 'matched_source_ids': 1,
                                      'unmatched_source_ids': 0, 'ambiguous_source_ids': 1}
    assert report['catalogues'][1]['distinct_source_ids'] == 2
    assert report['frequencies'][1]['frequency_hz'] == 887.5e6
    assert json.loads(output.with_suffix('.diagnostics.json').read_text()) == report
    assert t.meta['XM2_GRAIN'] == 'component'
    assert 'DR1' in ''.join(t.meta['HISTORY'])
    before = output.read_bytes()
    with pytest.raises(FileExistsError):
        execute_survey_job(raw)
    assert output.read_bytes() == before


@pytest.mark.parametrize('radius', ['1.0', '0 arcsec', '-1 deg', 'nan arcsec', 'inf deg'])
def test_invalid_radii(radius):
    with pytest.raises(ValueError):
        parse_angle(radius)


def test_coordinate_units_masks_ranges_and_ra_wrap(tmp_path):
    a = Table()
    a['RA'] = MaskedColumn([359.9999, 12., np.nan, -1., 361., 0., 10.],
                           mask=[False, True, False, False, False, False, False], unit=u.deg)
    a['DEC'] = [0., 0., 0., 0., 0., 91., 0.]
    p1, p2 = tmp_path / 'a.fits', tmp_path / 'b.fits'
    a.write(p1)
    Table({'RA': [0.0001], 'DEC': [0.]}).write(p2)
    c1, c2 = [prepare_catalog(str(p), 'RA', 'DEC', 'icrs', 'T') for p in (p1, p2)]
    assert list(c1.valid_indices) == [0, 6]
    assert c1.skipped == 5
    t, matches = cross_match_catalogues(c1, c2, parse_angle('1 arcsec'), True, enrich=False)
    assert [(m[0], m[1]) for m in matches] == [(0, 0)]
    assert t['sep_arcsec'][0] == pytest.approx(.72)
    bad = Table({'RA': [1.], 'DEC': [1.]})
    bad['RA'].unit = u.m
    bad.write(tmp_path / 'bad.fits')
    with pytest.raises(u.UnitConversionError):
        prepare_catalog(str(tmp_path / 'bad.fits'), 'RA', 'DEC', 'icrs', 'T')


def test_empty_matches_and_all_invalid(tmp_path):
    raw = _fixtures(tmp_path)
    t = Table.read(raw['input1']['path'])
    t['DEC'] = [100.]
    t.write(raw['input1']['path'], overwrite=True)
    report = execute_survey_job(raw)
    assert report['associations'] == 0
    assert report['catalogues'][0]['valid'] == 0
    assert report['catalogues'][0]['invalid'] == 1
    assert report['coverage']['status'] == 'unknown'
    t = Table.read(tmp_path / 'Sky-CrossMatches' / 'matches.fits')
    assert len(t) == 0 and 'ambiguous' in t.colnames


def test_one_to_one_keeps_legacy_selection_but_flags_ambiguity(tmp_path):
    raw = _fixtures(tmp_path)
    raw['association_mode'] = 'one-to-one'
    report = execute_survey_job(raw)
    t = Table.read(tmp_path / 'Sky-CrossMatches' / 'matches.fits')
    assert len(t) == 1 and t['input_row_2'][0] == 0
    assert t['ambiguous'][0]
    assert report['catalogues'][1]['unmatched'] == 2
    assert report['catalogues'][1]['ambiguous'] == 2


def test_table_selection_required_columns_and_incomplete_files(tmp_path):
    raw = _fixtures(tmp_path)
    p = tmp_path / 'multi.fits'
    with fits.open(raw['input2']['path']) as hs:
        fits.HDUList([fits.PrimaryHDU(), hs[1].copy(), hs[1].copy()]).writeto(p)
    raw['input2']['path'] = str(p)
    # Legacy preparation continues to select the first FITS table.
    assert len(prepare_catalog(str(p), 'RAJ2000', 'DEJ2000', 'icrs', 'T').table) == 3
    spec = CatalogueSpec.from_config(raw['input2'])
    with pytest.raises(ValueError, match='Select an HDU'):
        spec.prepare('T')
    raw['input2']['hdu'] = 2
    assert len(CatalogueSpec.from_config(raw['input2']).prepare('T').table) == 3
    raw['input2']['source_id'] = 'missing'
    with pytest.raises(KeyError, match='missing'):
        CatalogueSpec.from_config(raw['input2']).prepare('T')
    raw['input2']['path'] = str(tmp_path / 'download.crdownload')
    with pytest.raises(ValueError, match='Incomplete'):
        CatalogueSpec.from_config(raw['input2']).prepare('T')


def test_survey_settings_do_not_change_legacy_pairs(tmp_path, monkeypatch):
    raw = _fixtures(tmp_path)
    cfg = Config('mixed', paths=PathsCfg(sky_xmatches_dir=str(tmp_path / 'Sky-CrossMatches')))
    cfg.extra = {'xmatch_pairs': [['a.fits', 'b.fits'], ['c.fits', 'd.fits', 'custom.fits']],
                 'survey_xmatch_jobs': [raw]}
    commands = []
    monkeypatch.setattr(step6_xmatch, '_run', lambda cmd, **kw: commands.append(cmd))
    step6_xmatch.run(cfg)
    assert len(commands) == 2
    for command in commands:
        assert command[command.index('--max-error') + 1] == '1.0 arcsec'
        assert command[command.index('--ra-col-2') + 1] == 'RA'
    assert len(Table.read(tmp_path / 'Sky-CrossMatches' / 'matches.fits')) == 2


def test_missing_requested_survey_is_audited_and_disabled_job_ignored(tmp_path):
    raw = _fixtures(tmp_path)
    raw['input2']['path'] = str(tmp_path / 'missing.fits')
    cfg = Config('audit', paths=PathsCfg(reports_dir=str(tmp_path / 'reports')))
    cfg.extra = {'survey_xmatch_jobs': [raw, {'enabled': False}]}
    with pytest.raises(StepInputsFailed):
        run_step_with_audit('xm', cfg.paths.reports_dir, lambda: step6_xmatch.run(cfg))
    log = next((tmp_path / 'reports' / 'pipeline_audits').glob('*.log')).read_text()
    assert 'Status: FAILED' in log and 'missing.fits' in log


def test_image_wiring_preserves_structured_jobs_and_deduplicates_symlinks(tmp_path):
    raw = _fixtures(tmp_path)
    alias = tmp_path / 'alias.fits'
    alias.symlink_to(raw['input2']['path'])
    duplicate = deepcopy(raw)
    duplicate['input2']['path'] = str(alias)
    duplicate['radius'] = '0.0005555555555555556 deg'
    cfg = Config('wiring', paths=PathsCfg(sky_xmatches_dir=str(tmp_path / 'Sky-CrossMatches')),
                 reference=Target('CMC1', images=[str(tmp_path / 'cmc.fits')]),
                 tests=[Target('GPU', images=[str(tmp_path / 'gpu.fits')])])
    cfg.extra = {'survey_xmatch_jobs': [raw, duplicate],
                 'survey_xmatches': [{'name': 'racs', 'catalogue': raw['input2'],
                                      'targets': ['CMC1', 'GPU'], 'radius': '10 arcsec'}]}
    plan = wire_image_pipeline(cfg)
    assert len(plan.survey_xmatch_jobs) == 3
    assert len(plan.xmatch_pairs) == len(plan.positions) == 1
    assert len(plan.flux) == 0
    assert cfg.extra['survey_xmatch_jobs'][0]['radius'] == '2 arcsec'
    assert cfg.extra['survey_xmatch_jobs'][0]['input2']['ra_col'] == 'RAJ2000'
    assert build_image_pipeline_plan(cfg).survey_xmatch_jobs == plan.survey_xmatch_jobs
    assert resolve_survey_jobs(cfg) == list(plan.survey_xmatch_jobs)


def test_custom_initial_profile_frequency_and_coverage(tmp_path):
    raw = _fixtures(tmp_path)
    profiles = {'initial': {**CatalogueSpec.from_config(raw['input2']).__dict__,
                            'release': 'low3 INITIAL', 'frequency_hz': 943.5e6}}
    raw['input2'] = {'profile': 'initial', 'path': raw['input2']['path']}
    job = SurveyJob.from_config(raw, profiles)
    assert 'Preliminary' in job.input2.qualification
    assert job.input2.frequency_hz == 943.5e6
    raw['coverage'] = 'absent'
    with pytest.raises(ValueError, match='coverage_evidence'):
        SurveyJob.from_config(raw, profiles)


def test_frequency_uses_fits_axis_units_and_never_band_midpoints(tmp_path):
    raw = _fixtures(tmp_path)
    image = tmp_path / 'image.fits'
    h = fits.Header({'CTYPE3': 'STOKES', 'CRVAL3': 1, 'CTYPE4': 'FREQ',
                     'CRVAL4': 1234.567, 'CUNIT4': 'MHz', 'WCSAXES': 4,
                     'RESTFRQ': 949e6})
    fits.PrimaryHDU(header=h).writeto(image)
    raw['input1'].update(image=str(image), frequency_hz=949e6)
    spec = CatalogueSpec.from_config(raw['input1'])
    freq = frequency_provenance(spec, Table())
    assert freq['frequency_hz'] == pytest.approx(1234.567e6)
    assert 'CRVAL4' in freq['origin']
    raw['input1'].pop('image')
    raw['input1'].pop('frequency_hz')
    assert frequency_provenance(CatalogueSpec.from_config(raw['input1']), Table())['frequency_hz'] is None


def test_custom_profile_survives_yaml_path_expansion(tmp_path):
    yaml = tmp_path / 'config.yaml'
    yaml.write_text('''project_name: survey
path_vars: {root: /tmp/catalogues}
extra:
  catalogue_profiles:
    custom: {survey: RACS-mid, release: synthetic, grain: source, ra_col: x, dec_col: y, frame: icrs, coordinate_unit: deg, frequency_hz: 1367500000.0}
  survey_xmatch_jobs:
    - input1: {profile: pybdsf_source, path: '${root}/mk.fits'}
      input2: {profile: custom, path: '${root}/survey.xml', format: votable, table_id: catalogue}
      output: '${root}/match.fits'
      radius: 10 arcsec
''')
    jobs = resolve_survey_jobs(load_config(str(yaml)))
    assert jobs[0]['input2']['ra_col'] == 'x'
    assert jobs[0]['input2']['table_id'] == 'catalogue'
    assert jobs[0]['input2']['path'] == str(Path('/tmp/catalogues/survey.xml').resolve())


def test_pybdsf_gaul_grain_and_parent_source_counts(tmp_path):
    raw = _fixtures(tmp_path)
    t = Table({'RA': [10., 10.00005], 'DEC': [-20., -20.],
               'Source_id': [7, 7], 'Gaus_id': [8, 9]})
    t.write(raw['input1']['path'], overwrite=True)
    with pytest.raises(ValueError, match='Repeated source IDs'):
        CatalogueSpec.from_config(raw['input1']).prepare('T1')
    raw['input1']['profile'] = 'pybdsf_gaul'
    report = execute_survey_job(raw)
    assert report['associations'] == 4
    counts = report['catalogues'][0]
    assert counts['grain'] == 'component'
    assert counts['input'] == counts['matched'] == counts['ambiguous'] == 2
    assert counts['distinct_source_ids'] == counts['matched_source_ids'] == counts['ambiguous_source_ids'] == 1
    assert counts['unmatched_source_ids'] == 0
    out = Table.read(tmp_path / 'Sky-CrossMatches' / 'matches.fits')
    assert set(out['Gaus_id_1']) == {8, 9}


def test_no_candidates_with_valid_coordinates(tmp_path):
    raw = _fixtures(tmp_path)
    t = Table.read(raw['input1']['path'])
    t['RA'] = [200.]
    t.write(raw['input1']['path'], overwrite=True)
    report = execute_survey_job(raw)
    assert report['associations'] == 0
    assert report['catalogues'][0]['unmatched'] == report['catalogues'][0]['unmatched_source_ids'] == 1
    assert report['coverage']['status'] == 'unknown'
    assert report['coverage']['unmatched_input1'] == 1


def test_votable_multi_table_selection(tmp_path):
    from astropy.io.votable import from_table, writeto
    raw = _fixtures(tmp_path, xml=True)
    first = from_table(Table({'RAJ2000': [1.], 'DEJ2000': [1.], 'ID': ['x'], 'GID': ['a']}))
    second = from_table(Table({'RAJ2000': [2., 3.], 'DEJ2000': [2., 3.], 'ID': ['y', 'z'], 'GID': ['b', 'c']}))
    first.resources[0].tables[0].ID = 'first'
    second.resources[0].tables[0].ID = 'second'
    first.resources.append(second.resources[0])
    writeto(first, raw['input2']['path'])
    with pytest.raises(ValueError, match='Select table_id'):
        CatalogueSpec.from_config(raw['input2']).prepare('T')
    for selector in (1, 'second'):
        raw['input2']['table_id'] = selector
        assert len(CatalogueSpec.from_config(raw['input2']).prepare('T').table) == 2
    raw['input2']['table_id'] = 'missing'
    with pytest.raises(ValueError, match='not found'):
        CatalogueSpec.from_config(raw['input2']).prepare('T')


def test_low_high_inputs_explicit_and_template_opt_in(tmp_path):
    from meerkat_corr_imaging.config import LowHighSliceCfg
    raw = _fixtures(tmp_path)
    cfg = Config('bands', reference=Target('CMC1', images=[str(tmp_path / 'mfs.fits')],
                                         cuboid=str(tmp_path / 'cube.fits')),
                 low_high_slice=LowHighSliceCfg(enabled=True, add_to_source_finding=True))
    template = {'name': 'racs', 'catalogue': raw['input2'], 'radius': '10 arcsec'}
    cfg.extra = {'survey_xmatch_jobs': [raw], 'survey_xmatches': [template]}
    assert len(resolve_survey_jobs(cfg)) == 2  # explicit plus generated MFS
    template['bands'] = ['low', 'high']
    jobs = resolve_survey_jobs(cfg)
    assert len(jobs) == 3
    assert all('band' in job['input1']['path'] for job in jobs[1:])


def test_initial_qualification_survives_repeated_wiring(tmp_path):
    raw = _fixtures(tmp_path)
    raw['input2']['release'] = 'synthetic INITIAL'
    cfg = Config('initial', extra={'survey_xmatch_jobs': [raw]})
    first = wire_image_pipeline(cfg)
    second = wire_image_pipeline(cfg)
    assert first.survey_xmatch_jobs == second.survey_xmatch_jobs


def test_conflicting_outputs_rejected(tmp_path):
    raw = _fixtures(tmp_path)
    other = deepcopy(raw)
    other['radius'] = '20 arcsec'
    cfg = Config('conflict', extra={'survey_xmatch_jobs': [raw, other]})
    with pytest.raises(ValueError, match='same output'):
        resolve_survey_jobs(cfg)


@pytest.mark.parametrize('explicit_output', [False, True])
def test_legacy_wrapper_runs_real_matcher_with_angular_units(tmp_path, explicit_output):
    a, b = tmp_path / 'a.fits', tmp_path / 'b.fits'
    for path, offset in ((a, 0.), (b, .0001)):
        t = Table({'RA': [10. + offset], 'DEC': [-20.], 'Source_id': [1],
                   'Total_flux': [2.], 'Peak_flux': [1.]})
        t['Total_flux'].unit = u.Jy
        t.write(path)
    out = tmp_path / 'Sky-CrossMatches' / ('explicit.fits' if explicit_output else 'a_X_b.fits')
    pair = [str(a), str(b)]
    if explicit_output:
        pair.append(str(out))
    cfg = Config('legacy', paths=PathsCfg(sky_xmatches_dir=str(out.parent)),
                 extra={'xmatch_pairs': [pair]})
    run_step_with_audit('xm', str(tmp_path / 'reports'), lambda: step6_xmatch.run(cfg))
    t = Table.read(out)
    assert len(t) == 1
    assert 'total_flux_ratio' in t.colnames
    assert t['Total_flux_1'].unit == u.Jy
    assert t['sep_arcsec'][0] == pytest.approx(.338289, abs=.000001)


def test_legacy_collision_selection_not_reassigned(tmp_path):
    a, b = tmp_path / 'a.fits', tmp_path / 'b.fits'
    Table({'RA': [10., 10.0001], 'DEC': [0., 0.]}).write(a)
    Table({'RA': [10., 10.0005], 'DEC': [0., 0.]}).write(b)
    cats = [prepare_catalog(str(p), 'RA', 'DEC', 'icrs', 'T') for p in (a, b)]
    t, matches = cross_match_catalogues(*cats, parse_angle('2 arcsec'), False, enrich=False)
    assert [(m[0], m[1]) for m in matches] == [(0, 0)]
    assert t['ambiguous'][0]
    assert t.meta['HIERARCH XM1_AMBIGUOUS'] == 2
    assert t.meta['HIERARCH XM1_UNMATCHED'] == 1


def test_non_degree_angular_columns(tmp_path):
    p = tmp_path / 'angular.fits'
    t = Table({'RA': [np.pi / 12], 'DEC': [-1200.]})
    t['RA'].unit = u.rad
    t['DEC'].unit = u.arcmin
    t.write(p)
    cat = prepare_catalog(str(p), 'RA', 'DEC', 'icrs', 'T')
    assert cat.coords.ra.deg[0] == pytest.approx(15.)
    assert cat.coords.dec.deg[0] == pytest.approx(-20.)
    assert cat.table['RA'].unit == u.rad
