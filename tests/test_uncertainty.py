import json
import numpy as np
import pytest
from astropy import units as u
from astropy.io import fits
from astropy.table import Table, MaskedColumn
from scipy.stats import binomtest, norm

from meerkat_corr_imaging.uncertainty import (
    UncertaintyConfig, ONE_SIGMA, ratio_error, column, enrich_matches,
    position_covariance, radial_errors, bootstrap, wilson, measurement,
    decision, format_uncertainty, linear_fit, write_json)
from meerkat_corr_imaging.verification_report import (
    fit_rigid_transform, _robust_rms, _quality_rows, _score_band, _overall_status)

FAST = UncertaintyConfig(200, ONE_SIGMA, 3)


def test_ratio_jacobian_zero_and_covariance():
    ratio, error = ratio_error([4, 0, 1, 1, 1], [2, 2, 0, 2, 2],
                               [.4, .4, .1, 0, -1], [.1]*5)
    assert ratio[:2] == pytest.approx([2, 0])
    assert error[:2] == pytest.approx([np.sqrt(.05), .2])
    assert np.isnan(error[2:]).all()
    assert np.isnan(ratio[2])
    _, correlated = ratio_error(4, 2, .4, .1, covariance=.02)
    assert correlated == pytest.approx(np.sqrt(.03))
    assert np.isnan(ratio_error(4, 2, .4, .1, covariance=1)[1])


def test_positional_cos_dec_and_covariance():
    errors = np.array([[1, 2, 3, 4]])/3600
    cov = position_covariance([30], [-60], [30], [-60], errors)
    assert cov[0, 0, 0] == pytest.approx(2.5, rel=1e-6)
    assert cov[0, 1, 1] == pytest.approx(20, rel=1e-6)
    assert abs(cov[0, 0, 1]) < 1e-6
    full = np.diag(errors[0]**2)[None]
    full[0, 0, 2] = full[0, 2, 0] = .5/3600**2
    changed = position_covariance([30], [-60], [30], [-60], errors, full)
    assert changed[0, 0, 0] == pytest.approx(2.25, rel=1e-6)


def test_radial_linear_and_origin_monte_carlo():
    cov = np.array([np.diag([.04, .09]), np.eye(2)])
    sd, lo, hi, theta = radial_errors([[3, 4], [0, 0]], cov, UncertaintyConfig())
    assert sd[0] == pytest.approx(np.sqrt((3/5)**2*.04+(4/5)**2*.09))
    # At the origin radius is Rayleigh, not a zero-error scalar measurement.
    assert sd[1] == pytest.approx(np.sqrt((4-np.pi)/2), rel=.04)
    assert 0 < lo[1] < hi[1]
    assert np.isnan(theta[1])


def test_normalization_units_masks_and_missing_errors():
    t = Table({'RA_1': [1., 2.], 'DEC_1': [-60., -60.], 'RA_2': [1., 2.], 'DEC_2': [-60., -60.],
               'Total_flux_1': [1000., 0.], 'Total_flux_2': [2., 2.]})
    t['Total_flux_1'].unit = u.mJy
    t['Total_flux_2'].unit = u.Jy
    t['E_RA_1'] = MaskedColumn([1., 999.], mask=[False, True], unit=u.arcsec)
    t['E_Total_flux_1'] = [100., -1.]*u.mJy
    t['E_Total_flux_2'] = [.2, np.inf]*u.Jy
    result = enrich_matches(t, FAST)
    assert result['total_flux_ratio'][0] == 2
    assert result['total_flux_ratio_err'][0] == pytest.approx(np.sqrt(.08))
    assert result['ra_err_1'][0] == pytest.approx(1/3600)
    assert np.isnan(result['ra_err_1'][1])
    assert np.isnan(result['separation_err_arcsec']).all()
    assert len(result) == len(t)
    assert np.isnan(result['total_flux_ratio'][1])


def test_bootstrap_reproducible_and_small_samples():
    x = np.arange(1., 21.)
    for statistic in [np.median, lambda a: np.percentile(a, 95)]:
        a = bootstrap(x, statistic, FAST)
        b = bootstrap(x, statistic, FAST)
        assert a == b
        assert a['ci_low'] < a['ci_high']
        assert a['standard_error'] > 0
    assert bootstrap([5], config=FAST)['ci_low'] is None
    assert bootstrap([], config=FAST)['n'] == 0
    assert bootstrap([3]*10, config=FAST)['ci_low'] == 3


@pytest.mark.parametrize('k,n', [(0,10),(10,10),(3,10),(1,1)])
def test_wilson_matches_scipy(k, n):
    actual = wilson(k, n)
    expected = binomtest(k, n).proportion_ci(confidence_level=ONE_SIGMA, method='wilson')
    assert actual['ci_low'] == pytest.approx(expected.low)
    assert actual['ci_high'] == pytest.approx(expected.high)


def test_fraction_empty():
    assert wilson(0, 0)['ci_high'] is None
    with pytest.raises(ValueError):
        wilson(5, 4)


def test_rigid_gls_recovers_params_and_formal_floor():
    rng = np.random.default_rng(10)
    reference = rng.normal(size=(40, 2))*500
    angle = np.deg2rad(.012)
    rot = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
    test = reference @ rot.T+[.3, -.2]
    cov = np.tile(np.diag([.01, .04]), (40,1,1))
    result = fit_rigid_transform(reference, test, bootstrap_samples=100,
                                 reference_covariance=cov, test_covariance=cov)
    assert result.rotation_deg == pytest.approx(.012, abs=1e-8)
    assert result.translation_east_arcsec == pytest.approx(.3, abs=1e-8)
    assert result.translation_north_arcsec == pytest.approx(-.2, abs=1e-8)
    assert result.rotation_uncertainty_deg > 0
    assert result.uncertainties['translation_east_arcsec']['formal']['standard_error'] > 0
    assert np.linalg.eigvalsh(result.parameter_covariance).min() > 0
    assert len(result.residual_err_arcsec) == 40


def test_rigid_downweights_bad_source_and_falls_back():
    rng = np.random.default_rng(10)
    ref = rng.normal(size=(40,2))*100
    test = ref+[.2, -.1]+rng.normal(0, .01, ref.shape)
    test[0] += [1, 1]
    cov = np.tile(np.eye(2)*.01, (40,1,1))
    cov[0] *= 1000
    fit = fit_rigid_transform(ref, test, bootstrap_samples=50,
                              reference_covariance=cov, test_covariance=cov)
    assert fit.translation_east_arcsec == pytest.approx(.2, abs=.01)
    cov[:] = np.nan
    fallback = fit_rigid_transform(ref, test, bootstrap_samples=50,
                              reference_covariance=cov, test_covariance=cov)
    assert 'unweighted' in fallback.fit_method
    with pytest.raises(ValueError, match='unidentifiable'):
        fit_rigid_transform(np.zeros((10,2)), np.zeros((10,2)), bootstrap_samples=5)


def test_linear_fit_uses_both_axes_and_preserves_covariance():
    x = np.arange(1., 11.)
    fit = linear_fit(x, 2*x+.5, np.ones(10)*.1, np.ones(10)*.2, FAST)
    assert fit['slope'] == pytest.approx(2, abs=1e-7)
    assert fit['intercept'] == pytest.approx(.5, abs=1e-7)
    assert fit['slope_err'] > 0
    assert np.shape(fit['covariance']) == (2,2)
    fallback = linear_fit(x, 2*x+.5, np.full(10, np.nan), np.ones(10)*.2, FAST)
    assert 'unweighted' in fallback['method']


def test_acceptance_interval_boundary_missing_and_flux_sides():
    assert decision(.9, measurement(.9,.05), 1)['status'] == 'Pass'
    assert decision(.99, measurement(.99,.02), 1)['status'] == 'Concern'
    assert decision(.9, None, 1)['status'] == 'Not assessed'
    for value in (.96, 1.04):
        result = decision(value, measurement(value,.02), .05, target=1.)
        assert result['point_status'] == 'Pass'
        assert result['status'] == 'Concern'
    assert decision(1., measurement(1., .01), .05, target=1.)['status'] == 'Pass'


def test_formatting_and_strict_json(tmp_path):
    assert format_uncertainty(1, measurement(1,.1)) == '1.00 ± 0.10'
    assert '[' in format_uncertainty(1, {'ci_low': .99, 'ci_high':1.5})
    assert '[' in format_uncertainty(0, {'ci_low': .3, 'ci_high':1.5})
    assert 'unavailable' in format_uncertainty(1)
    assert format_uncertainty(np.nan) == '—'
    path = tmp_path/'result.json'
    write_json(path, {'value': np.nan, 'array': np.array([np.inf,1])})
    assert json.loads(path.read_text()) == {'value':None, 'array':[None,1]}


def _image(path, beam=True):
    rng = np.random.default_rng(8)
    header = fits.Header({'CTYPE1':'RA---TAN', 'CTYPE2':'DEC--TAN',
        'CRVAL1':30., 'CRVAL2':-60., 'CRPIX1':101., 'CRPIX2':101.,
        'CDELT1':-.006, 'CDELT2':.006})
    if beam:
        header['BMAJ'], header['BMIN'] = .02, .01
    fits.writeto(path, rng.normal(0, 1e-5, (201,201)), header)


def test_rms_uses_mad_constant_effective_beams(tmp_path):
    path = tmp_path/'image.fits'
    _image(path)
    value, info = _robust_rms(path, stride=1, return_details=True)
    assert value == _robust_rms(path, stride=1)
    assert 1 < info['effective_independent_beams'] < info['n_pixels']
    expected = value*1.4826/(4*norm.pdf(norm.ppf(.75))*np.sqrt(info['effective_independent_beams']))
    assert info['standard_error'] == pytest.approx(expected)
    missing = tmp_path/'missing.fits'
    _image(missing, beam=False)
    assert _robust_rms(missing, return_details=True)[1]['ci_low'] is None


def test_score_missing_peak_errors_does_not_crash(tmp_path):
    image = tmp_path/'image.fits'
    _image(image)
    t = Table({'RA_1':[30.], 'RA_2':[30.], 'DEC_1':[-60.], 'DEC_2':[-60.],
               'Peak_flux_1':[1.], 'Peak_flux_2':[1.]})
    cat = tmp_path/'match.fits'
    t.write(cat)
    assert len(_quality_rows(t)) == 0
    result = _score_band({'band':'MFS','reference_image':str(image), 'test_image':str(image), 'xmatch_table':str(cat)}, (1,1), FAST)
    assert result.quality_count == 0
    assert result.position_status == 'Not assessed'
    assert result.flux_status == 'Not assessed'
    assert result.rms_ratio == 1
    assert _overall_status([result], 'Reviewed') == 'Partial'


@pytest.mark.parametrize('kwargs', [{'bootstrap_samples':0}, {'bootstrap_samples':3.5},
                                   {'confidence_level':1}, {'random_seed':-1}])
def test_config_validation(kwargs):
    with pytest.raises(ValueError):
        UncertaintyConfig(**kwargs)


def test_cluster_comparison_preserves_whole_scans():
    from meerkat_corr_imaging.uncertainty import cluster_comparison
    result = cluster_comparison([np.ones(10), np.ones(10)*3], [np.ones(10)*2, np.ones(10)*4], FAST)
    assert result['difference']['value'] == 1
    assert result['reference']['n'] == 2
    assert result['reference']['ci_high']-result['reference']['ci_low'] == 2
    single = cluster_comparison([np.ones(1000)], [np.ones(1000)*2], FAST)
    assert single['difference']['ci_low'] is None


def test_catalogue_decision_retains_measurement_floor():
    from meerkat_corr_imaging.uncertainty import catalogue_intervals
    n=20
    t=Table({'RA_1': np.arange(n)*.01+30, 'DEC_1':np.full(n,-60.),
             'RA_2': np.arange(n)*.01+30, 'DEC_2':np.full(n,-60.),
             'Total_flux_1':np.ones(n), 'Total_flux_2':np.full(n,1.04),
             'E_Total_flux_1':np.full(n,.15), 'E_Total_flux_2':np.full(n,.15)})
    enriched=enrich_matches(t, FAST)
    info=catalogue_intervals(enriched, FAST)['total_flux_ratio_median']
    assert info['sampling']['ci_low'] == info['sampling']['ci_high'] == 1.04
    assert info['measurement']['ci_high'] > 1.05
    assert decision(1.04,info,.05,target=1.)['status'] == 'Concern'
    t.remove_column('E_Total_flux_1')
    missing=catalogue_intervals(enrich_matches(t,FAST),FAST)['total_flux_ratio_median']
    assert missing['measurement']['ci_high'] is None
    assert missing['sampling']['ci_high'] == 1.04


def test_postfit_interval_retains_measurement_uncertainty():
    rng=np.random.default_rng(32)
    xy=rng.normal(size=(20,2))*100
    cov=np.tile(np.eye(2)*.01,(20,1,1))
    fit=fit_rigid_transform(xy,xy+[.2,-.1],bootstrap_samples=100,
                           reference_covariance=cov,test_covariance=cov)
    info=fit.uncertainties['postfit_p95_arcsec']
    assert info['sampling']['ci_high'] < 1e-8
    assert info['measurement']['ci_high'] > .1
    assert info['ci_high'] > .1
