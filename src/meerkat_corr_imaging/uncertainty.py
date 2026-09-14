"""Reusable numerical uncertainties; no report-format calculations.

Formal errors are one standard deviation. Bootstrap CIs describe sampling of the
selected population, not calibration systematics. Missing errors stay NaN in
arrays and null in JSON. Independent catalogues are assumed unless covariance is
explicitly supplied to the numerical helpers.
"""
from dataclasses import dataclass
import json
from pathlib import Path
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from scipy.stats import norm

ONE_SIGMA = 0.6826894921370859


@dataclass(frozen=True)
class UncertaintyConfig:
    bootstrap_samples: int = 5000
    confidence_level: float = ONE_SIGMA
    random_seed: int = 20260911

    def __post_init__(self):
        if isinstance(self.bootstrap_samples, bool) or not isinstance(self.bootstrap_samples, (int, np.integer)) or self.bootstrap_samples < 2:
            raise ValueError('bootstrap_samples must be an integer >= 2')
        if not 0 < self.confidence_level < 1:
            raise ValueError('confidence_level must lie between 0 and 1')
        if isinstance(self.random_seed, bool) or not isinstance(self.random_seed, (int, np.integer)) or self.random_seed < 0:
            raise ValueError('random_seed must be a non-negative integer')

    @classmethod
    def from_mapping(cls, mapping=None):
        return cls(**(mapping or {}))


def column(table, name, unit=None, error=False):
    if name not in table.colnames:
        return np.full(len(table), np.nan)
    col = table[name]
    values = np.asarray(np.ma.filled(np.ma.asarray(col, dtype=float), np.nan))
    if unit is not None and getattr(col, 'unit', None) is not None:
        values = (values * col.unit).to_value(unit)
    values = np.where(np.isfinite(values), values, np.nan)
    return np.where(values > 0, values, np.nan) if error else values


def ratio_error(a, b, sigma_a, sigma_b, covariance=0.0):
    """A/B Jacobian, valid even for A=0; B=0 and invalid errors give NaN.

    Covariance is Cov(A,B) in product units, if known. Zero/negative formal
    standard errors are treated as unavailable, not exact measurements.
    """
    a, b, sa, sb, cov = np.broadcast_arrays(*map(lambda x: np.asarray(x, float), (a, b, sigma_a, sigma_b, covariance)))
    with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
        value = a / b
        variance = (sa / b)**2 + (a * sb / b**2)**2 - 2 * a * cov / b**3
    valid = np.isfinite(a) & np.isfinite(b) & (b != 0)
    good = valid & np.isfinite(sa) & np.isfinite(sb) & (sa > 0) & (sb > 0) & np.isfinite(cov) & (np.abs(cov) <= sa * sb) & (variance >= 0)
    return np.where(valid, value, np.nan), np.where(good, np.sqrt(np.maximum(variance, 0)), np.nan)


def interval(draws, value, config=None, method='paired-source percentile bootstrap', n=None):
    cfg = config or UncertaintyConfig()
    draws = np.asarray(draws, float)
    draws = draws[np.isfinite(draws)]
    result = dict(value=float(value), ci_low=None, ci_high=None, standard_error=None,
                  confidence_level=cfg.confidence_level, method=method, n=n)
    if draws.size >= 2:
        tail = (1 - cfg.confidence_level) / 2
        result.update(ci_low=float(np.quantile(draws, tail)), ci_high=float(np.quantile(draws, 1-tail)),
                      standard_error=float(np.std(draws, ddof=1)))
    else:
        result['reason'] = 'fewer than two independent samples or valid resamples'
    return result


def bootstrap(values, statistic=np.median, config=None):
    cfg = config or UncertaintyConfig()
    values = np.asarray(values, float)
    values = values[np.isfinite(values)]
    value = float(statistic(values)) if values.size else np.nan
    if values.size < 2:
        return interval([], value, cfg, n=int(values.size))
    rng = np.random.default_rng(cfg.random_seed)
    draws = []
    # Bound memory for large catalogues; support arbitrary scalar statistics.
    for start in range(0, cfg.bootstrap_samples, 128):
        indices = rng.integers(0, values.size, (min(128, cfg.bootstrap_samples-start), values.size))
        draws.extend(float(statistic(values[index])) for index in indices)
    return interval(draws, value, cfg, n=int(values.size))


def describe(values, config=None):
    functions = {'mean': np.mean, 'std': np.std, 'median': np.median,
                 'p16': lambda a: np.percentile(a, 16), 'p84': lambda a: np.percentile(a, 84),
                 'p95': lambda a: np.percentile(a, 95),
                 'nmad': lambda a: 1.4826*np.median(np.abs(a-np.median(a)))}
    return {key: bootstrap(values, function, config) for key, function in functions.items()}


def wilson(k, n, config=None):
    cfg = config or UncertaintyConfig()
    if n < 0 or k < 0 or k > n or int(k) != k or int(n) != n:
        raise ValueError('Wilson requires integer 0 <= k <= n')
    result = interval([], k/n if n else np.nan, cfg, method='Wilson binomial sampling interval', n=int(n))
    if n:
        z = norm.ppf((1 + cfg.confidence_level)/2)
        p = k/n
        center = (p+z*z/(2*n))/(1+z*z/n)
        radius = z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/(1+z*z/n)
        result.update(ci_low=max(0., center-radius), ci_high=min(1., center+radius))
        result.pop('reason', None)
    return result


def measurement(value, error, config=None, method='propagated independent measurement errors'):
    cfg = config or UncertaintyConfig()
    result = interval([], value, cfg, method=method)
    if np.isfinite(value) and error is not None and np.isfinite(error) and error >= 0:
        z = norm.ppf((1+cfg.confidence_level)/2)
        result.update(ci_low=float(value-z*error), ci_high=float(value+z*error), standard_error=float(error))
        result.pop('reason', None)
    else:
        result['reason'] = 'measurement uncertainty unavailable'
    return result


def decision(value, uncertainty, limit, *, target=None):
    """Strict interval containment; preserve the original point decision too."""
    finite = value is not None and np.isfinite(value)
    deviation = abs(value-target) if finite and target is not None else value
    point = 'Not assessed' if not finite else ('Pass' if deviation < limit else 'Concern')
    lo, hi = (uncertainty or {}).get('ci_low'), (uncertainty or {}).get('ci_high')
    if not finite or lo is None or hi is None or not np.isfinite([lo, hi]).all():
        return dict(status='Not assessed', point_status=point, reason='decision uncertainty unavailable')
    bound = max(abs(lo-target), abs(hi-target), deviation) if target is not None else max(value, hi)
    return dict(status='Pass' if bound < limit else 'Concern', point_status=point,
                reason='entire interval satisfies limit' if bound < limit else 'estimate or interval reaches/exceeds limit',
                decision_bound=float(bound), ci_low=float(lo), ci_high=float(hi))


def _offsets(coords):
    ra1, dec1, ra2, dec2 = coords.T
    a = SkyCoord(ra1*u.deg, dec1*u.deg)
    b = SkyCoord(ra2*u.deg, dec2*u.deg)
    east, north = a.spherical_offsets_to(b)
    return np.column_stack((east.arcsec, north.arcsec))


def position_covariance(ra1, dec1, ra2, dec2, errors, covariance=None):
    """Exact spherical-offset numerical Jacobian, including RA cos(dec).

    errors: N x 4 coordinate SDs in degrees. Optional N x 4 x 4 input
    covariance replaces diagonal independence; output is N x 2 x 2 arcsec².
    """
    coords = np.column_stack((ra1, dec1, ra2, dec2)).astype(float)
    errors = np.asarray(errors, float)
    good = np.isfinite(coords).all(axis=1) & np.isfinite(errors).all(axis=1) & (errors >= 0).all(axis=1)
    out = np.full((len(coords), 2, 2), np.nan)
    if not good.any():
        return out
    x = coords[good]
    jac = np.empty((len(x), 2, 4))
    for j in range(4):
        step = 1.e-5 if j in (0, 2) else 1.e-6
        plus, minus = x.copy(), x.copy()
        plus[:, j] += step
        minus[:, j] -= step
        jac[:, :, j] = (_offsets(plus)-_offsets(minus))/(2*step)
    cov = np.zeros((len(x), 4, 4))
    cov[:, np.arange(4), np.arange(4)] = errors[good]**2
    if covariance is not None:
        cov = np.asarray(covariance, float)[good]
        if not np.allclose(cov, cov.swapaxes(1, 2)) or np.any(np.linalg.eigvalsh(cov) < -1.e-20):
            raise ValueError('coordinate covariance must be symmetric positive semidefinite')
    out[good] = jac @ cov @ jac.swapaxes(1, 2)
    return out


def radial_errors(xy, cov, config=None):
    cfg = config or UncertaintyConfig()
    xy, cov = np.asarray(xy, float), np.asarray(cov, float)
    radius = np.linalg.norm(xy, axis=1)
    sd = np.full(len(xy), np.nan)
    low, high = sd.copy(), sd.copy()
    theta_sd = sd.copy()
    rng = np.random.default_rng(cfg.random_seed)
    for i in range(len(xy)):
        if not np.isfinite(cov[i]).all() or not np.isfinite(xy[i]).all():
            continue
        scale = np.sqrt(max(np.linalg.eigvalsh(cov[i])))
        if radius[i] > 3*scale:
            grad = xy[i]/radius[i]
            sd[i] = np.sqrt(max(0, grad @ cov[i] @ grad))
            m = measurement(radius[i], sd[i], cfg)
            low[i], high[i] = m['ci_low'], m['ci_high']
        else:
            draws = rng.multivariate_normal(xy[i], cov[i], size=cfg.bootstrap_samples)
            m = interval(np.linalg.norm(draws, axis=1), radius[i], cfg, method='Gaussian coordinate Monte Carlo')
            sd[i], low[i], high[i] = m['standard_error'], m['ci_low'], m['ci_high']
        if radius[i] > 3*scale:
            g = np.array([-xy[i, 1], xy[i, 0]])/radius[i]**2
            theta_sd[i] = np.rad2deg(np.sqrt(max(0, g @ cov[i] @ g)))
    return sd, low, high, theta_sd


def enrich_matches(table, config=None):
    """Preserve original columns and append normalized measurements/errors."""
    t = table.copy()
    for suffix in ('1', '2'):
        for source, name, unit in [('RA', 'ra', u.deg), ('DEC', 'dec', u.deg),
                                   ('Total_flux', 'total_flux', u.Jy), ('Peak_flux', 'peak_flux', u.Jy/u.beam)]:
            for prefix, tail in [('', ''), ('E_', '_err')]:
                key = f'{name}{tail}_{suffix}'
                t[key] = column(t, f'{prefix}{source}_{suffix}', unit, error=bool(prefix))
                t[key].unit = unit
    for label in ('total', 'peak'):
        a, b = (column(t, f'{label}_flux_{s}') for s in ('2', '1'))
        sa, sb = (column(t, f'{label}_flux_err_{s}') for s in ('2', '1'))
        value, error = ratio_error(a, b, sa, sb)
        t[f'{label}_flux_ratio'] = value
        t[f'{label}_flux_ratio_err'] = error
        t[f'{label}_flux_fractional_difference'] = value-1
        t[f'{label}_flux_fractional_difference_err'] = error
    coords = [column(t, name) for name in ('ra_1', 'dec_1', 'ra_2', 'dec_2')]
    good = np.isfinite(coords).all(axis=0)
    good &= (np.abs(coords[1]) <= 90) & (np.abs(coords[3]) <= 90)
    xy = np.full((len(t), 2), np.nan)
    if good.any():
        xy[good] = _offsets(np.column_stack(coords)[good])
    errors = np.column_stack([column(t, name) for name in ('ra_err_1', 'dec_err_1', 'ra_err_2', 'dec_err_2')])
    cov = np.full((len(t), 2, 2), np.nan)
    cov[good] = position_covariance(*(c[good] for c in coords), errors[good])
    radius = np.full(len(t), np.nan)
    if good.any():
        a = SkyCoord(coords[0][good]*u.deg, coords[1][good]*u.deg)
        b = SkyCoord(coords[2][good]*u.deg, coords[3][good]*u.deg)
        radius[good] = a.separation(b).arcsec
    sd, lo, hi, theta_sd = radial_errors(xy, cov, config)
    for key, data, unit in [('east_offset_arcsec', xy[:, 0], u.arcsec), ('north_offset_arcsec', xy[:, 1], u.arcsec),
                            ('east_offset_err_arcsec', np.sqrt(cov[:, 0, 0]), u.arcsec),
                            ('north_offset_err_arcsec', np.sqrt(cov[:, 1, 1]), u.arcsec),
                            ('east_north_cov_arcsec2', cov[:, 0, 1], u.arcsec**2),
                            ('separation_arcsec', radius, u.arcsec), ('separation_err_arcsec', sd, u.arcsec),
                            ('separation_ci_low_arcsec', lo, u.arcsec), ('separation_ci_high_arcsec', hi, u.arcsec),
                            ('offset_angle_err_deg', theta_sd, u.deg)]:
        t[key] = data
        t[key].unit = unit
    return t


def format_uncertainty(value, info=None, *, scale=1., digits=3):
    if value is None or not np.isfinite(value):
        return '—'
    value *= scale
    info = info or {}
    lo, hi = info.get('ci_low'), info.get('ci_high')
    if lo is None or hi is None or not np.isfinite([lo, hi]).all():
        return f'{value:.{digits}g} (unc. unavailable)'
    lo, hi = lo*scale, hi*scale
    lowerr, higherr = value-lo, hi-value
    width = max(abs(lowerr), abs(higherr))
    decimals = max(0, min(10, 1-int(np.floor(np.log10(width))))) if width > 0 else digits
    if lowerr >= 0 and higherr >= 0 and width > 0 and abs(lowerr-higherr) <= .1*width:
        return f'{value:.{decimals}f} ± {(lowerr+higherr)/2:.{decimals}f}'
    return f'{value:.{decimals}f} [{lo:.{decimals}f}, {hi:.{decimals}f}]'


def write_json(path, payload):
    def clean(obj):
        if isinstance(obj, dict):
            return {str(k): clean(v) for k, v in obj.items()}
        if isinstance(obj, (list, tuple, np.ndarray)):
            return [clean(v) for v in obj]
        if isinstance(obj, (float, np.floating)):
            return float(obj) if np.isfinite(obj) else None
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.bool_):
            return bool(obj)
        return obj
    Path(path).write_text(json.dumps(clean(payload), indent=2, sort_keys=True, allow_nan=False)+'\n')


def linear_fit(x, y, xerr=None, yerr=None, config=None, through_origin=False):
    """Orthogonal-distance fit if both formal axes exist; otherwise paired OLS.

    Complete-error rows are used for weighted fits; missing errors never receive
    fabricated weights. Bootstrap resamples the adopted fitted population.
    """
    from scipy import odr
    cfg = config or UncertaintyConfig()
    x, y = np.asarray(x, float), np.asarray(y, float)
    finite = np.isfinite(x) & np.isfinite(y)
    xe = np.full(x.shape, np.nan) if xerr is None else np.asarray(xerr, float)
    ye = np.full(y.shape, np.nan) if yerr is None else np.asarray(yerr, float)
    complete = finite & np.isfinite(xe) & np.isfinite(ye) & (xe > 0) & (ye > 0)
    p = 1 if through_origin else 2
    weighted = complete.sum() > p
    mask = complete if weighted else finite
    x, y, xe, ye = x[mask], y[mask], xe[mask], ye[mask]
    if len(x) <= p or (not through_origin and np.ptp(x) == 0) or np.sum(x*x) == 0:
        return {'slope': np.nan, 'intercept': 0. if through_origin else np.nan,
                'slope_err': np.nan, 'intercept_err': np.nan, 'covariance': None,
                'uncertainties': {}, 'method': 'unidentifiable fit', 'n_fit': len(x)}
    def solve(ids):
        xx, yy = x[ids], y[ids]
        design = xx[:, None] if through_origin else np.column_stack((xx, np.ones(len(xx))))
        if np.linalg.matrix_rank(design) < p:
            raise ValueError('degenerate resample')
        beta = np.linalg.lstsq(design, yy, rcond=None)[0]
        if weighted:
            model = odr.Model((lambda b, xx: b[0]*xx) if through_origin else (lambda b, xx: b[0]*xx+b[1]))
            fit = odr.ODR(odr.RealData(xx, yy, sx=xe[ids], sy=ye[ids]), model, beta0=beta).run()
            if fit.info > 4:
                raise ValueError('ODR did not converge')
            return fit.beta, fit.cov_beta
        residual = yy-design @ beta
        return beta, np.linalg.inv(design.T @ design)*np.sum(residual**2)/(len(xx)-p)
    beta, cov = solve(np.arange(len(x)))
    rng = np.random.default_rng(cfg.random_seed)
    draws = []
    for _ in range(cfg.bootstrap_samples):
        try:
            b, _ = solve(rng.integers(0, len(x), len(x)))
            draws.append(b)
        except (ValueError, np.linalg.LinAlgError):
            continue
    draws = np.asarray(draws).reshape(-1, p)
    summaries = {}
    for i, key in enumerate(('slope', 'intercept')[:p]):
        info = interval(draws[:, i], beta[i], cfg, n=len(x))
        formal = measurement(beta[i], np.sqrt(max(0, cov[i, i])), cfg,
            method='absolute ODR measurement covariance' if weighted else 'residual-scaled OLS covariance')
        if weighted:
            info['sampling'] = info.copy()
            info['formal'] = formal
            info['ci_low'] = min(formal['ci_low'], info['ci_low']) if info['ci_low'] is not None else formal['ci_low']
            info['ci_high'] = max(formal['ci_high'], info['ci_high']) if info['ci_high'] is not None else formal['ci_high']
            info['standard_error'] = max(formal['standard_error'], info['standard_error'] or 0)
            info['method'] = 'envelope of formal ODR and paired-source bootstrap'
        summaries[key] = info
    return {'slope': float(beta[0]), 'intercept': 0. if through_origin else float(beta[1]),
            'slope_err': summaries['slope']['standard_error'] if summaries['slope']['standard_error'] is not None else np.nan,
            'intercept_err': 0. if through_origin else (summaries['intercept']['standard_error'] if summaries['intercept']['standard_error'] is not None else np.nan),
            'covariance': cov.tolist(), 'covariance_parameter_order': ['slope'] if through_origin else ['slope','intercept'],
            'uncertainties': summaries, 'method': 'both-axis ODR' if weighted else 'unweighted OLS; sampling uncertainty', 'n_fit': len(x)}


def paired_summary(arrays, function, keys, config=None, periods=None, minimum=2):
    """Resample complete source rows; unwrap periodic parameters at the estimate."""
    cfg = config or UncertaintyConfig()
    arrays = [np.asarray(a, float) for a in arrays]
    mask = np.isfinite(arrays).all(axis=0)
    arrays = [a[mask] for a in arrays]
    n = len(arrays[0])
    central = function(*arrays)
    draws = {key: [] for key in keys}
    rng = np.random.default_rng(cfg.random_seed)
    if n >= minimum:
        for _ in range(cfg.bootstrap_samples):
            ids = rng.integers(0, n, n)
            estimate = function(*(a[ids] for a in arrays))
            for key in keys:
                value = estimate[key]
                period = (periods or {}).get(key)
                if period and np.isfinite(central[key]):
                    value = central[key]+(value-central[key]+period/2)%period-period/2
                draws[key].append(value)
    return {key: interval(draws[key], central[key], cfg, n=n) for key in keys}


def cluster_comparison(reference, test, config=None):
    """Independent scan-cluster bootstrap of two dataset medians and B-A.

    Inputs are lists of scan arrays; whole scans retain baseline dependence.
    Independence between scans/datasets is an explicit sampling assumption.
    A single scan cannot establish scan-to-scan sampling uncertainty.
    """
    cfg = config or UncertaintyConfig()
    rng = np.random.default_rng(cfg.random_seed)
    values, distributions, summaries = [], [], []
    for groups in (reference, test):
        groups = [np.asarray(g, float)[np.isfinite(g)] for g in groups]
        groups = [g for g in groups if len(g)]
        value = float(np.median(np.concatenate(groups))) if groups else np.nan
        draws = []
        if len(groups) >= 2:
            for _ in range(cfg.bootstrap_samples):
                ids = rng.integers(0, len(groups), len(groups))
                draws.append(np.median(np.concatenate([groups[i] for i in ids])))
        values.append(value)
        distributions.append(np.asarray(draws))
        summaries.append(interval(draws, value, cfg, method='whole-scan cluster bootstrap', n=len(groups)))
    difference = distributions[1]-distributions[0] if all(len(d) for d in distributions) else []
    delta = interval(difference, values[1]-values[0], cfg, method='independent whole-scan cluster bootstrap difference B-A')
    return {'reference': summaries[0], 'test': summaries[1], 'difference': delta,
            'measurement_uncertainty': 'unavailable: exported spectral summaries contain no measurement covariance'}


def envelope(sampling, propagated):
    """Keep both components, use their envelope without double-counting noise.

    This is a conservative sensitivity envelope of two nominal-confidence
    intervals, not an assertion of exact frequentist coverage for a combined CI.
    """
    result = sampling.copy()
    result['sampling'] = sampling.copy()
    result['measurement'] = propagated
    available = [info for info in (sampling, propagated)
                 if info.get('ci_low') is not None and info.get('ci_high') is not None]
    if available:
        result['ci_low'] = min(info['ci_low'] for info in available)
        result['ci_high'] = max(info['ci_high'] for info in available)
        errors = [info['standard_error'] for info in available if info.get('standard_error') is not None]
        result['standard_error'] = max(errors) if errors else None
        result.pop('reason', None)
    result['method'] = ('envelope of sampling and propagated measurement intervals'
                        if len(available) == 2 else available[0]['method'] if available else 'uncertainty unavailable')
    return result


def catalogue_intervals(table, config=None):
    """Catalogue aggregates with separate source-sampling and measurement CIs.

    Formal-error Monte Carlo perturbs the fixed matched sample; bootstrap
    resamples observed sources. Their envelope retains a measurement floor
    without adding measurement noise to an already noisy source bootstrap.
    Missing formal errors make the whole-population measurement CI unavailable;
    they never remove valid sources from the sampling statistic.
    """
    cfg = config or UncertaintyConfig()
    rng = np.random.default_rng(cfg.random_seed)
    output = {}
    for kind, prefix in [('total', 'total_flux'), ('peak', 'peak_flux')]:
        a, b = column(table, kind+'_flux_2'), column(table, kind+'_flux_1')
        ea, eb = column(table, kind+'_flux_err_2'), column(table, kind+'_flux_err_1')
        valid = np.isfinite(a) & np.isfinite(b) & (a > 0) & (b > 0)
        a, b, ea, eb = (v[valid] for v in (a,b,ea,eb))
        ratio = a/b
        specs = [(prefix+'_ratio_median', lambda x, axis=None: np.median(x, axis=axis))]
        if kind == 'total':
            specs.append(('total_flux_abs_deviation_p95', lambda x, axis=None: np.percentile(np.abs(x-1), 95, axis=axis)))
        complete = len(a) > 0 and np.isfinite(ea).all() and np.isfinite(eb).all() and (ea > 0).all() and (eb > 0).all()
        distributions = {key: [] for key, _ in specs}
        if complete:
            for start in range(0, cfg.bootstrap_samples, 128):
                size = (min(128, cfg.bootstrap_samples-start), len(a))
                aa, bb = rng.normal(a, ea, size), rng.normal(b, eb, size)
                with np.errstate(divide='ignore', invalid='ignore'):
                    simulated = aa/bb
                for key, func in specs:
                    distributions[key].extend(func(simulated, axis=1))
        for key, func in specs:
            sampled = bootstrap(ratio, func, cfg)
            propagated = interval(distributions[key], sampled['value'], cfg, method='independent Gaussian flux-error Monte Carlo of fixed matched sample', n=len(a))
            if not complete:
                propagated['reason'] = 'whole-sample measurement propagation unavailable: missing/invalid formal flux errors'
            output[key] = envelope(sampled, propagated)
        if kind == 'total':
            output['total_flux_fraction_beyond_limit'] = wilson(int(np.sum(np.abs(ratio-1) > .05)), len(ratio), cfg)
    xy = np.column_stack((column(table, 'east_offset_arcsec'), column(table, 'north_offset_arcsec')))
    radius = column(table, 'separation_arcsec')
    valid = np.isfinite(radius)
    xy, radius = xy[valid], radius[valid]
    cov = np.zeros((len(table),2,2))
    cov[:,0,0] = column(table, 'east_offset_err_arcsec')**2
    cov[:,1,1] = column(table, 'north_offset_err_arcsec')**2
    cov[:,0,1] = cov[:,1,0] = column(table, 'east_north_cov_arcsec2')
    cov = cov[valid]
    specs = [('separation_median_arcsec', lambda x, axis=None: np.median(x, axis=axis)),
             ('separation_p95_arcsec', lambda x, axis=None: np.percentile(x, 95, axis=axis)),
             ('separation_max_arcsec', lambda x, axis=None: np.max(x, axis=axis))]
    distributions = {key: [] for key, _ in specs}
    complete = len(radius) > 0 and np.isfinite(cov).all()
    if complete:
        eigenvalues, eigenvectors = np.linalg.eigh(cov)
        root = eigenvectors * np.sqrt(np.maximum(eigenvalues, 0))[:, None, :]
        for start in range(0, cfg.bootstrap_samples, 128):
            size = (min(128, cfg.bootstrap_samples-start), len(radius), 2)
            draws = xy + np.einsum('nij,bnj->bni', root, rng.normal(size=size))
            radii = np.linalg.norm(draws, axis=2)
            for key, func in specs:
                distributions[key].extend(func(radii, axis=1))
    for key, func in specs:
        sampled = (interval([], np.max(radius) if len(radius) else np.nan, cfg,
                   method='observed maximum; population-tail sampling CI unavailable', n=len(radius))
                   if key.endswith('max_arcsec') else bootstrap(radius, func, cfg))
        propagated = interval(distributions[key], sampled['value'], cfg,
            method='Gaussian positional-covariance Monte Carlo of fixed matched sample', n=len(radius))
        if not complete:
            propagated['reason'] = 'whole-sample measurement propagation unavailable: missing/invalid formal position errors'
        output[key] = envelope(sampled, propagated)
    return output
