"""Rigid astrometry with full E/N covariance and paired-source resampling."""
import numpy as np
from .uncertainty import UncertaintyConfig, interval, measurement, radial_errors, envelope


def _rotation(angle):
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[c, -s], [s, c]])


def solve_rigid(reference, test, ref_cov=None, test_cov=None):
    """Iterated GLS for exact rotation; covariance order east, north, radians.

    Positional covariance is rotated at each iteration. Formal covariance is
    absolute (not reduced-chi-square scaled); known catalogue error floors must
    not disappear when a small sample happens to fit perfectly.
    """
    ref_center, test_center = reference.mean(axis=0), test.mean(axis=0)
    left, _, right = np.linalg.svd((reference-ref_center).T @ (test-test_center))
    rotation = right.T @ left.T
    if np.linalg.det(rotation) < 0:
        right[-1] *= -1
        rotation = right.T @ left.T
    angle = np.arctan2(rotation[1, 0], rotation[0, 0])
    translation = test_center-rotation @ ref_center
    weighted = ref_cov is not None and test_cov is not None
    for _ in range(12):
        rotation = _rotation(angle)
        rotated = reference @ rotation.T
        h = np.zeros((len(reference), 2, 3))
        h[:, 0, 0] = h[:, 1, 1] = 1
        h[:, 0, 2], h[:, 1, 2] = -rotated[:, 1], rotated[:, 0]
        cov = rotation @ ref_cov @ rotation.T + test_cov if weighted else np.broadcast_to(np.eye(2), (len(reference), 2, 2))
        weights = np.linalg.inv(cov)
        normal = np.einsum('nai,nab,nbj->ij', h, weights, h)
        if np.linalg.matrix_rank(normal) < 3:
            raise ValueError('rigid rotation is unidentifiable: degenerate source geometry')
        parameter_cov = np.linalg.inv(normal)
        residual = test-rotated-translation
        update = parameter_cov @ np.einsum('nai,nab,nb->i', h, weights, residual)
        translation += update[:2]
        angle += update[2]
        if np.max(np.abs(update)) < 1.e-10:
            break
    if not weighted:
        dof = 2*len(reference)-3
        parameter_cov *= np.sum(residual**2)/dof if dof > 0 else np.nan
    return angle, translation, parameter_cov, h, cov


def fit_rigid(reference, test, config=None, ref_cov=None, test_cov=None, minimum=10):
    cfg = config or UncertaintyConfig()
    reference, test = np.asarray(reference, float), np.asarray(test, float)
    if reference.shape != test.shape or reference.ndim != 2 or reference.shape[1] != 2:
        raise ValueError('reference_xy and test_xy must be matching N x 2 arrays')
    finite = np.isfinite(reference).all(axis=1) & np.isfinite(test).all(axis=1)
    if finite.sum() < minimum:
        raise ValueError(f'at least {minimum} matches are required')
    weighted = False
    eligible = finite.copy()
    if ref_cov is not None and test_cov is not None:
        ref_cov, test_cov = np.asarray(ref_cov, float), np.asarray(test_cov, float)
        valid_cov = np.isfinite(ref_cov).all(axis=(1, 2)) & np.isfinite(test_cov).all(axis=(1, 2))
        indices = np.flatnonzero(valid_cov)
        valid_cov[indices] &= (np.linalg.eigvalsh(ref_cov[indices]) > 0).all(axis=1) & (np.linalg.eigvalsh(test_cov[indices]) > 0).all(axis=1)
        if (finite & valid_cov).sum() >= minimum:
            eligible &= valid_cov
            weighted = True
    def solve(indices):
        return solve_rigid(reference[indices], test[indices], ref_cov[indices] if weighted else None, test_cov[indices] if weighted else None)
    keep = eligible.copy()
    for _ in range(6):
        angle, translation, _, _, _ = solve(keep)
        residual = np.linalg.norm(test-reference @ _rotation(angle).T-translation, axis=1)
        median = np.median(residual[keep])
        sigma = 1.4826*np.median(np.abs(residual[keep]-median))
        new = eligible & (residual <= median+4*max(sigma, 1.e-6))
        if np.array_equal(new, keep) or new.sum() < minimum:
            break
        keep = new
    angle, translation, cov, h, obs_cov = solve(keep)
    residual_xy = test-reference @ _rotation(angle).T-translation
    residual = np.linalg.norm(residual_xy, axis=1)
    factors = np.array([1., 1., 180/np.pi])
    formal_cov = factors[:, None]*cov*factors[None, :]
    parameters = np.r_[translation, np.rad2deg(angle)]
    rng = np.random.default_rng(cfg.random_seed)
    source_ids = np.flatnonzero(keep)
    draws, postmedian, postp95 = [], [], []
    measured_median, measured_p95 = [], []
    if weighted:
        root_ref = np.linalg.cholesky(ref_cov[keep])
        root_test = np.linalg.cholesky(test_cov[keep])
    for _ in range(cfg.bootstrap_samples):
        idx = rng.choice(source_ids, size=len(source_ids), replace=True)
        try:
            a, t, _, _, _ = solve(idx)
        except (ValueError, np.linalg.LinAlgError):
            continue
        a = angle + np.arctan2(np.sin(a-angle), np.cos(a-angle))
        draws.append(np.r_[t, np.rad2deg(a)])
        r = np.linalg.norm(test[idx]-reference[idx] @ _rotation(a).T-t, axis=1)
        postmedian.append(np.median(r))
        postp95.append(np.percentile(r, 95))
        if weighted:
            # Perturb the fixed inlier sample and refit: preserve dependence of
            # residuals on the same measurements used to estimate the transform.
            rr = reference[keep] + np.einsum('nij,nj->ni', root_ref, rng.normal(size=(len(source_ids),2)))
            tt = test[keep] + np.einsum('nij,nj->ni', root_test, rng.normal(size=(len(source_ids),2)))
            try:
                ma, mt, _, _, _ = solve_rigid(rr, tt, ref_cov[keep], test_cov[keep])
                mr = np.linalg.norm(tt-rr @ _rotation(ma).T-mt, axis=1)
                measured_median.append(np.median(mr))
                measured_p95.append(np.percentile(mr,95))
            except (ValueError, np.linalg.LinAlgError):
                continue
    draws = np.asarray(draws).reshape(-1, 3)
    names = ('translation_east_arcsec', 'translation_north_arcsec', 'rotation_deg')
    uncertainties = {}
    for j, name in enumerate(names):
        sampled = interval(draws[:, j], parameters[j], cfg, n=len(source_ids))
        formal = measurement(parameters[j], np.sqrt(max(0, formal_cov[j, j])), cfg,
                             method='absolute GLS positional covariance' if weighted else 'residual-scaled least-squares covariance')
        # Envelope retains measurement floor without adding noise twice to the
        # source bootstrap (which already includes observed measurement scatter).
        adopted = sampled.copy()
        if weighted:
            adopted['ci_low'] = min(formal['ci_low'], sampled['ci_low']) if sampled['ci_low'] is not None else formal['ci_low']
            adopted['ci_high'] = max(formal['ci_high'], sampled['ci_high']) if sampled['ci_high'] is not None else formal['ci_high']
            adopted['method'] = 'envelope of formal GLS and paired-source bootstrap intervals'
            adopted['standard_error'] = max(formal['standard_error'], sampled['standard_error'] or 0)
        adopted['formal'] = formal
        adopted['sampling'] = sampled
        uncertainties[name] = adopted
    uncertainties['postfit_median_arcsec'] = interval(postmedian, np.median(residual[keep]), cfg, n=len(source_ids), method='paired-source bootstrap with rigid refit, fixed inlier selection')
    uncertainties['postfit_p95_arcsec'] = interval(postp95, np.percentile(residual[keep], 95), cfg, n=len(source_ids), method='paired-source bootstrap with rigid refit, fixed inlier selection')
    for name, simulations in [('postfit_median_arcsec', measured_median), ('postfit_p95_arcsec', measured_p95)]:
        propagated = interval(simulations, uncertainties[name]['value'], cfg,
            method='Gaussian positional-error Monte Carlo with rigid refit; fixed inlier sample', n=len(source_ids))
        if not weighted:
            propagated['reason'] = 'measurement propagation unavailable for unweighted fit'
        uncertainties[name] = envelope(uncertainties[name], propagated)
    # Residuals of fitted observations are correlated with fitted parameters.
    # C_residual = C_observation - H Cov(beta) H^T, not a sum.
    residual_cov = np.full((len(reference), 2, 2), np.nan)
    if weighted:
        residual_cov[keep] = obs_cov-h @ cov @ h.swapaxes(1, 2)
    sd, low, high, _ = radial_errors(residual_xy, residual_cov, cfg)
    return dict(rotation_deg=float(parameters[2]), rotation_uncertainty_deg=uncertainties['rotation_deg']['standard_error'],
                translation_east_arcsec=float(parameters[0]), translation_north_arcsec=float(parameters[1]),
                postfit_median_arcsec=float(np.median(residual[keep])), postfit_p95_arcsec=float(np.percentile(residual[keep], 95)),
                n_inliers=int(keep.sum()), uncertainties=uncertainties,
                parameter_covariance=formal_cov.tolist(),
                bootstrap_parameter_covariance=np.cov(draws, rowvar=False).tolist() if len(draws) > 1 else None,
                covariance_parameter_order=list(names),
                fit_method='positional-covariance weighted GLS' if weighted else 'unweighted rigid fit; sampling uncertainty only',
                inlier_mask=keep.tolist(), residual_arcsec=residual.tolist(), residual_err_arcsec=sd.tolist(),
                residual_ci_low_arcsec=low.tolist(), residual_ci_high_arcsec=high.tolist())
