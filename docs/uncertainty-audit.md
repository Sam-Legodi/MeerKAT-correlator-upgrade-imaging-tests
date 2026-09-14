# Uncertainty audit and implementation contract

Audited 2026-09-11 before editing analysis code. Reference output: existing
`draft_1788421870_1kS4_imaging_verification_1788421870_1kS4.pdf` (six pages),
under the S4 processing reports directory. The supplied attachment is instructions,
not a report. Its structure matches `verification_report.py`.

## Provenance and current calculations

- `cli.py`, `image_pipeline.py`, `steps/step[1-8]*`, `config.py`: orchestrate
  calibration, image slicing, source finding, matching, detailed analyses and
  consolidated reporting. `low_high_slice.py` selects image planes; it does not
  infer independent measurement errors. `pybdsf_srcfind.py` exports the fitted
  catalogue and RMS/residual/model images.
- `xmatch_pybdsf.py`: spherical nearest-neighbour associations with the configured
  gate; all input columns are copied with `_1`/`_2` suffixes. `sep_arcsec` is
  spherical separation. Gate truncation and association ambiguity are not
  represented by catalogue formal errors.
- Actual S4 GAUL tables have E_RA/E_DEC (deg), E_Total_flux (Jy), E_Peak_flux
  (Jy/beam), E_Xposn/E_Yposn (pix), shape/PA/deconvolved-shape errors, and
  E_Isl_Total_flux. No RA/DEC or cross-dataset covariance columns are present.
  Shape/pixel/island fields are preserved but not quantitative comparison metrics.
- Consolidated report: isolated S_Code=S components, peak/error >=10 on both
  sides, at least 10 matches. Raw separation median/p95/max use SkyCoord.
  Rigid fit uses phase-centred spherical offsets, SVD translation/rotation,
  iterative radial MAD clipping; only rotation had a 400-draw bootstrap SD.
  Postfit median/p95 use retained sources. Total/peak ratios require positive
  finite fluxes. Report includes median ratio, p95 absolute total-ratio deviation,
  fraction above 5%. RMS is 1.4826*MAD of finite annular pixels, strided by 3.
  Ratios and all aggregates except rotation were point estimates.
- `positions_analysis.py`: spherical east/north offsets, beam-normalized offsets,
  radius, direction, pixel field angles; tables give mean/std/median/p16/p84.
  Harmonic least squares gives quadrupole amplitude/phase and R-squared changes;
  circular summaries give resultant length/mean direction/Rayleigh p-value.
  Ellipse widths and displayed SDs describe scatter, not errors of the mean.
- `flux_analysis.py`: per-source test/ref and fractional differences; mean/std,
  median/p16/p84/NMAD; iterated through-origin gain and Huber line fits. Existing
  line covariance is residual-scaled and ignores reference-axis error; missing
  y errors incorrectly receive maximum weight. Per-scan bars use scatter despite
  a mean-uncertainty label. Detailed reports only store formatted summaries.
- `vis_amp_analyze.py`: exact observed flag counts, scan/baseline mean amplitudes,
  detrended spectral SD and SD/mean, group medians and reference/test differences.
  These have correlated channels/baselines, no independent measurement covariance
  in exported summaries. Exact flag census is not a binomial measurement sample.
- `calibration_checks.py`: CASA script quality gates on exact solution coverage
  and sampled model residuals; operational gates, not CMC1/CMC2 report metrics.
- Notebook-produced radial-noise/artefact panels are embedded external PNGs;
  report generation has no source pixel samples for those panel annotations.
  They cannot acquire measurement errors in the reporting layer. Existing user
  edits to the notebook must be preserved.

## Adopted methods

| Quantity | Method |
|---|---|
| Individual total/peak ratio and fractional difference | Independent analytic Jacobian; valid positive formal errors only |
| East/north and phase-centred coordinates | Numerical spherical-coordinate Jacobian, including cos(dec) and induced covariance |
| Radial separation / direction | Jacobian away from zero; Gaussian coordinate Monte Carlo for low-S/N radii; direction undefined at zero |
| Catalogue mean/SD/median/percentiles/NMAD | Reproducible paired-row percentile bootstrap; consolidated separation/flux metrics also retain fixed-sample measurement Monte Carlo and adopt an interval envelope |
| Fraction above threshold | Wilson interval, including boundary fractions |
| Rigid E/N translation and rotation | Iterated GLS using positional covariance when sufficient; formal full covariance plus paired-source bootstrap, circular unwrapping |
| Post-fit residual statistics | Refit each bootstrap sample; conditional on robust inlier set |
| Image MAD scale | Gaussian MAD asymptotic SE with effective independent beams, not pixel-count SD formula |
| RMS ratio | Analytic propagation of independent image-scale SEs |
| Gain and linear-fit parameters | Both flux axes contribute formal variance; paired-source bootstrap for conditional fitted sample |
| Sample maximum | Observed extremum; no population-tail CI from ordinary bootstrap; fixed-sample measurement Monte Carlo includes switching of the maximum source |
| Counts, beam/header values, configured limits | Exact bookkeeping/metadata; no invented measurement error |
| Rayleigh p-values / external panel annotations | Diagnostic statistics, not estimates with invented measurement errors |
| Visibility summary differences | Scan-cluster sampling intervals where independent scan groups exist; measurement covariance unavailable |

Default confidence is 0.6826894921370859 (Gaussian +/-1 sigma), 5000 resamples,
seed 20260911. Absolute CI endpoints remain numeric; asymmetric CIs are retained.
The user's latest instruction supersedes attachment item 11: existing numerical
limits stay fixed, but Pass requires the complete 1-sigma decision interval below
the upper limit (inside 0.95--1.05 for flux). Crossing intervals are Concern;
missing decision uncertainty is Not assessed. Store point-estimate status too.
Sampling intervals do not include common calibration systematics; no quadrature
addition of population scatter and measurement noise that would double count it.

## Numerical conventions and source references

- PyBDSF catalogue error meanings and units:
  https://pybdsf.readthedocs.io/en/stable/write_catalog.html
- Percentile bootstrap convention:
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.bootstrap.html
- Wilson intervals are checked against `scipy.stats.binomtest(...).proportion_ci`.
- Gaussian MAD derivation: with q=Phi^-1(0.75), the density of |X| at
  sigma*q is 2*phi(q)/sigma. The median quantile variance is 1/(4*N*f²),
  hence SE(1.4826*MAD) = RMS * 1.4826/(4*phi(q)*sqrt(N_eff)), approximately
  1.166*RMS/sqrt(N_eff). A Gaussian restoring-beam area is an approximation
  to noise correlation area, not a measurement of the complete noise covariance.
- Covariance matrices preserve parameter order and units in JSON. Position
  matrices use arcsec²; rigid parameters are E arcsec, N arcsec, rotation degrees.
- Formal and sampling intervals are stored separately. Their envelope is a
  conservative sensitivity rule and not an exact combined-coverage claim.
- Flux errors assume independent Gaussian fitted flux errors; low-S/N reference
  denominators make ratio distributions non-Gaussian. Asymmetric bounds are
  retained. Common flux-scale errors are not available from these products.
