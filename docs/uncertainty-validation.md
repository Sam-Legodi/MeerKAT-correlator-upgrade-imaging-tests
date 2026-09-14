# Validation completed 2026-09-14

## Results

- Full suite: **151 passed** (15.92 seconds), using the existing
  `/Users/samuel/miniforge3/envs/meerkat-ci/bin/python` environment and `PYTHONPATH=src`.
- Analytic ratio and spherical positional propagation; correlated coordinate
  inputs; masked/missing/zero errors; rigid weighted fit, covariance and
  measurement floors; bootstrap, Wilson, cluster resampling; formatting and
  strict JSON; RMS effective-beam calculation; acceptance boundary cases covered.
- Independent calculation of all 110 MFS total-flux ratio errors agrees to
  relative tolerance 1e-12 with the exported FITS values.
- Real data: CMC2 1788421870 1kS4 versus CMC1 1787557158 4kS4, MFS/low/high,
  using the historical report's exact image and cross-match paths. Current YAML
  references renamed directories; validation recovers paths from saved report
  provenance and redirects outputs to `data/uncertainty_validation/`.
- Consolidated, detailed positional and detailed flux reports ran with 5000
  resamples, confidence 0.6826894921370859 and seed 20260911.
- All 12 non-fit scalar/count fields per band match the previous report exactly.
  Rigid fits deliberately change because positional-covariance weighting replaces
  equal weighting. Historical reports and source images were not overwritten.
- MFS and high-band positional point verdicts were Pass; their 1-sigma decision
  bounds now cross 1 arcsec and give Concern. Low-band position remains Concern.
  Flux and RMS verdicts remain unchanged for this dataset.
- Consolidated DOCX rendered with the bundled `render_docx.py` (6 pages).
  Native pipeline PDFs render as consolidated 7 pages, positions 5 pages and
  flux 11 pages. Page overviews checked for overflow, missing content and table
  continuity; uncertainty tables remain readable and intervals are retained.
- `git diff --check` passes. The pre-existing notebook modifications remain
  untouched. No raw visibility/CASA rerun was performed for this image-only
  configuration; scan-cluster uncertainty is covered by deterministic tests.

## Files and review points

Core helpers: `uncertainty.py`, `astrometry.py`. Consumers: matcher, position,
flux and consolidated report modules; visibility comparison; configuration and
position/flux step wrappers. Documentation: README and uncertainty audit.
Tests: `tests/test_uncertainty.py` plus existing suite.

The audit maps each report metric to its inputs and method. Review assumptions:
independent catalogue/scan populations, Gaussian formal errors, fixed selection
and inlier cuts, fixed beam metadata, approximate Gaussian noise correlation
area. Interval envelopes are conservative sensitivity bounds, not exact combined
confidence coverage. No input cross-dataset covariance or calibration-systematic
error model exists.

Uncertainty remains unavailable where formal input errors are missing; sampling
intervals are separately retained. Population extrema, per-spectrum visibility
measurement errors without channel covariance, and externally generated notebook
panel annotations are not assigned invented errors. Exact counts and operational
flag/CASA gates retain their existing meaning. Per-scan positional overlay
ellipses remain descriptive scatter diagnostics; the main positional report
carries numerical propagated source errors and sampling summaries.

## Reproduction artifacts

`data/uncertainty_validation/validation_config.yaml` records the resolved inputs.
`old_new_comparison.json` records unchanged values, fit changes and both verdicts.
The consolidated report and matched-source FITS sidecars are under
`reports_dir/imaging_verification/`. Detailed reports and numerical summaries are
under `crossmatched-positions/` and `crossmatched-fluxes/`.
