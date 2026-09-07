# Imaging verification implementation log

Date: 2026-09-01

## Scope implemented

- Added a consolidated, image-domain continuum verification report to the `meerkat_corr_imaging.cli` pipeline as the `report` step.
- Added the report step to the `images` and `all` workflows.
- Generated one five-page draft DOCX report and one machine-readable JSON metrics sidecar for each of the 32kL, 8kL, 4kL and 1kL test observations.
- Used the supplied primary writing sample for voice: short declarative statements, evidence before interpretation, numerical results in the main text, and explicit limits on what the evidence supports.
- Retained the required `draft_` prefix for every generated DOCX basename.

## Evidence reviewed

### Calibration reports

| Test | PDFs reviewed | Visual finding | Report disposition |
|---|---:|---|---|
| 1784437450_32kL | Four spectral-chunk reports | Broad persistent-RFI masks and structured channel flagging in both polarisations. No antenna was identified as rejected for all channels or all gain-solution times. | Concern; no visibility-derived aggregate fraction was calculated. |
| 1784443274_8kL | Four spectral-chunk reports | Broad persistent-RFI masks and structured channel flagging in both polarisations. No antenna was identified as rejected for all channels or all gain-solution times. | Concern; no visibility-derived aggregate fraction was calculated. |
| 1785059546_4kL | One report | Broad, highly structured flagging. Antenna m005 was rejected for all channels in the second bandpass solution and for all times in the gain solution. | Concern. |
| 1785065779_1kL | One report | Broad persistent-RFI blocks and structured channel flagging in both polarisations. Antenna m027 was rejected for all channels in the second bandpass solution; no antenna was rejected for all gain-solution times. | Concern. |

The calibration PDFs were treated as contextual visual evidence. The `<20%` flagging and `<1%` amplitude-oscillation acceptance thresholds were not scored because the required visibility-domain arrays were unavailable.

### Image products

| Test | Visual finding | Report disposition |
|---|---|---|
| 1784437450_32kL | Concentric deconvolution residuals around the phase-centre source in the MFS image and every inspected subband; no subband exceeded the 10 mJy/beam severe-corruption diagnostic. | Concern. |
| 1784443274_8kL | Concentric deconvolution residuals around the phase-centre source in the MFS image and every inspected subband; no subband exceeded the 10 mJy/beam severe-corruption diagnostic. | Concern. |
| 1785059546_4kL | Strong concentric and grid-like residual structure. Five of nine independent subbands (1157, 1232, 1411, 1543 and 1613 MHz) were classified as severely corrupted; outer-field RMS values span 11.3--239 mJy/beam. | Concern. |
| 1785065779_1kL | Concentric deconvolution residuals around the phase-centre source in the MFS image and every inspected subband; no subband exceeded the 10 mJy/beam severe-corruption diagnostic. | Concern. |

## Measurement definitions implemented

- Position acceptance: Pass when the raw radial source-separation p95 is below 1 arcsec. Median, maximum, fitted rigid rotation, translation and post-fit residuals are reported as secondary diagnostics.
- Flux acceptance: Pass when the absolute deviation of the robust median integrated-flux test/reference ratio from unity is below 5%. Peak-flux ratio, p95 absolute deviation and the fraction of sources beyond 5% are secondary diagnostics.
- RMS acceptance: Pass when the robust test-image RMS divided by the corresponding CMC1 RMS is below 1.2. Both values use `1.4826 × MAD` in the same 0.25--0.50 degree phase-centred annulus.
- Catalogue selection: isolated PyBDSF components (`S_Code=S`) with peak-flux signal-to-noise ratio at least 10 in both catalogues.
- Astrometric model: a scale-fixed east/north translation plus one global rigid rotation, with robust clipping and bootstrap rotation uncertainty.
- Cross-match search gate: increased from 1 to 5 arcsec. This prevents the source-selection gate from truncating the p95 statistic at the 1-arcsec decision boundary.
- Primary-beam state: MFS uses PB-corrected reference and test products. Low- and high-band comparisons use non-PB reference and test products. Unlike correction states are never combined.

## Acceptance results

| Test | MFS: position / flux / RMS | Low: position / flux / RMS | High: position / flux / RMS | Overall |
|---|---|---|---|---|
| 32kL | Concern / Pass / Pass | Concern / Concern / Pass | Concern / Pass / Concern | Concern |
| 8kL | Concern / Concern / Concern | Concern / Pass / Pass | Concern / Pass / Concern | Concern |
| 4kL | Concern / Concern / Concern | Concern / Pass / Concern | Not assessed / Not assessed / Concern | Concern |
| 1kL | Concern / Concern / Concern | Concern / Concern / Pass | Concern / Concern / Concern | Concern |

The 4kL high-band test image produced only two PyBDSF catalogue rows and no quality-selected matches. Catalogue-derived astrometry and flux metrics are therefore not assessed for that product. Its image RMS ratio is 893.11 and remains a Concern.

## Code and configuration patched

- `src/meerkat_corr_imaging/verification_report.py`: report analysis, decisions, plots, DOCX generation and JSON sidecar.
- `src/meerkat_corr_imaging/steps/step8_verification_report.py`: pipeline wrapper.
- `src/meerkat_corr_imaging/cli.py`: new `report` command and report integration into image workflows.
- `tests/test_verification_report.py`: rigid-transform recovery, insufficient-match and FITS PB-state parsing tests.
- `README.md`: documented Step 8, image/all workflow order, report outputs and the `draft_` DOCX convention.
- Four run configurations under `GPU Correlator Commisioning/mci_pipeline_processing/config_files`: 5-arcsec cross-match gate, calibration-report findings, image-inspection findings and report provenance.
- Existing output-path and wording patches were retained. Their details are in `REPORT_WORDING_PATCH_LOG.md`.

The final repository audit also hardened `PBCOR` parsing so string-valued FITS booleans such as `"F"` are not treated as truthy by Python. Explicit non-PB filename forms are recognised before the PB filename fallback. This did not change the current report results.

## Checklist limits retained in the reports

- No raw visibilities: auto/cross flagging fractions, time/channel/correlation breakdowns, scan-averaged visibility mean/RMS, and amplitude oscillations are not assessed.
- No per-scan images: time/scan-dependent image-domain position and relative-flux comparisons are not assessed.
- No new PB correction in this cycle: low/high results are valid only as non-PB like-for-like comparisons.
- The planned visibility-spectrum method is recorded as a 51-channel, third-order Savitzky--Golay model over an RFI-free flag-derived channel mask.

## Verification performed

- Rebuilt all four reports through `meerkat_corr_imaging.cli report`.
- Successfully compiled the modified package.
- Ran the full test suite successfully: 45 tests passed.
- Rendered all four final DOCX files to PDF and PNG.
- Visually inspected all 20 final pages. Each report contains five populated pages with no clipping, unintended blank pages or broken figures.
- Reopened each DOCX package and parsed each JSON sidecar. Audited all 12 pipeline DOCX files: every basename starts with `draft_`, and none contains the targeted informal or writer-directed wording.
