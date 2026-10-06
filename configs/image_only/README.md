# Image-only correlator comparison configs

These configs use the PB-corrected continuum images for source finding,
catalogue cross-matching and astrometry. Visibility QA and CASA imaging are
intentionally skipped because no MeasurementSets are supplied.

The non-PB MFImage cuboids contain 7–10 independent subband images. The
`low_high_slice` step selects one maximum-overlap plane for each configured
low/high frequency range and writes non-PB 2-D FITS products. The cuboids do
not contain the original 1k/4k/8k/32k visibility channels, so correlator data
cannot be re-averaged to 1k channels from these FITS files. The configs compare
the delivered continuum image products as they are. This tests end-product
imaging performance, not an isolated channel-count effect.

Run a complete image-domain comparison from the repository root with:

```bash
PYTHONPATH=src python -m meerkat_corr_imaging.cli images --config CONFIG.yaml
```

The `images` mode runs `low_high_slice -> src -> xm -> pos -> flux`. It passes
the resolved low/high images into source finding, the deterministic PyBDSF
catalogues into cross-matching, and the cross-match outputs into position and
flux analysis. Explicit MFS pairs and position entries in these configs are
retained; missing low/high handoffs are added automatically.

For inspection or recovery, every stage remains independently runnable:

```bash
PYTHONPATH=src python -m meerkat_corr_imaging.cli low_high_slice --config CONFIG.yaml
PYTHONPATH=src python -m meerkat_corr_imaging.cli src --config CONFIG.yaml
PYTHONPATH=src python -m meerkat_corr_imaging.cli xm  --config CONFIG.yaml
PYTHONPATH=src python -m meerkat_corr_imaging.cli pos --config CONFIG.yaml
```

`l_n10732k.yaml`, `l_wide.yaml`, `u.yaml` and `s4.yaml` have usable CMC1
references. `s1_source_find_only.yaml` has no valid full-band CMC1 reference
and must not be used for cross-matching against the S4 image: S1 covers about
2007–2801 MHz, while S4 covers about 2676–3455 MHz.

Flux analysis runs automatically only when a test has separate low-band,
high-band and MFS match tables. The extracted low/high products are explicitly
non-PB, while the delivered MFS comparison images are PB-corrected; downstream
interpretation must preserve that distinction.

## SUMSS catalogue astrometry

The paired `l_wide.yaml`, `l_n10732k.yaml`, `u.yaml` and `s4.yaml` configs now
include independent SUMSS V2.1r matches for the CMC1 reference and every GPU test,
using existing MFS catalogues and a 10 arcsecond radius. `xm` resolves these jobs
without imaging/source finding. Low/high-band matching remains opt-in. Gaussian
components, parent IDs, original coordinate units and ambiguity flags are retained.

For a SUMSS-only run with an available CMC1 reference and all four wide-L GPUs:

```bash
meerkat-ci --config configs/survey_xmatch/sumss_astrometry.yaml xm
```

That config writes separate results under
`data/sumss_astrometry/processed/Sky-CrossMatches`, without executing the legacy
CMC1–GPU pair jobs. Choose a fresh output directory for subsequent runs. The
paired configs retain their configured references; a missing reference catalogue
is reported as a failure rather than replaced with a different observation.

The SUMSS tables provide MeerKAT `RA_1/DEC_1`, SUMSS `_RAJ2000_2/_DEJ2000_2`,
`sep_arcsec`, stable input row IDs and ambiguity flags for absolute catalogue
astrometry. SUMSS uses FK5 J2000, while MeerKAT uses ICRS; signed east/north offsets
must be evaluated in a common frame, with an explicit sign convention. These
survey tables are not inputs to the existing `pos` reports, which require two
MeerKAT image/catalogue schemas. See the
[SUMSS astrometry configuration guide](../../docs/local_survey_crossmatching.md#sumss-astrometry-configurations)
for the scope and association qualifications.

Read-only preflight on 6 October 2026 found missing existing MFS catalogues for
CMC1 1786962370 (`l_wide`), GPU 1783266356 (`l_n10732k`), CMC1 1785990374 and
both GPU 1784017048/1784022889 (`u`), and CMC1 1787034555 (`s4`). Their configured
SUMSS jobs will fail clearly until those catalogues are supplied. The dedicated
SUMSS-only configuration and all eight report-refresh MFS pairs have readable
catalogue inputs; no missing catalogue was regenerated during this config patch.
