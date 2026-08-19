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

Run each comparison from the repository root in three inspectable stages:

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

The flux wrapper remains disabled until source finding and cross-matching have
produced separate low-band, high-band and MFS match tables. The extracted
low/high products are explicitly non-PB, while the delivered MFS comparison
images are PB-corrected; downstream analyses must preserve that distinction.
