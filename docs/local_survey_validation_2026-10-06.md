# Local catalogue validation — 6 October 2026

Validation used the existing `meerkat-ci` Python 3.10 environment (Astropy 6.1.7).
All catalogue inputs were read locally from the resolved
`GPU Correlator Commisioning/catalogues` directory. No downloaded catalogues
were added to the repository, and no imaging or source finding was run.
All generated match tables, JSON diagnostics and audits were written to fresh
temporary directories. Existing science products were not overwritten.

## Real catalogue ingestion

| Input | Format/HDU | Rows | Valid coordinates | Distinct parent IDs | Nearest row to J2147 centre |
| --- | --- | ---: | ---: | ---: | ---: |
| summs212.fit | FITS ASCII / 1 | 211,050 | 211,050 | no parent ID mapped | 0.0000056° |
| racs-2deg-nearJ2147.fit | FITS ASCII / 1 | 177 | 177 | 133 | 0.7996190° |
| racs-10deg-nearJ2147.fit | FITS ASCII / 1 | 21,934 | 21,934 | 18,190 | 0.7996190° |
| racs.fit | FITS ASCII / 1 | 2,462,693 | 2,462,693 | 2,123,638 | 0.7996190° |

The completed `racs.fit` appeared during validation. Its extension is
`J_other_PASA_38_58_gausscut` and has the same inspected RAJ2000/DEJ2000,
GID/ID Gaussian schema as the two regional downloads. It is RACS-low DR1,
887.5 MHz, and is not low3, mid, or high. It was fully ingested once and matched
against all five selected MeerKAT catalogues through the shared matcher;
each yielded zero pairs within 10 arcseconds. No full-file subset was substituted
for this ingestion check.

SUMSS numeric coordinates and St/e_St, Sp/e_Sp were verified against the local
schema. Peak-brightness semantics were verified against the linked catalogue
documentation in the [usage guide](local_survey_crossmatching.md). Original FITS
unit labels were retained and the peak-brightness qualification recorded.

## Real MFS associations

The legacy PyBDSF launcher exports `catalog_type="gaul"`. Consequently, all five
files named `-source-cat.fits` are Gaussian-component tables, not distinct-source
tables. The validated profile preserves Gaus_id and Source_id without discarding
rows. The table below separates Gaussian row counts from parent-source counts.
Matching used a 10 arcsecond radius with all candidates retained. SUMSS
coordinates used FK5 J2000; MeerKAT/RACS coordinates used ICRS.

| MeerKAT product | Input/valid components | Distinct sources | Maximum extent | SUMSS pairs / matched components | SUMSS matched sources | SUMSS ambiguous components / sources | SUMSS unmatched components / sources | RACS pairs (2°, 10°, full) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| CMC1 1785941476 4kL | 508 | 432 | 0.770992° | 77 | 47 | 53 / 23 | 431 / 385 | 0 / 0 / 0 |
| GPU 1784437450 32kL | 697 | 614 | 0.772445° | 90 | 50 | 67 / 27 | 607 / 564 | 0 / 0 / 0 |
| GPU 1784443274 8kL | 621 | 537 | 0.772466° | 85 | 52 | 64 / 31 | 536 / 485 | 0 / 0 / 0 |
| GPU 1785059546 4kL | 250 | 203 | 0.780653° | 68 | 44 | 38 / 14 | 182 / 159 | 0 / 0 / 0 |
| GPU 1785065779 1kL | 447 | 385 | 0.784520° | 62 | 47 | 26 / 11 | 385 / 338 | 0 / 0 / 0 |

All MeerKAT coordinates were valid. For the RACS comparisons, every valid
MeerKAT row/parent source was unmatched and none had ambiguous positional
candidates. The small regional query, larger query and full catalogue all have
the same central row gap. These facts do not independently establish survey
noncoverage or failed MeerKAT astrometry; coverage remains `unknown` in outputs.
SUMSS matches remain positional candidates, with many ambiguous associations.
No secure counterpart selection or cross-frequency flux analysis was performed.

FITS CRVAL3 spectral-reference frequencies were read from the actual MFS images:

| Product | FITS frequency (MHz) |
| --- | ---: |
| CMC1 1785941476 | 1283.895507812 |
| GPU 1784437450 | 1283.986938477 |
| GPU 1784443274 | 1283.947753906 |
| GPU 1785059546 | 1283.059570312 |
| GPU 1785065779 | 1283.582031250 |

These are labelled spectral reference frequencies, not weighted effective MFS
frequencies. No nominal low/high midpoint was substituted.

## Checks and reproduction

The repository suite passed **295 tests**, including synthetic binary/ASCII FITS
and VOTable tables, multi-table selection, angular units, masked/invalid
coordinates, RA wraparound, empty matches, shared/ambiguous counterparts,
legacy nearest-candidate collision behavior, source/component identities,
parent-source counts, non-overwrite behavior, audited missing inputs,
INITIAL qualification, and preservation/deduplication through pipeline wiring.
Legacy two- and three-path wrapper jobs were also run through the real packaged
matcher on temporary synthetic files, verifying the explicit `1.0 arcsec` fix
and retained MeerKAT measurement enrichment. `git diff --check` and compilation
of the Python sources passed.

Run the suite in the repository's configured environment:

```sh
PYTHONPATH=src python -m pytest -q
```

The validation helper described in the [usage guide](local_survey_crossmatching.md)
reproduces the 15 regional jobs in a new directory and records input SHA-256
hashes, positions relative to J2147, full job specifications, diagnostics and
audits. The full export can be passed to `--racs-wide` for a slower full-file run.
The four original catalogue inputs and all science image/catalogue paths remain
external to this repository.

## Remaining real-data requirements

`RACS-low3_INITIAL_sources.fits`, `RACS-mid1_sources.xml` and
`RACS-high_sources.xml` were absent from the local catalogue directory. Their
real table selection, coordinate/ID schemas, units and release provenance have
not been validated. VOTable ingestion and custom INITIAL profile behavior were
validated using explicitly synthetic fixtures. New profiles must be based on
inspection of the actual files; low3 INITIAL results must remain preliminary.

The CMC1 1786962370 image exists, but its expected PyBDSF catalogue in the wide-L
configuration is currently absent. The example therefore uses the available
CMC1 1785941476 catalogue. No catalogue was regenerated to fill this requirement.
Low/high-band explicit job configuration is supported and tested with synthetic
fixtures, but real low/high survey associations were outside this MFS validation.
