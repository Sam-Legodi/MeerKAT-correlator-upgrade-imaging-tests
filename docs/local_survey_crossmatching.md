# Local SUMSS and RACS cross-matching

The `xm` step accepts optional `extra.survey_xmatch_jobs` and
`extra.survey_xmatches` in addition to legacy two/three-path
`extra.xmatch_pairs`. Existing pairs continue to use the global `xmatch` settings.
Survey jobs have independent settings and produce positional associations and
counts, without automatic cross-frequency flux calibration or position analysis.
No catalogue downloads, source finding, or imaging are needed for standalone `xm`.

The opt-in [local MFS example](../configs/survey_xmatch/local_mfs.yaml) selects an
available CMC1 catalogue and four GPU L-band catalogues against local SUMSS and
RACS-low DR1. Run from the repository with its installed Python environment:

```sh
meerkat-ci --config configs/survey_xmatch/local_mfs.yaml xm
```

Set `paths.sky_xmatches_dir` and `paths.reports_dir` to a fresh output location
before each run. Survey outputs refuse overwrites. No survey jobs are added to
existing configurations automatically. An enabled, explicitly requested missing,
incomplete or invalid catalogue fails the step with an audited diagnostic;
`enabled: false` skips an optional job/template. No wildcard download discovery is
performed; `.crdownload`, `.part`, and `.tmp` inputs are rejected.
The `GPU_Correlator_Commissioning` directory is a symlink to
`GPU Correlator Commisioning`; inputs are resolved before job deduplication.

## Explicit jobs and profiles

```yaml
extra:
  survey_xmatch_jobs:
    - input1:
        profile: pybdsf_source
        path: /path/to/CMC1-mfs-source-cat.fits
        image: /path/to/CMC1-mfs-image.fits
        hdu: 1
      input2:
        profile: racs_low_dr1_gausscut
        path: /path/to/racs-2deg-nearJ2147.fit
        hdu: 1
      output: /path/to/new/Sky-CrossMatches/CMC1_X_RACS_DR1.fits
      radius: 10 arcsec
      association_mode: all-candidates
      coverage: unknown
```

Every catalogue resolves to a structured specification: path, optional `format`
(`fits` or `votable`), FITS `hdu` (index or extension name), VOTable `table_id`
(zero-based index or TABLE ID), `survey`, `release`, `grain` (`source` or
`component`), `ra_col`, `dec_col`, `frame`, `coordinate_unit`, optional
`frequency_hz`, `image`, `source_id`, `component_id`, `required_columns`,
`qualification`, and `documentation`. Profiles supply defaults; each input can
override them. Multi-table files require explicit selection. Both FITS ASCII and
binary tables and VOTable XML are supported. The original columns and units are
retained with `_1` and `_2` suffixes. No flux fields are renamed or summed.

Built-in profiles:

| Profile | Coordinates/frame | Grain and identifiers | Survey frequency |
| --- | --- | --- | --- |
| `pybdsf_gaul` | RA, DEC; ICRS | component; Gaus_id, parent Source_id | image metadata if supplied, otherwise unknown |
| `pybdsf_source` | RA, DEC; ICRS | actual source rows; unique Source_id | image metadata if supplied, otherwise unknown |
| `sumss_v21r` | _RAJ2000, _DEJ2000; FK5 J2000 | source rows | 843 MHz |
| `racs_low_dr1_gausscut` | RAJ2000, DEJ2000; ICRS | component; GID, parent ID | 887.5 MHz |

SUMSS `St/e_St` are integrated flux/error in mJy. `Sp/e_Sp` represent peak
brightness/error in mJy/beam, despite the local VizieR FITS columns' `mJy` label;
the profile records this qualification without changing the original unit.
See the [SUMSS catalogue documentation](https://cdsarc.cds.unistra.fr/viz-bin/ReadMe/VIII/81B?format=html).
RACS Gaussian components keep their positions, shapes, `GID` and parent `ID`.
Repeated `Ftot` values belong to the parent source and must not be summed across
components. Input row counts and distinct parent-source counts are separate.

The existing PyBDSF launcher writes `catalog_type="gaul"` even though filenames
end in `-source-cat.fits`. Generated jobs therefore use `pybdsf_gaul`, preserving
`Gaus_id` and parent `Source_id`. The inspected local inputs have repeated source
IDs and must not be labelled source-grain. Use `pybdsf_source` only for actual
source tables; repeated source IDs in a source-grain profile are rejected.

Use `extra.catalogue_profiles` for newly inspected schemas. A custom profile is
a mapping containing the specification fields above except `path`; reference it
by name under either input's `profile`. No built-in coordinate or ID schema is
assumed for missing newer RACS files. An INITIAL release or filename automatically
adds a preliminary measurement qualification. Explicit configurations should also
include the catalogue's documented qualifications and frequency.

## Pipeline generation and low/high-band inputs

`extra.survey_xmatches` holds templates with `name`, `catalogue`, `radius`,
`association_mode`, optional `input_profile` (default `pybdsf_gaul`), optional `targets` (exact target names), and optional `bands`
(default `[mfs]`). It generates one job per selected image product from
`reference.images` and `tests[].images`. Product catalogues are located using the
existing PyBDSF naming rules. More than one arbitrary image per target is labelled
`image_1`, `image_2`, etc.; it is not silently assumed to be MFS.

Standalone `xm` and the `images` command resolve the same jobs. Image pipeline
wiring retains explicit jobs and settings and adds deterministic generated jobs;
equivalent jobs are deduplicated using resolved input paths and settings, including
radius in arcseconds. A conflicting output path for different jobs is rejected.
Generated survey outputs include a digest of settings. Survey tables are kept
outside the automatic MeerKAT–MeerKAT `positions` and `flux` handoffs.

For already produced low/high-band source catalogues, add explicit
`survey_xmatch_jobs` as above with their catalogue and image paths. No slice step
is required for explicit inputs. Alternatively, when using existing configured
slice products (`low_high_slice.enabled` and `add_to_source_finding`), add
`bands: [mfs, low, high]` to a template. Run only `xm` if you intend to reuse those
catalogues. `images` also runs the image/source-finding pipeline as documented
elsewhere and should only be used when that work is intended.

The configured wide-L ranges 898–1,000 and 1,460–1,700 MHz have nominal
midpoints 949 and 1,580 MHz. These are not substituted for observing or effective
image frequencies. For a supplied image, frequency provenance reads a FITS FREQ
axis and its declared units, or explicit FREQ/CFREQ header cards. RESTFRQ and
filename guesses are not observing-frequency substitutes. The spectral reference
frequency in an MFS header is reported as such; no weighted effective frequency
is inferred. An explicit frequency can be configured for inputs without usable
image metadata, and its origin remains `configured`.

Desired future comparisons require the actual RACS-low3 943.5 MHz
`RACS-low3_INITIAL_sources.fits`, RACS-mid 1,367.5 MHz
`RACS-mid1_sources.xml` (for approximately 1.2 GHz MeerKAT products), and RACS-high
1,655.5 MHz `RACS-high_sources.xml`. Inspect their actual schemas and table IDs
before defining profiles. The low3 INITIAL measurements must remain preliminary.
The existing 887.5 MHz DR1 gausscut downloads are distinct from these releases.
The completed `racs.fit` was inspected and is also that DR1 gausscut table,
with 2,462,693 Gaussian rows; it can use the same inspected profile. It does
not fill the J2147 catalogue gap. Use the smaller local export for routine
field matching to avoid repeatedly loading the full ASCII table.

## Associations, diagnostics and coverage

`all-candidates` retains every pair within a positive finite angular `radius`.
`one-to-one` preserves the legacy algorithm: choose each input-1 row's nearest
input-2 row, sort by separation, and accept unused rows. A row lost to a collision
is not reassigned to its second-nearest candidate. Both modes evaluate all nearby
candidates when marking ambiguity, so a selected one-to-one pair may be ambiguous.
These are positional candidates; neither mode establishes secure counterparts.

Each output FITS table resides in `Sky-CrossMatches` and includes:

- Original columns and units, `sep_arcsec` with an angular unit, and zero-based
  `input_row_1`/`input_row_2` identifiers for the original input table rows.
- `candidate_count_1`, `candidate_count_2`, and `ambiguous` (either row has multiple
  candidates), including shared counterparts and repeated parent IDs.
- Input, valid-coordinate, invalid-coordinate, matched, unmatched and ambiguous
  row counts per input in metadata. Ambiguous counts include candidates discarded
  by one-to-one selection. Counts use the declared source/component grain. Where parent IDs are available,
  separate valid, matched, unmatched and ambiguous distinct-source counts are
  also reported; duplicate component geometry remains in the output.
- Catalogue provenance, frame, release, grain, observing frequency and origin,
  radius, association mode and qualifications in FITS metadata/HISTORY and an
  adjacent `.diagnostics.json` with full paths and structured specifications.

Coordinates must have angular units (unitless columns use the explicitly
specified `coordinate_unit`), finite RA in [0, 360] degrees and DEC in [-90, 90].
Invalid and masked coordinates are excluded and counted separately from valid
unmatched rows. Matching is spherical and supports RA wraparound. Empty matches
produce valid zero-row tables and diagnostics with the original column schema.

Coverage defaults to `unknown`. `covered` or `absent` requires an independent
`coverage_evidence` reference in the job; the supplied assertion is recorded.
The pipeline does not infer a footprint from the absence of catalogue rows. Both
local RACS DR1 downloads have no entries inside about 0.7996 degrees of J2147's
centre, while the inspected GPU MFS rows extend only about 0.78 degrees. A larger
query radius does not fill that gap. Zero associations neither establish failed
MeerKAT astrometry nor justify increasing the match radius to force counterparts.
Cross-frequency flux calibration, spectral corrections and secure association
selection need a separate validated analysis.

## Reproduce local validation

The validation helper only ingests and matches existing catalogues, relocates all
outputs to a fresh directory, records input hashes and the J2147 radial extents,
and produces `validation.json` plus the normal FITS/diagnostic/audit files:

```sh
PYTHONPATH=src python scripts/validate_local_survey_matches.py \
  --config configs/survey_xmatch/local_mfs.yaml \
  --racs-wide '/Users/samuel/SARAO/MK+/Correlator/GPU_Correlator_Commissioning/catalogues/racs-10deg-nearJ2147.fit'
```

To repeat full-file validation, pass `catalogues/racs.fit` to `--racs-wide`.
`--output-root` must name a new directory; otherwise a temporary directory is
created. See [the recorded validation results](local_survey_validation_2026-10-06.md)
for the real inputs, synthetic coverage and remaining requirements.

## SUMSS astrometry configurations

`configs/image_only/{l_wide,l_n10732k,u,s4}.yaml` enable SUMSS matching for both
the configured CMC1 MFS reference and every GPU MFS test. All eight
`configs/report_refresh/imaging/*.yaml` explicitly select their original CMC1
and GPU MFS catalogues, so their multiple-image lists do not accidentally omit
SUMSS matching or include low/high products. The radius is 10 arcseconds with
all candidates retained and marked for ambiguity; existing MeerKAT-pair settings
remain independent.

For only SUMSS cross-matches, the ready-to-use
[`sumss_astrometry.yaml`](../configs/survey_xmatch/sumss_astrometry.yaml) selects
the available CMC1 1785941476 4kwide reference used by the L-band report-refresh
configs and all four GPU MFS products. It contains no legacy pair jobs or RACS
jobs and writes to a dedicated `data/sumss_astrometry` directory:

```sh
meerkat-ci --config configs/survey_xmatch/sumss_astrometry.yaml xm
```

The catalogue tables contain MeerKAT `RA_1/DEC_1` in ICRS, SUMSS
`_RAJ2000_2/_DEJ2000_2` in FK5 J2000, separations, original position errors,
component/source identifiers and ambiguity flags. To interpret signed offsets,
transform both coordinates to one frame and specify whether the sign means
MeerKAT minus SUMSS or SUMSS minus MeerKAT. Separations and coordinate errors
remain distinct quantities. Retained ambiguous candidates and Gaussian
components are not independent secure counterparts; use a stated association
selection before estimating an astrometric offset.

Configuration enables catalogue matching for astrometry comparison; it does
not add SUMSS plots or offset statistics to the existing `pos`/verification
reports. Their schema requires two MeerKAT images and pixel/beam fields absent
from SUMSS. No substitute image or unsupported schema is passed to those reports.
The CMC1 1786962370 reference remains configured in `l_wide.yaml`; its catalogue
must exist for that comparison, and the dedicated config avoids that missing
input using the explicitly documented 1785941476 reference.
