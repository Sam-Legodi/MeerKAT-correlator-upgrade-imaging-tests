# CMC1–GPU/CMC2 frequency-matched astrometry imaging

This is an **opt-in CASA 6 imaging mode** for already SDP-corrected MSs. It images
both reference and test. The default `experiment_mode: calibrator_and_target`
images `J1619-8418` (gain calibrator) and `J2147-8132` (target) separately.
Explicit `experiment_mode: target_only` images only the configured target.
Each field can have its own reference/test MS, as in the
samples; a multi-field MS may also be supplied for both field entries. Names
must match exactly and identify one FIELD row in each MS. Scans are discovered
independently per field and MS; they are never inferred from filenames or shared
between roles. Calibration is not run. Corrected samples in these exports are in
`DATA`, so use `datacolumn: data`, not `corrected`.

Configs without `casa.paired_astrometry.enabled: true` retain their existing
Step 3 behaviour, nominal low/high selections, MS list, names, and handoffs.
A legacy `reference.ms_paths` input alone does not enable paired mode. Paired
mode uses its explicit field input mapping and launches one CASA process for the
whole experiment. It does not automatically wire these new products into the
older source-finding/cross-match pipeline; inspect their manifests before choosing
images to compare. Use `cal`, not `all` or `images`, for this imaging experiment.

## Explicit target-only experiment

Set `casa.paired_astrometry.experiment_mode: target_only` and provide exactly one
`fields` entry with `kind: target`, its actual FIELD name, `reference_ms`, and
`test_ms`. Remove gain-calibrator entries from this mapping. The default remains
`calibrator_and_target`, which requires one gain calibrator and one target.
Unknown modes, duplicate/missing targets, or calibrator entries in target-only
mode fail validation before MS path resolution. There is no automatic fallback
when a calibrator has no usable channels. Legacy `reference.ms_paths` / test MS
lists are unused in either paired mode.

Target-only mode never requires, opens, surveys, calibrates, images or modifies a
calibrator MS. It preserves the same positive-finite-weight/FLAG criteria,
reference-defined native quartiles, gap-preserving selectors and common interval
overlap rules. Defaults remain 80% occupancy and 90% common coverage. Zero weights
remain unusable; this option does not repair an export or substitute weights.

Use `configs/paired_astrometry/l_band_target_only.yaml` or
`configs/paired_astrometry/s4_target_only.yaml`. Replace placeholder CBIDs/modes
and verify the two target MS paths. `field: J2147-8132` selects the actual FIELD
name; the independent `target_file_tag: J2147` builds filenames such as
`1785059546_J2147_SDPflags+cal.4kL.ms`. Both samples use already corrected `DATA`
and a separate `_target_only` output directory. Do not change existing MS data,
weights or flags to make eligibility pass.

The samples set `casa.executable: /opt/casa-6.6.5-31-py3.10.el8/bin/casa`
for bruce. Change this YAML setting for other servers; no shell export is needed.
An existing `$CASA` environment variable takes precedence, so `unset CASA` once
if it points to an old installation. Configs without `executable` retain the
legacy `casa` command on PATH.

From the repository root inside tmux/screen, run the configured executable. The new samples default to `dry_run: true`, so this exact command
performs the selection-only survey and writes the planned target image matrix:

```bash
meerkat-ci --config configs/paired_astrometry/l_band_target_only.yaml cal
```

After inspecting the plan, set `casa.paired_astrometry.dry_run: false` in that
same config. The exact imaging command is:

```bash
meerkat-ci --config configs/paired_astrometry/l_band_target_only.yaml cal
```

For S4 use `configs/paired_astrometry/s4_target_only.yaml` in both commands.
These commands require the package installed in the orchestration Python
environment; the equivalent `PYTHONPATH="$PWD/src" python -m
meerkat_corr_imaging.cli` invocation also works. Neither real-MS selection nor
CASA imaging was performed locally while implementing this mode.

Target-only results provide **no calibrator-image check of phase transfer or
astrometric systematics**. A successful target run is not evidence that the
gain-calibrator export or its calibration corrections were valid. Independently
verify calibration provenance and phase transfer before interpreting target
offsets. Existing frequency/scan, beam, UV-coverage and input-weight effective
frequency limitations below continue to apply.

## Remote run

Install this package in the orchestration Python (3.10+) environment, e.g.
`python -m pip install -e .`. Install/use the remote server's **CASA 6** executable;
its bundled Python must provide NumPy and Astropy. No casacore dependency is
needed for this step. The new helper uses only Python 3 standard libraries, NumPy,
CASA table/tasks, and Astropy for FITS provenance. CASA 5 remains a legacy option
but is rejected for paired experiments. No exact CASA 6 patch version has been
validated against real data here; every image records the running task version
and its fully expanded public task defaults. Use the same CASA installation for
both roles and repeat runs.

Copy and edit `configs/paired_astrometry/l_band.yaml` or
`configs/paired_astrometry/s4.yaml`. Replace **every PLACEHOLDER** CBID/mode and
check all four resolved `.ms` paths against the server layout. Paths use the
requested `remote_root`, `data_root`, `analysis_root`, `ref_stem`, `test_stem`
conventions; the calibrator has additional gain-field stems. The S4 sample uses
`S4_SDPflags+cal` as its analysis folder. The band label never defines frequency
selection: flags and weights in the MS do.

From the repository root, inside a tmux/screen session, run:

```bash
PYTHONPATH="$PWD/src" \
  python -m meerkat_corr_imaging.cli cal --config configs/paired_astrometry/l_band.yaml
```

For S4 change only the config argument to `configs/paired_astrometry/s4.yaml`.
The CLI uses the existing finite CASA batch runner and propagates failure status.
Increase `casa.timeout_seconds` if the complete matrix needs longer than the
sample's 86400 seconds. No local MS dry run or expensive imaging was performed
during implementation.

Before expensive imaging on the remote server, set
`casa.paired_astrometry.dry_run: true` and run the same command. The read-only MS
survey writes `selection_plan.json` and `native_channel_diagnostics.json` and
prints the matrix, channel counts, effective frequencies, and coverage. It does
not call `tclean` or export FITS, and does not replace completed image manifests.
Inspect these plans, then set `dry_run: false` and rerun.

The compatibility script also accepts `--pair-config` with a **resolved JSON
mapping of the paired_astrometry section**, not the master YAML. For example,
after resolving path variables with the normal config loader:

```python
import json
from meerkat_corr_imaging.config import load_config
cfg = load_config("configs/paired_astrometry/l_band.yaml")
with open("pair.json", "w") as handle:
    json.dump(cfg.casa.paired_astrometry, handle)
```

```bash
PYTHONPATH="$PWD/src" /path/to/casa-6/bin/casa --nologger --log2term --nogui \
  -c scripts/tclean_two_bands.py --pair-config pair.json
```

Do not pass positional MSs, `--scans`, legacy field exclusions, or
`force_calibrate: true` with paired mode. The pipeline clears inherited MS and
paired-config overrides to avoid accidentally selecting a different input.

## Eligibility and rounding

Only MAIN cross-correlation rows (`ANTENNA1 != ANTENNA2`) for the exact field ID
are counted; `antenna='*&*'` applies the same baseline policy in `tclean`.
Rows are queried by DATA_DESC_ID to handle variable channel/correlation counts
across SPWs, read in bounded chunks, and summed across all DDIDs for each physical
SPW/channel. `chunk_rows` and `max_chunk_bytes` cap the row count and estimated
working-array bytes (a single native row is the irreducible minimum). No complete
MAIN visibility/flag/weight column is loaded. SPW metadata uses per-row getcell.
The CASA selector follows the documented
[antenna and channel selection syntax](https://casadocs.readthedocs.io/en/v6.6.1/notebooks/visibility_data_selection.html).

* Select both XX/YY or both RR/LL; fail for a populated DDID without a complete
  parallel pair. Read the POLARIZATION table, not fixed correlation indices.
* Flags are CASA's boolean `FLAG OR FLAG_ROW`. To conservatively follow strict
  `stokes='I'`, all stored correlations must be unflagged in that row/channel.
  Both selected parallel hands must also have finite DATA and positive finite
  weights. Cross-hand weights do not enter the occupancy denominator.
* Use a defined, nonempty `WEIGHT_SPECTRUM` cell; otherwise broadcast the row's
  `WEIGHT` across channels. Zero, negative, infinite, and NaN spectrum weights
  are unusable and **never trigger fallback** to positive row weights.
* The per-channel denominator is the number of present selected cross-baseline
  rows times the number of selected parallel hands, including flagged and
  zero-weight opportunities. An opportunity contributes both parallel hands
  only when the strict conditions above hold. Occupancy is their joint usable
  count divided by that denominator. Keep occupancy >= `occupancy_threshold`
  (default 0.8), with at least one joint usable opportunity.
* Undefined required DATA/FLAG/WEIGHT cells are reported as absent and unusable
  within the present-row denominator. SPWs with no selected rows have zero
  denominator and are ineligible. Unscheduled rows are not invented; this metric
  does not measure loss relative to an expected observing schedule.
* Raw integer `BITFLAG` and `BITFLAG_ROW`, when present, are summarized as observed
  value histograms. They do not replace CASA FLAG or alter selection. Individual
  SDP reason bits may have been collapsed during export and cannot be recovered
  from boolean FLAG. Policy and selected correlations are recorded explicitly.
* Invalid/zero-width or nonfinite CHAN_FREQ/CHAN_WIDTH metadata is ineligible.

CASA's stricter Stokes-I flag handling is described in
[the tclean documentation](https://casadocs.readthedocs.io/en/v6.6.1/_modules/casatasks/imaging/tclean.html).
The input survey is deliberately conservative; it is not a measurement of the
final UV-grid contribution under every CASA minor release.

For each MS and field, sort usable physical channels by
`(CHAN_FREQ, SPW_ID, channel_ID)`. With N channels, use half-open Python slices:

| Product | Slice |
| --- | --- |
| low | `[0 : N//4]` |
| middle | `[(3*N)//8 : (5*N)//8]` |
| high | `[N-N//4 : N]` |

All boundaries use integer floor; no banker rounding. N must be >=4. For N=8,
indices are low=0–1, middle=3–4, high=6–7. For N=5 they are low=0, middle=1–2,
high=4. Small/odd N can yield unequal counts and a floor bias; the exact IDs are
always recorded. Gaps remain gaps in both interval diagnostics and CASA selectors
(e.g. `0:1~3;8;12~14,2:0~7`). Quartiles count channels, not MHz; different widths
and overlapping SPWs can therefore yield different native MHz spans. Both MSs'
native quartile diagnostics are saved even though the paired images use the
reference-defined bands.

## Paired interval selection

For each field, take the reference's native low/middle/high channel list as the
requested band, and all eligible test channels as candidates. Each channel's
interval is `centre +/- abs(CHAN_WIDTH)/2`, using its stored frequency frame.
Intervals are unions: overlapping SPWs never double count coverage. Different
frequency frames fail rather than silently comparing LSRK and TOPO values.
No Doppler/frame transformation is performed; TOPO values at different epochs
are compared numerically in their stored frame.

Iteratively remove channels in either MS with less than
`min_channel_overlap` (default 0.9) of their own interval covered by the other's
retained union, until stable. This excludes non-overlapping channels and protects
flag gaps. Require common intersection bandwidth to cover at least
`min_common_coverage` (default 0.9) of **each retained MS union and the original
requested reference band**. Fail clearly if any band or full-band selection
fails, before any imaging. This prevents retaining a tiny accidental intersection
and claiming adequate coverage. Test native quartiles are diagnostics, not the
paired band definition. Relaxing these thresholds changes the experiment contract.

Keep whole channels. No clipping, interpolation, or common-grid resampling is
performed. Identical interval unions have exact frequency support even when the
channel widths differ. Otherwise report the intersection, both actual native
unions, all unmatched retained edges, excluded reference intervals, and all three
coverage fractions. A 90% overlap **does not mean identical support**. Inspect
unmatched edges and effective frequency differences when interpreting astrometry.

The full-band common selection is calculated separately from all eligible
reference channels and the test candidate channels. Use it for every scan in both
roles. Per-scan actual support can differ because flags, weights, and data vary
with time. Each image records its contributing channel IDs, common coverage, and
`per_scan_support_incomplete`; all per-scan products are explicitly labelled as
independent scans with a common requested selection, not guaranteed equal actual
support or temporal pairing. No reference/test scan correspondence is inferred.

For each field and each MS there are three all-scan images and one full-band image
per discovered scan. `calibrator_and_target` totals
`12 + sum(scan counts across the four field/MS inputs)`;
`target_only` totals `6 + N_reference_target_scans + N_test_target_scans`.
All are single-term Stokes-I MFS. No calibrator and target scans are combined.

## Provenance, reruns, and parameters

Output names are `<output_dir>/images/<role>_<field-and-hash>_<band>` or
`..._fullband_scan<N>`, with `.image`, `.fits`, `_qa.txt`, and `.manifest.json`.
A shared reference FIELD PHASE_DIR sets both roles' imaging phase centre; original
MS directions and frames are preserved in the manifest. Moving, polynomial, and
ephemeris fields are rejected. Restoring beams remain independently fitted unless
a common `restoringbeam` is explicitly configured. The output includes:

* `native_channel_diagnostics.json`: both native channel inventories, denominator,
  flag/weight/absence counts, correlations, scans, native quartile IDs and fractions.
* `selection_plan.json`: every planned image and complete immutable input/task
  contract, including paired selections and experiment mode. This is also the
  dry-run deliverable.
* `images/*.manifest.json`: running/complete/failed status, input contract and
  fingerprint, timestamps, errors, FITS celestial WCS geometric centre/reference
  pixel/reference coordinates, restoring beam, and all resolved tclean parameters.
* `run_manifest.json`: the completed run's matrix, including failures and reused
  products. A failed task makes the batch exit unsuccessfully.
* `planning_failure.json`: field/product context and selection error if overlap
  planning fails. Native diagnostics for all configured inputs are preserved in this
  case; no imaging has started.

The plan, native diagnostics and each per-image contract record `experiment_mode`,
`calibrator_imaging_qa_status`, and `interpretation_limitations`. In target-only
mode the status is explicitly `not_performed_target_only`. This describes omitted
calibrator imaging/QA, not the target's calibration history. Combined mode records
calibrator products as requested; their individual processing status determines
whether they completed. The mode is part of every image fingerprint: switching
modes cannot reuse incompatible target products in the same output directory.
Old manifests lacking this provenance also fail the fingerprint check; use a
fresh output directory. No mode change deletes existing products.

Requested channel records come from the aggregate field selection. The separate
`selected_channels_with_actual_data` records and `selection.contributed_channel_ids`
come from the actual field/scan flags, finite data and positive input weights.
All requested channels are retained in the record, including zero-contribution
channels in individual scans. Frequency span is the outer interval endpoints;
actual bandwidth sums interval unions so gaps are not filled. The effective
frequency uses summed positive **input** weights on jointly valid opportunities,
not Briggs/robust or UV-grid weights. Final MFS gridding has no per-native-channel
output census; `contribution_status` explicitly identifies this input-level
estimate and unmeasured gridding contribution. Field-wide weight-source counts
and raw flag-bit histograms are survey summaries, while per-channel/per-image
flag and weight counts are selection-specific.

Reuse requires matching SHA256 contract fingerprints, an intact completed image
manifest and both image/FITS products. A missing/failed manifest or different
inputs, flags/weights, channel selection, phase centre, CASA version, or parameters
fails preflight for the whole run. There is no existence-only reuse or automatic
deletion/restart in paired mode. Choose a new `output_dir` (e.g. a run suffix), or
explicitly remove a stale product's entire family after review. The fingerprint
includes a streamed digest of selected DATA/flags/weights/TIME/UVW/baselines and
metadata, plus conservative MS file sizes/mtimes (excluding table.lock). Harmless
MS rewrites may invalidate reuse. Task-generated products must stay outside the
input MS. Visibility content hashing is independent of read chunk size. Existing
mask/voltage-pattern files are also fingerprinted, so changing their contents at
the same path does not authorize reuse. Keep MSs immutable throughout the run
and avoid simultaneous runs in
the same output directory. Dry-run plans never authorize reuse themselves.

Change `occupancy_threshold`, `min_channel_overlap`, `min_common_coverage`, and
`tclean` parameters in YAML, not Python constants. `cell`, `imsize`, `niter`,
`threshold`, `gridder`, `wprojplanes`, `weighting`, `robust`, `uvtaper`, `mask`,
`restoringbeam`, and other public CASA task parameters can be overridden. Selection
parameters, Stokes/specmode, start models/outlier fields, and task side effects are
controlled and cannot be overridden. `mtmfs` is rejected because this product
contract is single-term MFS. Unspecified task parameters are read from the running
CASA signature and recorded. The samples' cell sizes are starting points, not an
astrometric accuracy claim; select finer sampling and an appropriate common beam
for the planned experiment.

Subarcsecond interpretation still depends on calibration phase transfer,
ionosphere, observing epoch and baseline/UV coverage, source structure and spectral
index, primary beam/off-axis effects, restoring beam, S/N, deconvolution, and frame
conventions. Matching frequency overlap alone does not isolate correlator-induced
position shifts. A common image phase centre also does not calibrate away physical
or metadata differences. Inspect per-scan support and effective frequencies, verify
that the source is compact, and use consistent fitting and uncertainty methods.
The implementation was tested with mocks and synthetic FITS, not CASA imaging on
a representative MS. Validate the remote selection plan and actual CASA products
before drawing an astrometric conclusion.

## Automated verification

Run focused tests with `python -m pytest -q tests/test_paired_astrometry.py` and
the full repository suite with `python -m pytest -q`, using the normal project
environment (`pip install -e '.[dev]'`). Tests mock CASA tables/tasks and cover
eligibility, absence, weights, raw bits, variable SPWs, frequency intersections,
irregular gaps, all rounding boundaries, field/scan isolation, compatibility
dispatch, image provenance, failure status, chunk-independent input identity,
external-mask changes and stale-output refusal. Synthetic FITS tests verify WCS
centre and restoring-beam extraction. These tests do not validate CASA execution
or claim an observed astrometric accuracy.
Target-only tests also verify mode conflicts/defaults, independent scan counts,
frequency coverage, omitted-QA provenance, zero-weight rejection, mode-change
stale-output detection, compatibility dispatch, and a guarded run that refuses
any calibrator filesystem access or calibration call.
