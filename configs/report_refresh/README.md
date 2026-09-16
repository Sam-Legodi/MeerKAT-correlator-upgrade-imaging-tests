# Refresh the eight SDP comparison report sets

These drafts prepare future pipeline runs. No PDF/DOCX was created or modified
while preparing them. The image inputs were recovered from the existing metrics
JSON files, rather than guessing filenames from renamed configuration directories.
`input_manifest.json` records that provenance. All required local images and
catalogues exist. Output directories are new, leaving previous reports intact.

## Important capability boundary

The current `meerkat-ci report` is **image-only**. It does not consume visibility
CSVs, and its text explicitly excludes visibility assessment. No supported YAML
setting combines target/gaincal visibility results into that document.

The `vis` command produces separate visibility PDF/DOCX reports and diagnostic
plots. Whole-scan bootstrap comparison intervals are written to
`compare_dataset_summary.csv` and `compare_uncertainties.json`; they are **not yet
plotted or inserted into the visibility PDF/DOCX**. Visibility oscillation gates
still use point estimates; the 1-sigma image acceptance logic does not apply to
these visibility gates. Thus these configs refresh all supported products, but a
combined imaging-plus-visibility report, with visibility error bars and updated
visibility decisions, needs a subsequent pipeline code change. Do not regard
running `vis` followed by `report` as achieving that integration.

## Inputs

| Product | Required inputs |
|---|---|
| Consolidated imaging report | Matched FITS table for each band; corresponding CMC1/CMC2 images with celestial WCS and BMAJ/BMIN; explicit `extra.positions` |
| Detailed position report | Same matched table/images, including Xposn/Yposn columns |
| Detailed flux report | MFS, low-band and high-band matched FITS tables |
| Rebuild matches | CMC1/CMC2 PyBDSF FITS catalogues: RA, DEC, peak/total flux, pixel positions; E_RA/E_DEC, E_Total_flux/E_Peak_flux for formal errors; S_Code for isolation cuts |
| Rebuild catalogues, if absent | Existing continuum images and PyBDSF installation (`bdsf`) |
| Visibility QA | Complete, readable MS/MMS with its subtables, flags, scan/time, antenna/baseline, spectral-window and correlation metadata; populated corrected visibility values |
| Visibility comparison | Same-field CMC1 and CMC2 corrected MSs; two or more independent scans per dataset for scan-cluster sampling intervals |
| Optional contextual figures | Existing observation `mfimage_frequency_assessment/*_01_combined_plane_diagnostic.png` and `*_02_subband_common_scale.png`; calibration PDFs under `calreports/` are contextual, not quantitative visibility inputs |

The regenerated image report uses 5000 samples, seed 20260911 and central 68.27%
intervals. Missing formal errors stay unavailable; available sampling intervals
remain separate. Position and flux sampling are conditional on selection.
MFS PB-corrected images and non-PB subbands are compared only like-for-like.

## Config inventory

Each epoch has one `imaging/<epoch>.yaml` and two visibility drafts:
`visibility/<epoch>_J2147-8132.yaml` and
`visibility/<epoch>_J1619-8418.yaml` (24 YAML files total).

| CMC2 epoch/mode | CMC1 reference used in existing image reports |
|---|---|
| 1785059546_4kL | 1785941476_4kL |
| 1785065779_1kL | 1785941476_4kL |
| 1784443274_8kL | 1785941476_4kL |
| 1784437450_32kL | 1785941476_4kL |
| 1785659477_4kS4 | 1787557158_4kS4 |
| 1785652878_8kS4 | 1787557158_4kS4 |
| 1785567678_32kS4 | 1787557158_4kS4 |
| 1788421870_1kS4 | 1787557158_4kS4 |

Visibility reference epochs deliberately match image provenance. Older example
visibility configs use **1787034555** for S4, and one has no active reference.
Those are different comparisons. Confirm corrected MS availability for
**1787557158**, including both fields; if unavailable, state the different
visibility reference explicitly rather than silently replacing the image reference.

## Local image refresh commands

Run from the repository root in its Python environment. Editable installation
ensures the commands use the current checkout:

```bash
python -m pip install -e .
python configs/report_refresh/check_inputs.py imaging configs/report_refresh/imaging/*.yaml
```

The current catalogue inputs already exist, so source finding is unnecessary.
Rebuild matches into fresh output directories, then produce all report variants:

```bash
for cfg in configs/report_refresh/imaging/*.yaml; do
  meerkat-ci --config "$cfg" xm || break
  meerkat-ci --config "$cfg" pos || break
  meerkat-ci --config "$cfg" flux || break
  meerkat-ci --config "$cfg" report || break
done
```

For just the consolidated report, run `xm` and `report` for each config. For one
observation, for example:

```bash
meerkat-ci --config configs/report_refresh/imaging/1788421870_1kS4.yaml xm
meerkat-ci --config configs/report_refresh/imaging/1788421870_1kS4.yaml report
```

If catalogues have been deleted, install PyBDSF in the active environment and run
`meerkat-ci --config <imaging-config> src` before `xm`. Existing low/high images
are explicitly listed, so no cube slicing is required. If those images themselves
are absent, recover the exact historical products or use the original image
configuration's `low_high_slice` stage on its supplied MFImage cuboids; do not
substitute arbitrary frequency slices without updating the comparison provenance.

Use the individual stages above, **not `images` or `all`**: these drafts already
wire three explicit images per dataset; automatic image wiring would infer extra
pairs. `all` also invokes calibration. Neither recalibration nor reimaging is
needed for this SDP-product assessment.

Default new output base: `<repo>/data/report_refresh/<epoch>/`. Consolidated
reports/JSON go under `reports/imaging_verification/`; detailed position/flux
reports are alongside fresh matches under
`processed/Sky-CrossMatches/crossmatched-positions/` and `crossmatched-fluxes/`.
CLI `reports/produced_reports.log` lists exact filenames, including suffixes.

## Remote visibility preparation and commands

The remote root in each visibility draft comes from the existing server config:
`/home/slegodi/corr_upgrade_tests/2026`. Confirm it. Replace `reference_ms` and
`test_ms` with the **actual complete paths** to the separately exported field MSs;
placeholder basenames are labels, not claims that these files exist. Fill in
`extra.refresh_metadata.corrected_data_provenance` and `spectral_setup` for your
records. These metadata keys document provenance; they do not alter analysis.

The `vis` wrapper ignores `casa.field` and `casa.datacolumn`. It analyzes all MS
rows and prefers **CORRECTED_DATA**, then **DATA**, then MODEL_DATA. Therefore:

- Each configured MS must contain only its intended field's rows. Unused other
  FIELD-table entries are fine. Keep gaincal and target runs separate.
- For an SDP export with corrected values in DATA, confirm DATA is calibrated and
  no stale CORRECTED_DATA takes precedence. Never use a model-only MS as evidence.
- Do not run `cal`: the inputs are already SDP corrected. Recalibrating would
  change the dataset being evaluated. A filename containing `cal` alone does not
  establish SDP provenance; use the export history/logs.
- No splitting is needed for the separate field MSs you already have. If only a
  multi-field parent exists, export each field with its verified corrected column
  into a new MS before using the wrapper. For explicit field selection without
  exports, the underlying module supports `--field` (example below).

Transfer this config directory to the remote checkout (these new drafts have not
been committed/pushed automatically):

```bash
rsync -av configs/report_refresh/ USER@SERVER:/PATH/TO/REPO/configs/report_refresh/
```

On the server, activate the working analysis environment with `python-casacore`,
then from the updated repository root:

```bash
python -m pip install -e '.[radio]'
python configs/report_refresh/check_inputs.py visibility configs/report_refresh/visibility/*.yaml
# Proceed only when all needed paths/provenance/field checks pass.
for cfg in configs/report_refresh/visibility/*.yaml; do
  meerkat-ci --config "$cfg" vis || break
done
```

For one epoch's two fields:

```bash
meerkat-ci --config configs/report_refresh/visibility/1788421870_1kS4_J2147-8132.yaml vis
meerkat-ci --config configs/report_refresh/visibility/1788421870_1kS4_J1619-8418.yaml vis
```

No explicit uncertainty mapping is placed in visibility YAML: the wrapper does
not forward one. The comparison helper uses its existing 5000/68.27%/20260911
defaults. For specialized spectral settings or multi-field MS selection, invoke
the underlying module instead; it forwards `--field` to both datasets:

```bash
python -m meerkat_corr_imaging.vis_amp_analyze \
  --ms /ACTUAL/CMC2_SDP_corrected.ms \
  --ms-ref /ACTUAL/CMC1_SDP_corrected.ms \
  --field J2147-8132 \
  --outdir /NEW/OUTPUT/target --exact-outdir
```

Repeat with `J1619-8418` and a distinct output directory for gaincal. Default
Savitzky–Golay window/order are 51/3 and flag-derived channel exclusion is 20%.
Different channel widths mean 51 channels span different bandwidths: document
native-resolution comparisons, or establish matched physical frequency coverage,
channel widths and time averaging before claiming like-for-like noise/oscillation
comparisons. The module supports `--window`, `--order`, `--chan-min`, `--chan-max`
and `--rfi-flag-threshold`; YAML `frequency_ranges` does not restrict visibility QA.
Do not assume scan numbers identify the same time across independent epochs.
Target structure and gain-calibrator flux/variability also affect amplitudes.

Per-field default output base:
`<remote_root>/report_refresh/<epoch>/<field>/`. Actual per-MS reports, CSVs and
PNGs go under `interim/<MS-basename-without-final-suffix>/`; the two MS basenames
must differ. Comparisons are in the test MS output directory. Existing comparison
orientation is A=test, B=reference, hence B−A is **reference minus test**.

## Collect the visibility evidence locally

Copy the newly generated visibility results to a separate local directory:

```bash
rsync -av USER@SERVER:/home/slegodi/corr_upgrade_tests/2026/report_refresh/ \
  data/report_refresh_visibility/
```

Retain the configs/provenance, `scan_averaged_amp_stats.csv`,
`acceptance_summary.csv`, `flagging_by_channel.csv`, the other flag/scan/baseline
CSV tables, `compare_dataset_summary.csv`, `compare_uncertainties.json`, PNGs,
and produced visibility reports for **both fields and reference/test datasets**.
PDFs alone are not sufficient numerical inputs for future report integration.
The current image-report builder has no ingestion hook for this directory.

`check_inputs.py` is read-only and does not generate reports. Successful preflight
checks readability/schema/field selection, not the truth of SDP correction
provenance or full spectral comparability. Remote MS checks cannot be completed
on this local machine until the actual server paths and files are supplied.
