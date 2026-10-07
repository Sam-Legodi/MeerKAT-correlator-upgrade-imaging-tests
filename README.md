# MeerKAT Correlator Upgrade — Imaging Tests

Config-driven workflows for MeerKAT CMC1 reference and GPU/CMC2 test observations: visibility QA, optional CASA calibration/imaging, SDP image preparation, PyBDSF source finding, catalogue matching, astrometry, flux and sensitivity comparisons, and draft DOCX/PDF reports. Start with fresh visibilities and the corresponding SDP archive images, observation metadata and calibration reports; reuse existing products when refreshing results.

---

## Table of Contents

* [Overview](#overview)
* [Repository Layout](#repository-layout)
* [Prerequisites](#prerequisites)
* [Installation](#installation)
* [Manual Data Download (Step 1)](#manual-data-download-step-1)
* [Fresh dataset: required workflow](#fresh-dataset-required-workflow)
* [How to Run](#how-to-run)

  * [1) Install the project locally (once per machine)](#1-install-the-project-locally-once-per-machine)
  * [2) Prepare a master config](#2-prepare-a-master-config)
  * [3) Run individual steps (surgical control)](#3-run-individual-steps-surgical-control)
    * [3.1 Visibility QA (Step 2)](#31-visibility-qa-step-2)
    * [3.2 Calibrate & Image with CASA (Step 3)](#32-calibrate--image-with-casa-step-3)
    * [3.3 Low/high cuboid slices (Step 4)](#33-lowhigh-cuboid-slices-step-4)
    * [Sensitivity metadata](#sensitivity-metadata-before-flux-and-report)
    * [3.4 Source finding with PyBDSF (Step 5)](#34-source-finding-with-pybdsf-step-5)
    * [3.5 Cross-matching catalogues (Step 6)](#35-cross-matching-catalogues-step-6)
    * [3.6 Astrometry (positions) analysis (Step 7a)](#36-astrometry-positions-analysis-step-7a)
    * [3.7 Flux analysis (Step 7b)](#37-flux-analysis-step-7b)
    * [3.8 Consolidated verification report (Step 8)](#38-consolidated-verification-report-step-8)
  * [4) Run grouped pipeline stages](#4-run-grouped-pipeline-stages)
  * [5) Where things go (default)](#5-where-things-go-default)
  * [6) Quick verification checklist](#6-quick-verification-checklist)
  * [7) Common gotchas (and fixes)](#7-common-gotchas-and-fixes)
  * [8) Commit your config and results?](#8-for-dev-purposes-commit-your-config-and-results)
  * [TL;DR sequence](#tldr-sequence)
* [Refresh reports from existing products](#refresh-reports-from-existing-products)
* [Detailed behaviour and advanced workflows](#detailed-behaviour-and-advanced-workflows)
  * [Calibration controls](#calibration-application-and-optional-imaging)
  * [Paired CASA astrometry imaging](#paired-cmc1gpucmc2-astrometry-imaging-casa-6)
  * [Remote sensitivity and noise reports](#remote-sensitivity-and-existing-image-noise-reports)
  * [Sensitivity reuse contract](#sensitivity-artifact-and-reuse-contract)
* [Reproducibility & Provenance](#reproducibility--provenance)
* [Examples & Tests](#examples--tests)
* [Legacy Scripts](#legacy-scripts)
* [Development Practices](#development-practices)
* [Contributing](#contributing)
* [Citation](#citation)

---

## Overview

**Goal:** Compare like-for-like reference/test datasets using corrected-visibility diagnostics, image astrometry and flux consistency, measured image RMS, and product-specific theoretical sensitivity. The final consolidated report combines available numerical evidence with contextual SDP calibration reports and image diagnostics. Missing evidence remains explicitly unassessed.

The recommended order is **prepare inputs → calibrate if required → corrected-visibility QA → prepare images → sensitivity → source finding → matching → astrometry → flux → consolidated report**. Calibration and image creation are conditional; sensitivity must precede the flux/report stages when those stages consume sensitivity JSON. See the [fresh-dataset workflow](#fresh-dataset-required-workflow) for executable routes and their configuration requirements.

---

## Repository Layout

```
MeerKAT-correlator-upgrade-imaging-tests/
├─ README.md
├─ .gitignore
├─ .gitattributes               # Git LFS pointers for large artifacts
├─ configs/
│  ├─ example_cluster.yaml
│  └─ example_local.yaml
├─ data/
│  ├─ raw/                      # manual archive downloads
│  ├─ interim/                  # intermediate calibration/QA outputs
│  ├─ processed/                # final images, catalogues, xmatches
│  └─ reports/                  # ready-to-share DOCX/PDF/PNG artefacts
├─ examples/                   # space for example datasets
├─ legacy_scripts/              # preserved historical utilities (read-only)
│  ├─ corr_imanalysis.py
│  ├─ image-analyser.py
│  ├─ image-rotation.py
│  ├─ image_analysis.py
│  ├─ py3gokatsdpimager.py
│  ├─ set_pyenv_settings.py
│  ├─ srcfind.py
│  ├─ srcfind.pyc
│  └─ stats.py
├─ src/
│  └─ meerkat_corr_imaging/
│     ├─ __init__.py
│     ├─ cli.py                 # optional CLI wrapper (`python -m meerkat_corr_imaging.cli`)
│     ├─ config.py              # config loading/validation helpers
│     ├─ flux_analysis.py       # flux comparison and report module
│     ├─ low_high_slice.py      # installable low/high cuboid extractor
│     ├─ output_paths.py        # common output naming helpers
│     ├─ positions_analysis.py  # astrometric analysis module
│     ├─ pybdsf_srcfind.py      # PyBDSF launcher module
│     ├─ standalone_xxyy_solve.py # packaged CASA calibration script
│     ├─ tclean_two_bands.py    # packaged CASA imaging script
│     ├─ vis_amp_analyze.py     # visibility QA module
│     ├─ sensitivity_workflow.py # selected-MS sensitivity JSON and noise workflow
│     ├─ survey_xmatch.py        # local SUMSS/RACS catalogue associations
│     ├─ output_audit.py         # final per-step output inventory
│     ├─ verification_report.py # consolidated image-domain verification report
│     ├─ xmatch_pybdsf.py       # catalogue matching module
│     └─ steps/
│        ├─ step1_archive_docs.py     # documentation helpers (no code execution)
│        ├─ step2_vis_analysis.py     # wraps vis_amp_analyze.py
│        ├─ step3_calibrate_image.py  # wraps CASA scripts
│        ├─ step4_low_high_slice.py   # resolves/reuses/extracts cuboid planes
│        ├─ step5_srcfind.py          # wraps pybdsf_srcfind.py
│        ├─ step6_xmatch.py           # wraps xmatch_pybdsf.py
│        ├─ step7_positions.py        # wraps positions analysis
│        ├─ step7_flux.py             # wraps flux notebook export
│        └─ step8_verification_report.py # builds the consolidated report
├─ scripts/                     # backward-compatible thin entrypoints
│  ├─ vis_amp_analyze.py
│  ├─ standalone_xxyy_solve.py
│  ├─ tclean_two_bands.py
│  ├─ pybdsf_srcfind.py
│  ├─ xmatch_pybdsf.py
│  ├─ positions_analysis.py
│  └─ flux_analysis.ipynb
├─ tests/
│  └─ test_low_high_slice.py
```

> **Note:** Pipeline implementations live inside the installable package.
> Files under `scripts/` are compatibility entry points only; pipeline steps
> do not resolve them relative to the current working directory.

---

## Prerequisites

* **Python**: 3.10+
* **System tools**:

  * **CASA** (external binary; accessible on `$PATH`) for calibration/imaging
  * **Git LFS** if you intend to track large artifacts (`*.fits`, `*.ms*`, `*.image`, `*.docx`)
* **Python libraries** (installed during setup): `astropy`, `numpy`, `pandas`, `matplotlib`, `python-docx`, `PyYAML`, etc.
* **Radio-specific Python libraries**: `python-casacore` and `PyBDSF` are not published as universal wheels. Install them separately (see [Radio-specific dependencies](#radio-specific-dependencies)) using Conda/Conda-Forge, or build from source if you already have casacore.

---

## Installation

Create a dedicated virtual environment and install this repository in editable mode so CLI wrappers pick up your edits immediately.

```bash
# Get the repo
# git clone https://github.com/Sam-Legodi/MeerKAT-correlator-upgrade-imaging-tests.git
# cd MeerKAT-correlator-upgrade-imaging-tests
python3 -m venv .venv
source .venv/bin/activate
pip install -e .
# include tooling like pytest/pre-commit:
# pip install -e ".[dev]"
# linux-only casacore bindings:
# pip install -e ".[radio]"
```
OR (If a default Conda environment is active. You can still follow the steps—just make sure you’re in the repo folder and pick one environment strategy (either keep using Conda, or use python -m venv; don’t mix them).):

```bash
conda create -n meerkat-ci python=3.10 -y
conda activate meerkat-ci
pip install -e .
# or pip install -e ".[dev]"
# (linux) pip install -e ".[radio]"

```

This repository now ships with a standards-compliant `pyproject.toml`, so `pip install -e .` (or `pip install -e ".[dev]"` if you also want pytest/pre-commit helpers) registers the package `meerkat_corr_imaging` under `src/`. After installation you can invoke either the module (`python -m meerkat_corr_imaging.cli ...`) or the shortcut console script `meerkat-ci ...`. Leave the environment with `deactivate` and reactivate it later with `source .venv/bin/activate`. CASA remains a separate binary; make sure `casa` is on your `PATH` or set `CASA=/path/to/casa`.

Optional but recommended tooling:

```bash
pre-commit install
git lfs install  # large images/catalogues are already listed in .gitattributes
```

### Radio-specific dependencies

The vis-analysis and source-finding steps rely on `python-casacore` and `PyBDSF`. These packages are distributed primarily via conda-forge/kernsuite channels, so they are not part of the default `pip install -e .` dependency list.

With Conda (recommended on macOS and Linux):

```bash
conda install -n meerkat-ci -c conda-forge python-casacore pybdsf
```

On Linux systems with casacore already available, you can alternatively use pip:

```bash
pip install -e ".[radio]"            # installs python-casacore from PyPI
```

Older BDSF requires numpy<2
If this fails: 
```bash
python -m pip install bdsf    
```
Then do refer to: https://github.com/lofar-astron/PyBDSF.

If you maintain your own builds, ensure both packages are on the environment `PYTHONPATH` before running `python -m meerkat_corr_imaging.cli ...`.

---

## Manual Data Download (Step 1)

Collect both the CMC1 reference and GPU/CMC2 test datasets. Keep observation/epoch IDs, field names, spectral setup, scans, channel selections, time averaging, flags, calibration history and image weighting with the files.

* **Visibilities:** complete MS/MMS directories with their subtables. Prefer SDP-corrected exports when evaluating the SDP products. Record whether calibrated values reside in `DATA` or `CORRECTED_DATA`; a filename alone is insufficient.
* **Images:** the SDP MFS continuum FITS images, non-PB MFImage cuboids and available low/high products. Retain PB-corrected and corresponding non-PB images where available; compare matching PB states.
* **Context:** observation metadata, SDP calibration-report PDFs and any existing `mfimage_frequency_assessment` diagnostic figures. These provide provenance and visual context; calibration PDFs do not replace numerical visibility QA.

For example, an SDP export with calibration applied (only the target and gaincal have valid SDP corrections) and visibilities averaged to 1k channels:

```bash
mvftoms.py 1784022889_sdp_l0.full.8kU.rdb -p "HH,VV" --applycal "all" --target "J1619-8418" --flags "static,cam,data_lost,ingest_rfi" --chanbin 8 -o 1784022889_J1619.8kU.ms && \
mvftoms.py 1784022889_sdp_l0.full.8kU.rdb -p "HH,VV" --applycal "all" --target "J2147-8132" --flags "static,cam,data_lost,ingest_rfi" --chanbin 8 -o 1784022889_J2147.8kU.ms
```

Record the actual export options and resulting column contents in your dataset notes, for example `data/raw/README.md`. Download/export remains manual. If visibilities are uncorrected, follow the calibration branch below before claiming corrected-data QA.

---

## Fresh dataset: required workflow

Run from the repository root in the installed Python environment, preferably inside `tmux` or `screen`. Use one reference/test pair per consolidated assessment; the report selects the first configured test and `noise-report` requires exactly one pair. Keep target and gain-calibrator visibility comparisons separate.

| Order | Action | Command / condition |
| --- | --- | --- |
| 1 | Collect visibilities, SDP images, metadata and calibration reports; configure paths and provenance | Manual preparation and YAML editing |
| 2 | Calibrate uncorrected visibilities, optionally make CASA images | `cal`; skip calibration for verified SDP-corrected data |
| 3 | Assess corrected reference/test visibilities | `vis` after calibration, or directly for corrected exports |
| 4 | Select existing SDP/CASA images and create missing low/high products | `low_high_slice` only when cuboid extraction is needed |
| 5 | Calculate product-specific sensitivity from the final visibility selections | `sensitivity` |
| 6 | Generate/reuse PyBDSF catalogues | `src` |
| 7 | Match reference/test catalogues, optionally local SUMSS/RACS catalogues | `xm` |
| 8 | Generate detailed astrometry report | `pos` |
| 9 | Generate detailed flux/noise report | `flux` |
| 10 | Generate consolidated verification report | `report`; CLI exports newly written DOCX files to PDF at each invocation's end |

`pos` and `flux` both consume matched tables; neither requires the other's DOCX. `report` consumes configured tables, images and numerical evidence, rather than combining those detailed DOCX files.

### Choose calibration and imaging independently

| Visibility / image situation | `extra.force_calibrate` | `casa.imaging_enabled` | Action |
| --- | --- | --- | --- |
| Already SDP corrected; use archive images | `false` | `false` | Skip `cal`; run `vis` and analyse supplied images |
| Already corrected; create CASA images | `false` | `true` | Run `cal` for imaging only, then `vis` |
| Uncorrected; calibrate and create CASA images | `true` | `true` | Run `cal`, then corrected-data `vis` |
| Uncorrected; calibrate without creating images | `true` | `false` | Run `cal --no-imaging`, then `vis`; analyse archive images with their own documented SDP selections |

```yaml
extra:
  force_calibrate: false  # true only when calibration is required
casa:
  imaging_enabled: false # true when requesting CASA imaging
  datacolumn: data       # corrected SDP values exported into DATA
```

Use `casa.datacolumn: corrected` for imaging calibrated values in `CORRECTED_DATA`. This setting controls imaging, **not `vis`**: visibility QA prefers `CORRECTED_DATA`, then `DATA`, then `MODEL_DATA` and analyses all MS rows. Confirm that the preferred column actually contains calibrated observations. Use field-specific MS exports or the underlying visibility module's `--field` option for multi-field inputs. See [report-refresh visibility guidance](configs/report_refresh/README.md#remote-visibility-preparation-and-commands).

`force_calibrate` applies to all configured reference/test MSs in that invocation. For a mixture of corrected and uncorrected MSs, calibrate only the uncorrected subset with a separate config, then assemble the comparison config. Calibration modifies those MSs and may change flags. Verify calibrator fields and solve settings first. New local calibration does not establish the provenance of pre-existing SDP images; sensitivity must use the visibility selections underlying each actual image.

**Disabling calibration does not disable imaging.** `--no-imaging` overrides the imaging setting for `cal`/`all`; it does not suppress requested calibration or downstream image analysis.

### Recommended launch sequence

Use a dataset config for MS preparation/QA, then a noise-workflow config for the final supplied images:

```bash
tmux new -s mk-analysis
# Edit dataset.yaml from configs/example_local.yaml.
# Only for uncorrected inputs or when requesting new CASA images:
meerkat-ci cal --config dataset.yaml
# After correction; skip the previous command for corrected archive inputs:
meerkat-ci vis --config dataset.yaml
# Only if missing low/high archive images need cuboid extraction:
meerkat-ci low_high_slice --config dataset.yaml

# Edit a separate config from configs/noise_report/l_band.yaml (or s4.yaml).
# Supply the final mfs/low/high image paths and their exact MS selections.
meerkat-ci sensitivity --config noise.yaml
meerkat-ci noise-report --config noise.yaml --reuse-sensitivity
# Alternative to the previous two commands: noise-report generates sensitivity itself.
# meerkat-ci noise-report --config noise.yaml
```

`noise-report` runs **sensitivity → src → xm → pos → flux → report**, omitting sensitivity with `--reuse-sensitivity`. It schedules no calibration, imaging or slicing. `extra.sensitivity.products` supplies its reference/test associations and it derives downstream handoffs; leave `extra.xmatch_pairs`, `extra.positions` and `extra.flux` empty in this config. Include products named `mfs`, `low` and `high` for the detailed three-band flux report. Its automatic handoffs exist only for that invocation, so use a separately wired analysis config for individual report commands later.

For archive products, supply explicit field, data column, scans and native channel selections recovered from observation/imaging metadata. A generic archive metadata file or calibration PDF is not a completed paired-CASA manifest. Use `manifest` only for the repository's completed paired-astrometry manifest for that exact image. Never assume an image used every MS scan/channel from its filename. Missing selections produce unavailable sensitivity with a reason.

Carry visibility evidence into `noise.yaml` or your analysis config:

```yaml
extra:
  verification_report:
    reference_visibility_results: /path/to/interim/reference_ms
    test_visibility_results: /path/to/interim/test_ms
    calibration_status: Not assessed
    calibration_finding: Calibration PDFs supplied; visual review pending.
    image_status: Not assessed
    image_finding: Image diagnostic review pending.
```

Set findings/statuses to the actual review outcome. The consolidated report discovers test-observation calibration PDFs at `<test-image-parent>/../calreports/*.pdf` and diagnostic figures at `<test-image-parent>/../mfimage_frequency_assessment/`. Preserve that layout when transferring archive images or contextual products. Supplying those figures does not itself certify calibration quality. Explicit visibility-result paths are useful when QA and images live on different machines; copy numerical CSV/JSON products and figures, not just PDFs.

---

## How to Run

### 1) Install the project locally (once per machine)

* Install the project **IF** you are turning the repo into an **importable Python package** so wrappers and the CLI work.
* *Tip*: if you use CASA, it runs outside the venv as a separate binary. Just make sure `casa` is on your shell `PATH` (or set `CASA=/path/to/casa`).

---

### 2) Prepare a master config

Make a working copy of the example and fill in your paths (MS files, FITS images, and so on):

```bash
cp configs/example_local.yaml config.yaml
```

Edit `config.yaml` (see contents of `configs/example*_local.yaml config.yaml`):

* `reference.ms_paths` -> corrected MeasurementSets for your reference field
* `tests[].ms_paths` -> corrected MeasurementSets for each test field
* `reference.images` and `tests[].images` -> PB-corrected continuum and/or MFS, low-band, and high-band FITS image files
* `reference.cuboid` and `tests[].cuboid` -> non-PB MFImage cuboids used to
  create missing low/high products; optional `low_image` and `high_image`
  paths reuse existing products
* `frequency_ranges` -> the shared low/high ranges used by CASA and slicing
* `extra.sensitivity` -> per-image MS/field/scan/channel associations and the sensitivity JSON output (required for the new noise workflow)
* `extra.verification_report` -> numerical visibility-result directories and recorded visual findings
* `extra.xmatch_pairs`, `extra.positions`, `extra.flux` -> wire the files your wrappers need to perform cross-matching (cross-matching FITS catalogue pairs), and also cross-matched FITS catalogues for astrometry and flux analysis.
* `paths.*` -> where outputs will be written (defaults live under `data/` and are auto-created)

You can keep multiple configs (for example, one per dataset) and pass `--config path/to.yaml` to each command.

When many paths share directories, define composable values once under
`path_vars` and reference them as `${name}` anywhere else in the config. Path
variables may reference earlier or later path variables. Continue to use YAML
anchors for complete paths that are reused unchanged:

```yaml
path_vars:
  commissioning: "/data/GPU_Correlator_Commissioning"
  observations: "${commissioning}/Obsevations"
  ref_images: "${observations}/reference/images"
  ref_stem: "1785941476_continuum_image_J2147-8132_IClean"

reference:
  images:
    - &ref_image "${ref_images}/${ref_stem}_PB.fits"
  cuboid: "${ref_images}/${ref_stem}.fits"

extra:
  positions:
    - {ref_fits: *ref_image}
```

**`config.yaml` scaffold (excerpt):**

```yaml
project_name: "MeerKAT-correlator-upgrade-imaging-tests"
paths:
  raw_dir: "data/raw"
  interim_dir: "data/interim"
  processed_dir: "data/processed"
  reports_dir: "data/reports"
  sky_xmatches_dir: "data/processed/Sky-CrossMatches"

reference:
  name: "ref_field"
  ms_paths: []
  images: []
  cuboid: "/path/ref_IClean.fits"
# To reuse existing products instead, add low_image and high_image here.

tests:
  - name: "test_field_A"
    ms_paths: []
    images: []
    cuboid: "/path/test_IClean.fits"

frequency_ranges:
  lowband_hz: [8.98e+8, 1.00e+9]
  highband_hz: [1.46e+9, 1.70e+9]

low_high_slice:
  enabled: true
  overwrite: false
  add_to_source_finding: true

pybdsf:
  overwrite: false

extra:
  xmatch_pairs: []
  positions: []
  flux: {}
```

---

### 3) Run individual steps (surgical control)

Each CLI command wraps a helper in `src/meerkat_corr_imaging/steps/...`.
Python-backed helpers use the active interpreter with an installable
`meerkat_corr_imaging` module. CASA-backed helpers pass CASA an absolute path
to the packaged CASA script. No pipeline step resolves its implementation
relative to the launch directory.

#### 3.1 Visibility QA (Step 2)

```bash
python -m meerkat_corr_imaging.cli vis --config config.yaml
```

What happens:

* For each `reference.ms_paths` and each test `ms_paths`, it runs
  `python -m meerkat_corr_imaging.vis_amp_analyze ...`. Test comparisons reuse
  the first reference result, so the reference MeasurementSet is not analysed
  a second time.
* It derives separate auto- and cross-correlation RFI-free channel masks from
  `FLAG` and `FLAG_ROW`. A channel is retained only when its aggregate flagged
  fraction is strictly below 20%.
* It averages unflagged amplitudes by scan, baseline, spectral window and
  polarisation before calculating the mean and spectral RMS. The RMS is
  measured after subtracting a 51-channel, third-order Savitzky--Golay trend.
* It reports fractional amplitude oscillation as robust detrended RMS divided by
  median amplitude. The GPU result is compared with the paired CMC1 result at
  95% confidence. A `Concern` requires the lower GPU-minus-CMC1 interval bound
  to exceed zero. The 20% flagging criterion selects eligible spectral channels.
* Outputs are inspectable CSVs, diagnostic plots, and a concise draft DOCX
  report under the deterministic `data/interim/<msbase>/` directory. A rerun
  replaces same-named generated products in that directory.

Check after running:

* `data/interim/*/acceptance_summary.csv`
* `data/interim/*/rfi_free_channel_mask.csv`
* `data/interim/*/flagging_vs_time.csv`, `flagging_by_channel.csv`,
  `flagging_by_antenna.csv`, and `flagging_by_baseline.csv`
* `data/interim/*/scan_averaged_amp_stats.csv`
* `data/interim/*/scan_averaged_amp_by_antenna.csv` and
  `scan_averaged_amp_by_baseline.csv`
* `data/interim/*/mean_vs_scan_split.png`, `rms_vs_scan_split.png`,
  `oscillation_vs_scan_split.png`, and the remaining QA plots
* `data/interim/*/draft_vis_amp_summary.docx`

`perrow_amp_stats.csv` is still written for compatibility. The comparative
`acceptance_summary.csv` uses whole-scan resampling for aggregate flagging and
physical-baseline resampling for baseline-aggregated p95 oscillation. A lone
observation remains `Not assessed` until paired with CMC1.

#### 3.2 Calibrate & Image with CASA (Step 3)

```bash
python -m meerkat_corr_imaging.cli cal --config config.yaml
```

What happens:

* If `extra.force_calibrate: true`, it runs `standalone_xxyy_solve.py` for each
  distinct configured `reference.ms_paths` / `tests[].ms_paths` input before imaging it.
  Paths are resolved relative to the launch directory; `~` is expanded. Missing
  MS directories fail before CASA starts. Calibration modifies the configured MS.
* Standalone calibration requires `MCI_CAL_MSFILE=/path/to/input.ms` in the
  environment. The packaged script no longer falls back to an example dataset.
  Set `casa.flux_field` to the flux calibrator name or numeric field ID (default
  `J0408-6545`). `casa.refant` defaults to `auto`; set it to an antenna name to
  override automatic selection. Other calibrator fields and solve settings remain
  in the script's USER INPUTS section and must match the observation.
* Automatic reference selection uses the lowest flagged fraction over all
  flux-calibrator cross-correlations, counting both baseline ends and `FLAG_ROW`.
  Ties choose the lowest antenna ID. Missing fields or entirely flagged data fail
  calibration. `refant_stats.json` in each calibration output directory records
  the selected antenna and per-antenna statistics. Antennas above 80% flagged are
  reported, not automatically flagged. This adapts `legacy_scripts/calc_refant.py`
  without its external config-parser, bookkeeping or SLURM dependencies.
  Other CASA code can reuse `meerkat_corr_imaging.calc_refant.get_ref_ant`.
* Then it runs `casa -c <absolute-packaged-path>/tclean_two_bands.py [--scans=...] <all MS>`.
* Products go next to each MS, usually in `<msdir>/images/...` with FITS exported; QA text files are created.

Before running:

* Ensure `casa` is callable: `which casa`. If not, `export CASA=/full/path/to/casa` and rerun.

#### 3.3 Low/high cuboid slices (Step 4)

```bash
python -m meerkat_corr_imaging.cli low_high_slice --config config.yaml
```

The pipeline invokes the packaged extractor as
`python -m meerkat_corr_imaging.low_high_slice`, so this step does not depend
on the current working directory.

What happens:

* Reads the shared `frequency_ranges.lowband_hz` and `highband_hz` values.
  A step-specific value under `low_high_slice` is only needed when slicing
  intentionally uses a different range from CASA.
* For each configured reference/test cuboid, reads `NTERM`, `NSPEC`, and the
  `FRELnnnn`/`FEFFnnnn`/`FREHnnnn` plane metadata.
* Selects the single plane with the largest overlap with each requested range.
  Ties go to the lower-frequency plane for low band and the higher-frequency
  plane for high band.
* Writes non-PB 2-D `_lowband.fits` and `_highband.fits` products beside the
  cuboid. Configured `low_image`/`high_image` paths override these names; a
  missing configured output must still have the same parent directory as its
  cuboid. Existing images elsewhere may be validated and reused.
* Existing configured or deterministic files are validated and reused unless
  `overwrite: true`.

Run this step before `src`. With `add_to_source_finding: true`, the resolved
low/high files are automatically added to the PyBDSF image list.

#### Sensitivity metadata (before flux and report)

```bash
meerkat-ci sensitivity --config noise.yaml
```

This standalone step writes `extra.sensitivity.output_json` from the configured per-product reference/test MS selections. It requires python-casacore, but neither PyBDSF nor existing FITS images are required to calculate it. Prefer generating it after calibration/flagging and after selecting the final images so their identities are recorded. Recalculate after visibility, flag, selection or image changes.

The JSON records theoretical natural-weight Stokes-I RMS, usable parallel-hand exposure, antenna participation, bandwidth/time, band/SEFD assumptions, provenance and unavailable-result reasons. Flux and consolidated reports consume that artifact without reopening MSs. Theoretical sensitivity and measured annular image RMS are different quantities; imaging weights, taper, confusion and PB corrections are not silently applied. See [the sensitivity contract](#sensitivity-artifact-and-reuse-contract) for exact selection and reuse requirements.

For individually launched `flux`/`report`, configure explicit normal handoffs together with the same `extra.sensitivity.products` and JSON path. The older `extra.flux.thermal_noise` route remains available, but must not conflict with configured sensitivity JSON.

#### 3.4 Source finding with PyBDSF (Step 5)

```bash
python -m meerkat_corr_imaging.cli src --config config.yaml
```

What happens:

* Collects images from `reference.images`, each test `images`, resolved
  low/high slice products, and any `extra.images_globs`.
* Runs `python -m meerkat_corr_imaging.pybdsf_srcfind --images ... [--isl ... --pix ... --freq-* ...]`.
* Reuses inputs whose FITS PyBDSF catalogue already exists; the ASCII export is optional. Set
  `pybdsf.overwrite: true` to run source finding again and replace them.
* if input images do not have frequency information in their headers, run this step for each set of images that have the same reference frequency and specify that frequency via `freq_mhz` under the `pybdsf` config section.
* PyBDSF catalogues land near the images or wherever your script writes them (often under `data/processed/...`).

Sanity check:

* Look for generated catalogues (FITS tables), usually `*_gaul.fits` or `*_srl.fits`.

#### 3.5 Cross-matching catalogues (Step 6)

```bash
python -m meerkat_corr_imaging.cli xm --config config.yaml
```

What happens:

* Reads `extra.xmatch_pairs`: each item is `[input1, input2]` or `[input1, input2, output]`.
* Calls `python -m meerkat_corr_imaging.xmatch_pybdsf` for each pair with
  `--max-error` and related options.
* Without an explicit `output`, the wrapper creates one under `data/processed/Sky-CrossMatches/` with a sensible name.

* Optional `extra.survey_xmatch_jobs` and `extra.survey_xmatches` add independent
  local SUMSS/RACS associations with per-input profiles, table selection,
  radii, ambiguity flags, provenance and coverage diagnostics. They are retained
  through image-pipeline wiring and excluded from automatic MeerKAT flux/position
  analyses. See [local survey matching](docs/local_survey_crossmatching.md) and the
  [opt-in MFS configuration](configs/survey_xmatch/local_mfs.yaml).

Verify:

* `data/processed/Sky-CrossMatches/*.fits` exists and has match columns.

#### 3.6 Astrometry (positions) analysis (Step 7a)

```bash
python -m meerkat_corr_imaging.cli pos --config config.yaml
```

What happens:

* Loops over `extra.positions` entries, each with `xmatch_table`, `ref_fits`, `other_fits`, and optional `per_scan_glob`, `otherdatatag`.
* Calls `python -m meerkat_corr_imaging.positions_analysis ...` to generate
  the figures and DOCX interpretation.

Outputs:

* Plots and `draft_*.docx` files under `data/reports/crossmatched-positions/` (or wherever your script writes them).

#### 3.7 Flux analysis (Step 7b)

```bash
python -m meerkat_corr_imaging.cli flux --config config.yaml
```

What happens:

* Reads `extra.flux` with `ref_low_xmatch`, `ref_high_xmatch`, `ref_mfs_xmatch` (required), and optional `scans_glob`, `docx_name`.
* Calls `python -m meerkat_corr_imaging.flux_analysis ...` to produce plots
  plus a DOCX summary.

Outputs:

* Plots and `draft_*.docx` files under `data/reports/flux/` (or your configured location).

##### Theoretical thermal point-source noise

The `flux` report now includes natural-weighting Stokes-I thermal RMS in µJy/beam:

```text
sigma_I = SEFD / sqrt(2 * effective_bandwidth_Hz * on_source_seconds * N * (N-1))
```

Per-antenna mean SEFDs follow **ESDKB-Sensitivity calculators-051026-112758.pdf**,
pages 5–6: L 425 Jy, UHF 550 Jy, S0 365 Jy, S1 364 Jy, S2 365 Jy,
S3 366 Jy and S4 369 Jy. S1 is inferred from overlapping sub-bands in that
reference. A band-mean SEFD is an approximation to its frequency-dependent curve.
The calculator's default antenna counts and RFI losses are **not** imposed on data.

For a single flux comparison with one reference and one test, the step uses each
target's single `ms_paths` entry, `casa.field`, `casa.scans`, `casa.datacolumn`, and
the shared `frequency_ranges` for low/high products. MFS requests all channels.
This inference assumes those are the selections used for the catalogue images;
extra imaging cuts require explicit inputs. Multiple comparisons or multiple MSs
require explicit per-product associations. S0–S4 is taken from an explicit target
name label or `band`; overlapping S sub-bands cannot be identified from a sliced
frequency range alone. L/UHF is inferred only when the full MS frequency coverage
identifies one band unambiguously.

Override selections using `extra.flux.thermal_noise.low/high/mfs.reference/test`.
See `configs/flux_thermal_noise.example.yaml`. Each entry accepts:

* `manifest` and `band`: a **complete paired-astrometry image manifest**. The step
  reads its MS, field, scans, requested channels and data column, then re-reads
  the current MS. Keep the imaged MS flags/data unchanged when using this path.
* `ms`, `field`, optional `scans` (integer list or comma-separated IDs), `band`,
  `datacolumn` (default `data`), `frequency_range_hz` (inclusive channel centres),
  or `channels` (mapping SPW IDs to zero-based channel indices). No implicit
  union of separately imaged scans or selection based on source catalogue rows.
* `visibility_results`, `band`, and **`allow_approximate: true`**: an existing
  **field-specific** visibility QA directory. Automatic discovery from
  `extra.verification_report.reference_visibility_results/test_visibility_results`
  also requires `extra.flux.allow_approximate: true`. Without opt-in the result
  remains unavailable.
  This is explicitly approximate: widths come from native channel spacing,
  dump duration from median within-scan cadence, and channel flag fractions from
  `rfi_free_channel_mask.csv`. The spectral diagnostic's `RFI_FREE` gate is not
  applied as an imaging mask. Weight validity and EXPOSURE metadata were not
  retained. Scan/channel-index restrictions cannot be recovered this way.
* When neither MS nor summaries exist: explicit `antenna_count`,
  `effective_bandwidth_hz`, `on_source_integration_s`, `band`, and a nonempty
  `provenance` description; optionally `approximation`. Never substitute nominal
  bandwidth for unflagged bandwidth. Example metadata entry:
  `{antenna_count: 58, effective_bandwidth_hz: 385000000,
  on_source_integration_s: 3600, band: L, provenance: 'documented planning example'}`.

For MS inputs, N counts antennas with at least one usable selected cross sample;
fully flagged antennas and auto-correlations do not contribute. Time is the union
of selected on-source dump intervals, so baselines/SPWs do not multiply time and
scan gaps/calibrator time are excluded. Usable cells have FLAG/FLAG_ROW clear,
finite selected DATA, and positive finite WEIGHT_SPECTRUM (otherwise WEIGHT).
Both parallel hands must be stored; each valid hand contributes independently.
Positive weights gate eligibility; their magnitudes do not impose image weights.
With `C = sum(valid_hand_cell * abs(CHAN_WIDTH) * EXPOSURE)`, effective bandwidth
is `C / (N*(N-1)*on_source_seconds)`. This equals actual unflagged bandwidth for
a complete constant array, and accounts for partial flagging, missing baselines,
and changing participation. Overlapping selected spectral channels and duplicate
baseline/SPW/time records are rejected to prevent double-counting sensitivity.
The single bandwidth is an exposure-equivalent summary, not merely a frequency
span or a union of intermittently available channels.

Results appear in the normal flux DOCX/PDF section, the existing report JSON at
`thermal_noise.<product>.<reference|test>.theoretical_rms_ujy_beam`, and a new
`<report-stem>.thermal_noise.csv` beside that JSON. These outputs include N,
bandwidth, time, band, SEFD, provenance/assumptions and, where recovered, flag
fractions and channel selection. Missing or unreliable metadata is **unavailable**
with a reason, while existing catalogue flux results still run. Direct CLI calls
accept the same mapping through `--thermal-noise-json`.

This is a theoretical natural-weighting point-source expectation, not measured
image RMS or a replacement for the CMC1-scaled RMS comparison. No Briggs/robust,
tapering, confusion, calibration, dynamic-range or primary-beam corrections are
applied, and no Pass/Concern decision is added. The PDF's page-2 quick-look table
uses **robust −0.5**; its 1-hour L/S4 values (9.1/7.1 µJy/beam) are therefore not
direct checks of the natural-weighting calculation (about 4.44/3.28 µJy/beam with
the PDF's planning N and bandwidth assumptions).

#### 3.8 Consolidated verification report (Step 8)

```bash
python -m meerkat_corr_imaging.cli report --config config.yaml
```

What happens:

* Reads the MFS, low-band and high-band entries already wired through
  `extra.positions` and the catalogue pairs in `extra.xmatch_pairs`.
* Enforces like-for-like primary-beam correction states.
* Assesses positional offsets/p95, median integrated-flux parity and the
  measured/CMC1-scaled RMS comparison at 95% confidence. Missing decision
  uncertainty remains `Not assessed`.
* Consumes configured sensitivity JSON for product-specific noise comparisons;
  theoretical thermal ratios are diagnostic and add no acceptance threshold.
* Adds the available calibration-report PDF inventory and image-diagnostic
  figures, while marking visibility-only and per-scan checks as unassessed when
  those inputs are unavailable.

Outputs:

* `data/reports/imaging_verification/draft_<observation>_imaging_verification.docx`
* A matching PDF exported at CLI finalization, `_metrics.json` sidecar and
  acceptance-metrics figure. Detailed `pos`/`flux` reports are separate outputs.

Optional observation-specific visual findings can be supplied under
`extra.verification_report` using `calibration_status`, `calibration_finding`,
`image_status`, and `image_finding`.

---

### 4) Run grouped pipeline stages

Choose the grouped command by its actual stage list:

| Command | Stages | Configuration behaviour |
| --- | --- | --- |
| `noise-report` | sensitivity → src → xm → pos → flux → report | Derives handoffs from sensitivity products; `--reuse-sensitivity` omits generation |
| `images` | low_high_slice → src → xm → pos → flux → report | Automatically wires image products and catalogue/analysis handoffs |
| `all` | vis → cal → low_high_slice → src → xm → pos → flux → report | Uses explicit configured downstream handoffs |

**Neither `all` nor `images` automatically runs standalone sensitivity.** Generate JSON separately and configure its handoff when required. `all` runs `vis` before `cal`, so use the staged fresh-dataset sequence for initially uncorrected MSs when the final QA must assess corrected data. If calibration/flagging changes the inputs, regenerate sensitivity afterwards. `all` does not automatically discover new CASA image filenames: update image paths and downstream handoffs to the actual outputs.

To run only the image-domain stages, use `images`:

```bash
python -m meerkat_corr_imaging.cli images --config config.yaml
```

This mode runs exactly:

1. `low_high_slice`
2. `src`
3. `xm`
4. `pos`
5. `flux`
6. `report`

Before the first step, it builds a deterministic product plan. Low/high slice
images are included in source finding; each image's expected PyBDSF FITS
catalogue is included in cross-matching; generated cross-match tables are
included in position analysis; and each test with low, high and MFS tables is
included in flux analysis and the consolidated report. Existing
`extra.xmatch_pairs`, `extra.positions` and `extra.flux` entries are preserved,
and the planner fills in missing reference/test band pairs.

A single `reference.images` or `tests[].images` entry is treated as that
target's MFS image. If a target has multiple arbitrary full-band images, the
planner pairs them by configuration order but does not guess which one is the
MFS input required by flux analysis; configure those downstream handoffs
explicitly in that case.

If you have filled the config for every step:

```bash
python -m meerkat_corr_imaging.cli all --config config.yaml
```

Order:

1. `vis`
2. `cal`
3. `low_high_slice`
4. `src`
5. `xm`
6. `pos`
7. `flux`
8. `report`

Each sub-step logs its commands. Optional unconfigured work may skip; invalid or missing requested inputs can fail and stop downstream stages. Review the audit rather than assuming every stage ran.

Each requested step audits its combined stdout/stderr log and records succeeded,
failed, and not-run inputs. The CLI prints these summaries together once at the
end of the whole run, including failed runs and keyboard interrupts. Under each
`[AUDIT] Inputs:` line, `[AUDIT]  Outputs:` lists absolute output paths marked
`NEW`, `OVERWRITTEN`, or `UNCHANGED`. Logs are included; intermediate files in
the configured interim directory and CASA/PyBDSF work tables are grouped by
their parent folder. DOCX/PDF reports remain individually listed, including PDFs
exported during finalization. Unchanged files are listed only in the requested
step's output locations; unrelated historical reports are omitted. File metadata
before and after each step distinguishes rewrites from unchanged files.

Input summaries and complete command output are retained under
`<reports_dir>/pipeline_audits/<timestamp>_<step>.log`. Steps with independent
inputs attempt all of them before reporting failure, so one bad image, MS, or
catalogue pair does not hide the status of the remaining inputs. The `all`
command still stops before downstream steps when the completed step is failed.

---

### 5) Where things go (default)

* Raw archive downloads: `data/raw/`
* Interim QA and CSVs: `data/interim/<msbase>/...`
* Processed products (images, catalogues, cross-matches): `data/processed/...`
  * Cross-matches specifically: `data/processed/Sky-CrossMatches/`
* Sensitivity JSON: `extra.sensitivity.output_json` (choose an explicit path)
* Reports (DOCX/PDF, metrics JSON, figures): `data/reports/...`
  * Every pipeline-generated DOCX basename starts with `draft_`.
  * Per-step pipeline logs and audits: `data/reports/pipeline_audits/*.log`

You can change these in `config.yaml -> paths.*`. Directories are created automatically.

---

### 6) Quick verification checklist

* After `vis`: CSVs and PNGs under `data/interim/*/`
* After `cal`: CASA images in each MS `images/` directory with FITS and QA exports
* After `low_high_slice`: non-PB `_lowband.fits` and `_highband.fits` beside each cuboid
* After `sensitivity`: JSON with the expected products, selections and per-side status/reasons
* After `src`: PyBDSF catalogues (FITS) near images or in `data/processed/`
* After `xm`: matched FITS tables in `data/processed/Sky-CrossMatches/`
* After `pos` and `flux`: `draft_*.docx` files plus plots under `data/reports/`
* After `report`: draft verification DOCX/PDF, JSON metrics sidecar and acceptance figure under `data/reports/imaging_verification/`; page count varies with available evidence
* After every CLI run: inspect the final per-step input/output audit and `produced_reports.log` for absolute report paths

---

### 7) Common gotchas (and fixes)

* `casa` not found -> export `CASA=/full/path/to/casa` or add it to your `PATH`; retry `mci cal`
* Missing low/high files -> configure a valid non-PB MFImage `cuboid` and run `low_high_slice` before `src`
* if input images do not have frequency information in their headers, run the PYBDSF step for each set of images that have the same reference frequency and specify that frequency via `freq_mhz` under the `pybdsf` config section. 
* Permission errors writing under `data/` -> adjust `paths.*` to point at a writable location
* Wrong FITS paths -> update `config.yaml`; wrappers only forward paths
* No outputs appeared -> read the console; wrappers print the exact script command so you can rerun it by hand
* Long CASA runs -> expected; the wrapper streams CASA logs

---

### 8) For dev purposes: Commit your config and results?

* Commit `configs/example_local.yaml` (sanitised), but avoid committing personal `config.yaml` files with private paths.
* Do not commit bulky products unless you have Git LFS set up; prefer a tiny demo under `examples/tiny-demo/`.

---

### TL;DR sequence

```bash
# Once: install and activate an environment with the required radio dependencies.
python -m pip install -e .
tmux new -s mk-analysis
# Edit dataset.yaml and noise.yaml first; see the fresh-dataset workflow above.
# Conditional: cal for uncorrected inputs or new CASA imaging.
# meerkat-ci cal --config dataset.yaml
meerkat-ci vis --config dataset.yaml
# Conditional: extract missing low/high images before configuring final products.
# meerkat-ci low_high_slice --config dataset.yaml
meerkat-ci sensitivity --config noise.yaml
meerkat-ci noise-report --config noise.yaml --reuse-sensitivity
```

For explicit manual analysis instead of `noise-report`, after image preparation
run `sensitivity → src → xm → pos → flux → report` with an analysis config
containing all required handoffs. Do not run multiple grouped modes as consecutive
steps: each is an alternative launch route.

---

## Refresh reports from existing products

Choose the earliest stage affected by the change. Keep the original image/visibility selections and reference association; use fresh output paths for local survey matches, which refuse overwrites.

| Change / available products | Required rerun |
| --- | --- |
| Only report wording/config; valid matched tables and sensitivity JSON exist | `pos`, `flux`, and/or `report` as needed |
| Matching radius/settings or catalogue content changes | `xm` → `pos` → `flux` → `report` |
| PyBDSF settings or input images change | `src` (enable overwrite when replacing catalogues) → `xm` → `pos` → `flux` → `report`; regenerate sensitivity for changed images/selections |
| Visibilities, flags, calibration or selections change | Reassess image provenance; `vis`, sensitivity regeneration and affected image/catalogue stages before reporting |

For a manually wired analysis config with existing matches and valid sensitivity:

```bash
meerkat-ci pos --config analysis.yaml
meerkat-ci flux --config analysis.yaml
meerkat-ci report --config analysis.yaml
# Only the consolidated report is required? Run just the last command.
```

Keep `extra.positions`, `extra.flux`, image/catalogue pairs, the sensitivity JSON path and identical `extra.sensitivity.products` associations in that config. A standalone `report` does not require rerunning `pos` or `flux` to create their DOCX files. Each CLI invocation exports the DOCX files it writes to PDF.

For the automatic noise-workflow route:

```bash
meerkat-ci noise-report --config noise.yaml --reuse-sensitivity
```

This still runs `src`, `xm`, `pos`, `flux` and `report`; with `pybdsf.overwrite: false` it reuses existing deterministic PyBDSF FITS catalogues. Reuse does not validate their freshness or source-finding settings. Omit `--reuse-sensitivity` to regenerate the JSON first.

Reuse sensitivity only while its image/selection associations and underlying MS remain valid. Consumers check recorded image identities and configured associations but do not reopen MSs to detect later flag/visibility changes. Regenerate after those changes and retain the generating config with the JSON. Historical [report-refresh configurations](configs/report_refresh/README.md) use explicit handoffs; follow their individual-step instructions and add product-specific sensitivity associations when needed, rather than assuming they already contain sensitivity JSON wiring.


## Detailed behaviour and advanced workflows

### Draft reports and terminal sessions

The CLI checks `TMUX` and `STY` before loading the pipeline configuration. Outside
those sessions it offers commands with a deterministic session name (at most ten
characters, derived from the config filename and command), and asks whether to
continue. Only `y` or `yes` continues; Enter or unavailable stdin cancels.

Draft DOCX basenames include the observation/epoch suffix. Visibility reports use
their observation output directory; other reports use the nearest timestamped
observation directory, falling back to the configured test name(s). The consolidated
verification report uses its configured test label. Existing reports are not renamed.

At the end of a CLI run, each DOCX written during that run is rendered to a
same-basename PDF beside it using the existing python-docx and Matplotlib
dependencies. No external converter or additional installation is required.
Paragraphs, tables and embedded figures are retained in document order; PDF
pagination and styling are independent of Word layout. An export failure
retains the DOCX, prints the reason, and makes the CLI exit with status 1.
The final audit lists DOCX and PDF paths under the step that produced them.
A report-only list is also saved to `produced_reports.log` in the configured
reports directory (replaced on each run), whose path appears under
`[AUDIT] Step: run_finalization`. That log includes only reports written by that
invocation.

### Visibility plots and reference reuse

Visibility QA writes one baseline class per figure and omits classes and
polarization series with no finite plotted samples. Auto-baselines (same antenna)
are distinct from cross-hand products XY/YX; neither is invented when absent.
Mean, detrended RMS and flagging by baseline use physical antenna separation
from ANTENNA/POSITION in metres, not projected UV distance. Each figure in the
visibility draft has a caption describing its statistic, grouping and limits.

Reference QA in the configured interim directory is reused after a successful
run writes `reference_qa_complete.json`. Reuse requires unchanged MeasurementSet
file sizes/timestamps (excluding lock files), analyser code, and recorded outputs.
Legacy outputs without a completion record are regenerated once. Deleting the
record forces a fresh reference analysis. Test MeasurementSets continue to run
and compare against the reused reference statistics. Reused reports remain in
their original directory; they are not listed as newly produced reports.

### Calibration application and optional imaging

Calibration now applies solutions to every MS field by default. Configure:

```yaml
casa:
  imaging_enabled: false  # calibration-only; default true
  exclude_fields: []  # exact names or numeric IDs, excluded from both operations
  calibration_exclude_fields: []  # additional applycal exclusions
  imaging_exclude_fields: []  # additional imaging exclusions
  datacolumn: corrected  # use calibrated data when imaging
```

Exclusions affect application of calibration, not the fields used to solve gains.
Imaging retains the `casa.field` selection and subtracts excluded fields. Unknown
exclusions fail explicitly; excluding every imaging field skips imaging, while
excluding every calibration field fails. Existing corrected data in excluded
fields is not erased. Use `cal --no-imaging` to override the YAML and suppress
CASA imaging; the flag also works before the subcommand and with `all`. This
controls step 3 imaging only, not downstream analysis of existing image products.

For XX-only, YY-only, or XX/YY-only MeasurementSets, calibration automatically
sets `parang=False` on gain, bandpass and application tasks. The decision uses
all rows of `POLARIZATION.CORR_TYPE` separately for each MS. Leakage solving
and application are skipped for these inputs because cross-hands are absent.
Other correlation layouts retain the existing calibration behavior.

### Validated CASA batch calibration

Every requested calibration (`extra.force_calibrate: true`) now runs through a
batch bootstrap with stdin disconnected, streamed output and a 30-second
heartbeat. CASA exits explicitly after writing its result, including on failure.
Application checks cannot be disabled: task failures and relevant CASA log errors
stop dependent tasks; calibration tables must contain finite unflagged solutions;
selected fields must contain finite unflagged corrected samples. The known
leap-second-table warning is excluded from fatal task-log checks.

```yaml
casa:
  quality_check: false                 # optional; default off
  quality_min_solution_fraction: 0.95
  quality_max_residual: 0.10
  stage_timeout_seconds: 7200          # limit per calibration task
  timeout_seconds: 21600               # total CASA runtime, not inactivity
  shutdown_timeout_seconds: 30         # exit deadline after result appears
```

“Calibration applied successfully” means these execution/output checks passed,
not that scientific quality is certified. Corrected-data validation samples up to
roughly 256 rows per selected field, including all channels/correlations in each
sample. A column's existence alone never passes validation. Table reports record
usable solutions by antenna/SPW/field; partial coverage is reported.

The `cal` step uses the J0408-6545 epoch-2016 flux model from
`MeerCals/fluxcal/J0408_model.py` and `J0408_flux_model_comparison.ipynb`:
`log10(S/Jy) = -0.9790 + 3.3662*x - 1.1216*x^2 + 0.0861*x^3`,
where `x = log10(frequency/MHz)`. For `casa.flux_field` set to `J0408-6545`,
`0408-6545`, or its numeric MS field ID, `setjy` uses `standard="manual"`,
channel-dependent scaling, and zero Stokes Q/U/V. The cubic is converted
algebraically to CASA `spix` coefficients at the median channel frequency of
the selected SPWs, preserving the supplied spectrum across UHF, L and S bands.
Other flux calibrators use the configured CASA standard.

Enabling `quality_check` additionally checks the usable solution fraction in each
output table and the median fractional complex residual versus MODEL_DATA on
sampled flux-calibrator cross-correlations. It enables scratch models in setjy.
“Calibration quality passed” means only these configured checks passed, not an
independent flux-scale or target-image validation. Tune thresholds for your
observation. Keep the flux calibrator among the applied fields when enabling QA.

Run-specific JSON results live in `paths.reports_dir/calibration_results`.
The CLI requires a successful process exit and matching validated result before
starting imaging or marking the MS successful. A quality failure retains the
successful application status in JSON but fails the requested pipeline operation.
Calibration still requires `force_calibrate`; this change does not recalibrate
inputs when that option is disabled.

Validated application uses CASA `applymode='calflagstrict'`: samples lacking
applicable calibration solutions are flagged. This replaces the former `calonly`
behavior so uncalibrated samples cannot count as usable corrected samples. Review
flagging changes alongside the result. Imaging also runs in a finite batch
process with disconnected stdin, a runtime limit, and explicit exit.

### Measurement errors and CMC1 comparison

The comparison stages preserve PyBDSF formal errors, normalize their units, and
write numerical uncertainties into matched FITS tables and JSON sidecars. Existing
configurations use these defaults:

```yaml
extra:
  uncertainty:
    bootstrap_samples: 5000
    confidence_level: 0.6826894921370859  # central Gaussian +/-1 sigma
    random_seed: 20260911
```

`pos` and `flux` pass these settings to the packaged analysis commands; the same
commands accept `--bootstrap-samples`, `--confidence-level`, and `--random-seed`.
The consolidated report reads the resample count and seed from this mapping and
uses 95% intervals for performance decisions. All counts must be integers;
resample count must be at least two. Use at least 5000 for production reports.

Individual flux ratios propagate both catalogue errors. Spherical astrometry
propagates coordinate covariance, including the RA cos(dec) projection. At small
radial offsets Gaussian coordinate Monte Carlo replaces the singular linear
radial approximation. Rigid fits use positional covariance when at least ten
complete-error matches exist; otherwise they retain an unweighted fit with
sampling errors. Formal parameter covariance, bootstrap covariance, fit method,
inlier mask and residual errors are retained numerically. Residual errors include
the correlation with parameters estimated from those same sources.

Catalogue statistics use paired-source percentile bootstrap intervals. For
consolidated separation/flux aggregates, independent Gaussian measurement-error
Monte Carlo on the fixed matched population is also retained when all required
formal errors exist. The adopted bounds enclose both intervals: they preserve a
measurement-error floor without adding noise twice. These are conservative
sensitivity bounds, not exact combined-coverage confidence intervals. Missing
formal errors leave explicitly labelled sampling-only intervals. Proportions use
Wilson intervals, including all-zero/all-one samples. Fewer than two independent
samples cannot establish a sampling interval. P16/P84 and source standard
deviations describe the population's scatter; they are not errors on its median
or mean.

**Acceptance uses CMC1 as the nominal benchmark at 95% confidence.** Position
p95 is tested against a source-covariance noise-only distribution, with a joint
east/north translation test. The median integrated-flux ratio is assessed for
two-sided parity with one. Expected GPU image RMS is the CMC1 annular RMS scaled
by the square root of the reference/test effective unflagged exposure ratio;
the XX cross-correlation exposure proxies the imaged parallel hands. The RMS
judgment assumes comparable visibility and imaging weights. A one-sided
Concern requires the lower RMS-ratio interval bound to exceed one. Missing
decision uncertainty is Not assessed.

Image RMS is the existing annular `1.4826*MAD` estimator. Its sampling standard
error uses the Gaussian MAD asymptotic variance and effective independent beam
count from FITS BMAJ/BMIN and WCS pixel area, capped at the sampled pixel count.
Missing beam/area information or fewer than two effective samples makes its
uncertainty unavailable. This model assumes stationary, locally Gaussian noise;
PB variation, sidelobes and correlated calibration systematics can dominate its
formal error. The RMS ratio propagates independent image-scale errors.

Detailed flux reports use both-axis ODR plus paired-source bootstrap intervals
for fitted gains/slopes/intercepts. Detailed position reports include sampling
intervals for offset, harmonic and circular summaries. All matched-source fields
remain numerical; report formatting uses `value ± error` only for nearly symmetric
intervals and `value [lower, upper]` otherwise.

Visibility dataset comparison CSV/JSON products include whole-scan cluster
bootstrap intervals for diagnostic median amplitude and spectral RMS, plus
paired performance decisions for flagging and oscillation. Aggregate flagging
fractions are resampled by scan. Oscillation p95 is resampled by physical
baseline after taking scan medians. The confidence level is 95% and a
degradation Concern requires the lower GPU-minus-CMC1 bound to exceed zero.
Per-spectrum measurement errors remain unavailable without spectral covariance;
exact flag census and CASA operational checks retain their existing semantics.
External notebook-generated diagnostic figures are embedded as supplied, without
invented uncertainty annotations.

See [the metric audit and assumptions](docs/uncertainty-audit.md) for the full
metric-to-source map and limitations. Catalogue resampling conditions on matching,
quality cuts and inlier selection; it cannot correct association mistakes,
selection truncation or common calibration systematics. Input catalogues currently
contain no cross-source or cross-dataset covariance.

### Paired CMC1–GPU/CMC2 astrometry imaging (CASA 6)

An opt-in, flag/weight-based frequency selection mode images each gain-calibrator
and target scan plus reference-defined low/centred-middle/high all-scan bands.
It supports separate field MSs, whole-channel paired interval overlap, native
diagnostics, dry-run planning, and per-image provenance with stale-output checks.
Existing configurations retain legacy imaging behaviour. See
[paired astrometry instructions](docs/paired_astrometry.md) and the
[L-band sample](configs/paired_astrometry/l_band.yaml) /
[S4 sample](configs/paired_astrometry/s4.yaml).


### Remote sensitivity and existing-image noise reports

Use `configs/noise_report/l_band.yaml` or `configs/noise_report/s4.yaml` on the
server containing the MSs and FITS images. Edit paths, image PB states and exact
selections first. Outputs can be placed in a fresh directory; visibilities never
need to be downloaded. These commands schedule **no calibration, imaging or
low/high slicing**:

```bash
pip install -e '.[radio]'
# Install PyBDSF in the same Python environment if catalogues need to be generated.
tmux new -s noise
meerkat-ci sensitivity --config configs/noise_report/s4.yaml
meerkat-ci noise-report --config configs/noise_report/s4.yaml --reuse-sensitivity
# Or one launch, including sensitivity:
meerkat-ci noise-report --config configs/noise_report/s4.yaml
```

The sensitivity-only command needs NumPy and python-casacore, plus normal package
CLI dependencies, but no PyBDSF or catalogues. It does not require FITS images to
exist yet. Astropy checks explicit PBCOR headers when images exist. The full
workflow reuses the existing step audit and DOCX/PDF export. Source finding
reuses deterministic `pybdsf.results/<base>/<base>-source-cat.fits` when it exists
and `pybdsf.overwrite: false`; the matching ASCII export is optional. Otherwise
PyBDSF is required and processes the supplied image. When moving images to a
server, also copy their `pybdsf.results` directories, preserving relative paths.
Missing reusable FITS catalogues are reported with the exact expected path. Existing catalogue reuse is
not a validation of the source-finding parameters or catalogue freshness.

`extra.sensitivity.products` is the authoritative per-product reference/test
association. Products `mfs`, `low`, `high` feed the existing three-band flux step;
additional products feed positions and the consolidated report. One CMC1/GPU
pair is supported per launch. Do not also configure `extra.xmatch_pairs`,
`extra.positions`, or `extra.flux` for `noise-report`: the workflow creates those
handoffs. Each product specifies `image`, `band`, and boolean `pb_corrected`.
Prefer `manifest` pointing to the **completed paired-astrometry imaging manifest**
for that exact FITS image. Its MS identity must still match the current MS.
The manifest supplies MS, field, scans, channels and data column; conflicting
explicit selections are rejected. It must use the repository's paired-imaging
contract, including `output_paths.fits`, `ms_identity`, and requested channels.

Without a manifest supply the exact documented selection, for example:

```yaml
reference:
  image: /srv/meerkat/CMC1/images/low.fits
  ms: /srv/meerkat/CMC1/corrected.ms
  field: J2147-8132
  datacolumn: CORRECTED_DATA
  scans: [4, 7]
  channels: {'0': [110, 111, 112, 113]} # zero-based native channels; illustrative
  band: S4
  pb_corrected: false
  imaging_context: {weighting: briggs, robust: -0.5, uvtaper: []}
```

Use explicit `scans: all` / `channels: all` only if that image used all selected
field scans / all MS channels. These are not guessed from FITS headers or names.
The low/high channel lists above are examples, not recommended cuts. Copy the
actual selections from the imaging record. Missing/ambiguous selections remain
unavailable with a reason. Optional `frequency_range_hz` further restricts channel
centres; exact channel indices are preferable. Gaincal and target fields require
separate product associations. No MS concatenation or remote SSH connection is
performed by this command.

For an explicitly opted-in cache approximation add `visibility_results` and
`allow_approximate: true` to the association, retaining an explicit MS path,
field, band, PB state, `scans: all`, and `channels: all`. The MS is not opened in
this mode. Only full-field/all-scan cached summaries are admissible; frequency
cuts may be approximated but scan/channel-index cuts cannot. Cache weights,
finite-data validity and true EXPOSURE metadata were not retained. A paired
comparison cannot mix exact MS exposure with cached approximate exposure.

#### Sensitivity artifact and reuse contract

`extra.sensitivity.output_json` (or `extra.sensitivity_json` downstream) names the
atomic JSON output. Different configured JSON paths, or JSON plus explicit
`extra.flux.thermal_noise`, are rejected. With precomputed JSON the flux/report
steps never open an MS or recalculate sensitivity. For individually launched
`flux`/`report`, supply the same `extra.sensitivity.products` associations and
explicit normal positions/flux handoffs, including the image pair for each
crossmatch. The module-level flux CLI accepts `--sensitivity-json` and
`--sensitivity-products-json` for the equivalent handoff.

Schema version 1 / calculation `natural-stokes-i-exposure-1` stores:

* `products.<product>.<reference|test>.theoretical_rms_ujy_beam`, N/antenna IDs,
  effective Hz, on-source seconds, usable parallel-hand Hz s exposure;
* band, per-antenna SEFD, reference, formula, units and assumptions;
* requested association and its digest, resolved selection/actual channel indices
  and scan IDs, MS file-stat identity, image identities, PB state and imaging context;
* `status`, `approximation_status`, and failure `reason`.

MS identities exclude mutable lock files. These are file-stat provenance
fingerprints, not data-content checksums. Consumers validate schema, configured
image/selection associations, and recorded FITS identities. **They deliberately
cannot detect later MS changes without reopening it**: regenerate JSON when
visibilities/flags change. Keep the MS unchanged during generation. A completed
manifest verifies its stored MS identity. Explicit selections assert the user
has identified the visibilities underlying the image. JSON generated before an
image exists has no image-stat snapshot; downstream still checks paths and PB
headers. The JSON must remain with its selection config for reproducibility.

With C the summed eligible baseline/channel/parallel-hand exposure, each product
uses `sigma_expected,GPU = measured_CMC1 * sqrt(C_CMC1/C_GPU)` under equal
band/SEFD and comparable imaging-weight assumptions. The new workflow uses both
parallel hands consistently and never substitutes the older XX-only exposure
proxy. Effective bandwidth is C/[N(N−1)t]; the underlying radiometer formula and
SEFD values above remain unchanged. Positive weights are eligibility gates;
their magnitudes are not a measured imaging sensitivity correction.

The consolidated DOCX/PDF adds all six requested noise metrics, with values in
µJy/beam, product-specific exposure and input audit details. They are retained in
`*_metrics.json` under each band's `noise_comparison`. Flux DOCX/PDF, its normal
metrics JSON and `.thermal_noise.csv` retain the theoretical RMS and provenance.
Unavailable JSON results never fall back to an unrelated exposure convention.
Thermal ratios are diagnostics and do not add acceptance thresholds; the existing
measured/CMC1-scaled RMS decision retains its conditional imaging assumptions.

Measured image RMS remains 1.4826 × MAD in the phase-centred 0.25–0.50 degree
annulus. PB-corrected annular noise versus uncorrected on-axis thermal sensitivity
is explicitly **qualified**, without inventing a PB correction. Supply optional
`non_pb_image` on **both** sides for corresponding uncorrected images from the
same visibility selections; an additional column/JSON comparison is then
measured using the same annulus. Non-PB annular noise can still include imaging
weighting, taper, confusion and calibration effects. Known weighting/taper is
recorded as context, never silently corrected.

---

## Reproducibility & Provenance

* **Config-driven I/O**: all paths and parameters live in YAML configs under `configs/`
* **Version capture**: store `metadata.json` per run with:

  * Python & package versions, CASA version
  * Git commit hash
  * Effective configuration
* **Deterministic figures**: set random seeds where relevant; close Matplotlib figures (`plt.close('all')`) between steps
* **Large files**: track with Git LFS or reference them via `examples/tiny-demo` for quick tests

---

## Examples & Tests

* **Configs:** start with `configs/example_local.yaml`; use `configs/noise_report/` for sensitivity-driven analysis and `configs/paired_astrometry/` for opt-in CASA paired imaging.
* **Tests** (`pytest`):

  * `test_config_path_vars.py` — composable configuration paths
  * `test_image_pipeline.py` — deterministic image-stage handoffs
  * `test_sensitivity_workflow.py` — sensitivity generation, provenance and noise handoffs
  * `test_output_audit.py` / `test_report_runtime.py` — output inventory and report finalization
  * `test_survey_xmatch.py` — local catalogue matching

Run:

```bash
pytest -q
```

---

## Legacy Scripts

Historical code is preserved under **`legacy_scripts/`**. This keeps provenance while avoiding confusion with the current workflow.

---

## Development Practices

* **Pre-commit hooks**:

  ```bash
  pre-commit run --all-files
  ```

  Includes `black`, `ruff`, `isort`, trailing whitespace and EOF fixers.

* **Type hints & logging**:

  * Add type annotations for new/edited functions
  * Use structured logging: `%(asctime)s %(levelname)s %(name)s: %(message)s`
  * Prefer pure, testable functions in `src/…/steps/`; keep `scripts/` as thin entrypoints

* **Subprocess hygiene**:
  Use `subprocess.run([...])` (list form) to call CASA or other tools, not `os.system`.

* **Data handling**:
  Use `pathlib.Path`; avoid hard-coded paths in code; everything should flow from `config.yaml`.

---

## Contributing

1. Create a feature branch:

   ```bash
   git switch -c repo-reorg-YYYY-MM-DD
   ```
2. Make changes with tests where appropriate.
3. Run `pre-commit` and `pytest`.
4. Push and open a Pull Request for review.

For large artifacts, prefer small demo files and reproducible steps over committing entire datasets.

---

## Citation

If this work contributes to published research, please cite the repository.

```
@misc{MeerKAT-Correlator-Imaging-Tests,
  author       = {Legodi, L. S and collaborators},
  title        = {MeerKAT Correlator Upgrade — Imaging Tests},
  year         = {2025},
  howpublished = {\url{https://github.com/Sam-Legodi/MeerKAT-correlator-upgrade-imaging-tests}}
}
```

---
