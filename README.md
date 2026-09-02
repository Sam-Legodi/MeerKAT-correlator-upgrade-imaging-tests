# MeerKAT Correlator Upgrade — Imaging Tests

Reproducible workflows to download (manually), calibrate (if needed), image, source-find, cross-match, and analyze MeerKAT reference and test observations for correlator upgrade imaging verification. This repository standardizes the end-to-end process and makes it easy to share results with collaborators.
Reproducible workflows to download (manually), calibrate (if needed), image, source-find, cross-match, and analyze MeerKAT **reference** and **test** observations for correlator upgrade imaging verification.
This repository standardizes the end-to-end process and makes it easy to share results with collaborators.

---

## Table of Contents

* [Overview](#overview)
* [Repository Layout](#repository-layout)
* [Prerequisites](#prerequisites)
* [Installation](#installation)
* [Manual Data Download (Step 1)](#manual-data-download-step-1)
* [How to Run](#how-to-run)

  * [1) Install the project locally (once per machine)](#1-install-the-project-locally-once-per-machine)
  * [2) Prepare a master config](#2-prepare-a-master-config)
  * [3) Run individual steps (surgical control)](#3-run-individual-steps-surgical-control)
    * [3.1 Visibility QA (Step 2)](#31-visibility-qa-step-2)
    * [3.2 Calibrate & Image with CASA (Step 3)](#32-calibrate--image-with-casa-step-3)
    * [3.3 Low/high cuboid slices (Step 4)](#33-lowhigh-cuboid-slices-step-4)
    * [3.4 Source finding with PyBDSF (Step 5)](#34-source-finding-with-pybdsf-step-5)
    * [3.5 Cross-matching catalogues (Step 6)](#35-cross-matching-catalogues-step-6)
    * [3.6 Astrometry (positions) analysis (Step 7a)](#36-astrometry-positions-analysis-step-7a)
    * [3.7 Flux analysis (Step 7b)](#37-flux-analysis-step-7b)
    * [3.8 Consolidated verification report (Step 8)](#38-consolidated-verification-report-step-8)
  * [4) Run the whole pipeline (hands-off)](#4-run-the-whole-pipeline-hands-off)
  * [5) Where things go (default)](#5-where-things-go-default)
  * [6) Quick verification checklist](#6-quick-verification-checklist)
  * [7) Common gotchas (and fixes)](#7-common-gotchas-and-fixes)
  * [8) Commit your config and results?](#8-commit-your-config-and-results)
  * [TL;DR sequence](#tldr-sequence)
* [Reproducibility & Provenance](#reproducibility--provenance)
* [Examples & Tests](#examples--tests)
* [Legacy Scripts](#legacy-scripts)
* [Development Practices](#development-practices)
* [Contributing](#contributing)
* [Citation](#citation)

---

## Overview

**Goal:** Compare test observations against a reference by:

1. Inspecting corrected visibilities,
2. Calibrating and imaging any uncorrected fields in CASA,
3. Extracting configured low/high planes from MFImage cuboids when needed,
4. Running PyBDSF for source catalogues,
5. Cross-matching test vs reference catalogues,
6. Analyzing astrometric offsets and flux consistency, and
7. Producing ready-to-share figures, per-analysis DOCX files and a consolidated verification report.

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
├─ examples/
│  └─ tiny-demo/                # minimal runnable demo dataset
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

From the SARAO archive, request **corrected visibilities** by enabling `mvftoms` with `--applyall`. Download the following for both **reference** and **test** observations:

* Corrected visibility **MeasurementSets**: 

e.g: 
```bash
mvftoms.py 1757723806_sdp_l0.full.rdb -f --flags 'static, cam, data_lost, ingest_rfi, predicted_rfi, cal_rfi, postproc' --chanbin 4 --applycal 'all' -o 1757723806_sdp_l0.full.SDP_Allcorrected.ms
```

* **PB-corrected** continuum image(s) -- if comparing PB corrected images.
* **Multifrequency** image cubes (the default SARAO archive/Obit versions)

Document the request IDs, dates, and resulting file paths in `data/raw/README.md`.
This step remains manual due to authentication and archive UX.

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
* It reports fractional amplitude oscillation as detrended RMS divided by mean
  amplitude. Flagging below 20% and oscillation below 1% are `Pass`; values at
  or above either limit are `Concern`.
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

`perrow_amp_stats.csv` is still written for compatibility, but acceptance is
based on the scan-averaged products. The class/polarisation oscillation result
is conservative: any assessed scan/baseline spectrum at or above 1% makes that
class/polarisation a `Concern`.

#### 3.2 Calibrate & Image with CASA (Step 3)

```bash
python -m meerkat_corr_imaging.cli cal --config config.yaml
```

What happens:

* If `extra.force_calibrate: true`, it runs `standalone_xxyy_solve.py` once (rare).
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

#### 3.4 Source finding with PyBDSF (Step 5)

```bash
python -m meerkat_corr_imaging.cli src --config config.yaml
```

What happens:

* Collects images from `reference.images`, each test `images`, resolved
  low/high slice products, and any `extra.images_globs`.
* Runs `python -m meerkat_corr_imaging.pybdsf_srcfind --images ... [--isl ... --pix ... --freq-* ...]`.
* Reuses inputs whose FITS and ASCII PyBDSF catalogues already exist. Set
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

#### 3.8 Consolidated verification report (Step 8)

```bash
python -m meerkat_corr_imaging.cli report --config config.yaml
```

What happens:

* Reads the MFS, low-band and high-band entries already wired through
  `extra.positions` and the catalogue pairs in `extra.xmatch_pairs`.
* Enforces like-for-like primary-beam correction states.
* Scores positional p95, robust median integrated-flux error and the
  measured/CMC1 robust-RMS ratio against the configured verification policy.
* Adds the available calibration-report PDF inventory and image-diagnostic
  figures, while marking visibility-only and per-scan checks as unassessed when
  those inputs are unavailable.

Outputs:

* `data/reports/imaging_verification/draft_<observation>_imaging_verification.docx`
* A matching `_metrics.json` sidecar and acceptance-metrics figure.

Optional observation-specific visual findings can be supplied under
`extra.verification_report` using `calibration_status`, `calibration_finding`,
`image_status`, and `image_finding`.

---

### 4) Run the whole pipeline (hands-off)

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

Each sub-step logs the exact command it runs. Missing inputs cause a polite skip with a message.

Immediately after every requested step, the CLI audits the combined stdout/stderr
log and prints an input-level summary: succeeded, failed, and not-run inputs. The
same summary and complete command output are retained under
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
* Reports (DOCX, figures): `data/reports/...`
  * Every pipeline-generated DOCX basename starts with `draft_`.
  * Per-step pipeline logs and audits: `data/reports/pipeline_audits/*.log`

You can change these in `config.yaml -> paths.*`. Directories are created automatically.

---

### 6) Quick verification checklist

* After `vis`: CSVs and PNGs under `data/interim/*/`
* After `cal`: CASA images in each MS `images/` directory with FITS and QA exports
* After `low_high_slice`: non-PB `_lowband.fits` and `_highband.fits` beside each cuboid
* After `src`: PyBDSF catalogues (FITS) near images or in `data/processed/`
* After `xm`: matched FITS tables in `data/processed/Sky-CrossMatches/`
* After `pos` and `flux`: `draft_*.docx` files plus plots under `data/reports/`
* After `report`: a five-page draft verification DOCX, JSON metrics sidecar and acceptance figure under `data/reports/imaging_verification/`

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
# (once) setup
python3 -m venv .venv && source .venv/bin/activate
pip install -e .
# optionally: pip install -e ".[dev]" or ".[radio]"
cp configs/example_local.yaml config.yaml  # edit paths inside

# run a step (or all)
python -m meerkat_corr_imaging.cli vis  --config config.yaml
python -m meerkat_corr_imaging.cli cal  --config config.yaml
python -m meerkat_corr_imaging.cli src  --config config.yaml
python -m meerkat_corr_imaging.cli xm   --config config.yaml
python -m meerkat_corr_imaging.cli pos  --config config.yaml
python -m meerkat_corr_imaging.cli flux --config config.yaml
# or
python -m meerkat_corr_imaging.cli images --config config.yaml  # image-domain only
python -m meerkat_corr_imaging.cli all  --config config.yaml
```

> **To do:** add a `Makefile` with shortcuts (`make vis`, `make all`) and a tiny `examples/tiny-demo` config for a reproducible miniature run.

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

* **`examples/tiny-demo/`**: a minimal dataset and config to run the full pipeline quickly.
* **Tests** (`pytest`):

  * `test_config.py` — config loading/validation
  * `test_paths.py` — path conventions
  * `test_smoke_pipeline.py` — end-to-end smoke test on the tiny demo

Run:

```bash
pytest -q
```

---

## Legacy Scripts

Historical code is preserved (read-only) under **`legacy scripts/`**. This keeps provenance while avoiding confusion with the current workflow.

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
