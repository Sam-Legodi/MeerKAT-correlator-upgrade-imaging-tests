# MeerKAT Correlator Upgrade — Imaging Tests

Reproducible workflows to download (manually), calibrate (if needed), image, source-find, cross-match, and analyze MeerKAT **reference** and **test** observations.
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
    * [3.1 Visibility QA (Step 2 / Step 4 repeat)](#31-visibility-qa-step-2--step-4-repeat)
    * [3.2 Calibrate & Image with CASA (Step 3)](#32-calibrate--image-with-casa-step-3)
    * [3.3 Source finding with PyBDSF (Step 5)](#33-source-finding-with-pybdsf-step-5)
    * [3.4 Cross-matching catalogues (Step 6)](#34-cross-matching-catalogues-step-6)
    * [3.5 Astrometry (positions) analysis (Step 7a)](#35-astrometry-positions-analysis-step-7a)
    * [3.6 Flux analysis (Step 7b)](#36-flux-analysis-step-7b)
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
* [License](#license)
* [Citation](#citation)

---

## Overview

**Goal:** Compare test observations against a reference by:

1. Inspecting corrected visibilities,
2. Calibrating and imaging any uncorrected fields in CASA,
3. Running PyBDSF for source catalogues,
4. Cross-matching test vs reference catalogues,
5. Analyzing astrometric offsets and flux consistency, and
6. Producing ready-to-share figures and DOCX reports.

---

## Repository Layout

```
MeerKAT-correlator-upgrade-imaging-tests/
├─ README.md
├─ LICENSE
├─ .gitignore
├─ .gitattributes               # Git LFS for large artifacts (FITS, MS, DOCX)
├─ pyproject.toml               # Tooling config (black/ruff/mypy etc.)
├─ .pre-commit-config.yaml
├─ src/
│  └─ meerkat_corr_imaging/
│     ├─ __init__.py
│     ├─ cli.py                 # `mci` command-line interface (optional)
│     ├─ config.py              # config loader/validator
│     └─ steps/
│        ├─ step1_archive_docs.py     # documentation helpers (no code execution)
│        ├─ step2_vis_analysis.py     # wraps vis_amp_analyze.py
│        ├─ step3_calibrate_image.py  # wraps CASA scripts
│        ├─ step5_srcfind.py          # wraps pybdsf_srcfind.py
│        ├─ step6_xmatch.py           # wraps xmatch_pybdsf.py
│        ├─ step7_positions.py        # wraps positions analysis
│        └─ step7_flux.py             # wraps flux notebook export
├─ scripts/                     # runnable entrypoints (kept stable)
│  ├─ vis_amp_analyze.py
│  ├─ standalone_xxyy_solve.py
│  ├─ tclean_two_bands.py
│  ├─ pybdsf_srcfind.py
│  ├─ xmatch_pybdsf.py
│  ├─ 2025_positions_analysis.py
│  └─ 2025_Flux_analysis.ipynb
├─ configs/
│  ├─ example_local.yaml
│  └─ example_cluster.yaml
├─ data/
│  ├─ raw/                      # archive downloads (manual)
│  ├─ interim/                  # calibration outputs, QA CSVs
│  ├─ processed/                # images, catalogues, xmatches
│  └─ reports/                  # DOCX summaries, PNG/PDF figures
├─ examples/
│  └─ tiny-demo/                # minimal runnable demo (small files)
├─ tests/
│  ├─ test_config.py
│  ├─ test_paths.py
│  └─ test_smoke_pipeline.py
└─ legacy scripts/              # archived historical scripts (read-only)
```

> **Note:** Existing scripts continue to live under `scripts/`. New wrappers in `src/meerkat_corr_imaging/steps/` make the pipeline importable and testable.

---

## Prerequisites

* **Python**: 3.10+
* **System tools**:

  * **CASA** (external binary; accessible on `$PATH`) for calibration/imaging
  * **Git LFS** if you intend to track large artifacts (`*.fits`, `*.ms*`, `*.image`, `*.docx`)
* **Python libraries** (installed during setup): `astropy`, `numpy`, `pandas`, `matplotlib`, `python-docx`, `pybdsf`, etc.

---

## Installation

Create a dedicated virtual environment and install this repository in editable mode so CLI wrappers pick up your edits immediately.

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -e .
```
OR (If a default Conda environment is active. You can still follow the steps—just make sure you’re in the repo folder and pick one environment strategy (either keep using Conda, or use python -m venv; don’t mix them).):

```bash
conda create -n meerkat-ci python=3.10 -y
conda activate meerkat-ci
pip install -e .

```

`pip install -e .` registers the package `meerkat_corr_imaging` (under `src/`), so `python -m meerkat_corr_imaging.cli ...` works anywhere inside the repo. Leave the environment with `deactivate` and reactivate it later with `source .venv/bin/activate`. CASA remains a separate binary; make sure `casa` is on your `PATH` or set `CASA=/path/to/casa`.

Optional but recommended tooling:

```bash
pre-commit install
git lfs install  # large images/catalogues are already listed in .gitattributes
```

---

## Manual Data Download (Step 1)

From the SARAO archive, request **corrected visibilities** by enabling `mvftoms` with `--applyall`. Download the following for both **reference** and **test** observations:

* Corrected visibility **MeasurementSets**
* **PB-corrected** continuum image(s)
* **Multifrequency** image cubes

Document the request IDs, dates, and resulting file paths in `data/raw/README.md`.
This step remains manual due to authentication and archive UX.

---

## How to Run

### 1) Install the project locally (once per machine)

You are turning the repo into an importable Python package so wrappers and the CLI work.

```bash
# create and enter a clean virtualenv
python3 -m venv .venv
source .venv/bin/activate

# install this repo in editable mode (so edits reflect immediately)
pip install -e .
```

What this does:

* `pip install -e .` registers the package `meerkat_corr_imaging` (under `src/`), so `python -m meerkat_corr_imaging.cli ...` works.
* You can leave the venv with `deactivate` and reactivate it later with `source .venv/bin/activate`.

Tip: if you use CASA, it runs outside the venv as a separate binary. Just make sure `casa` is on your shell `PATH` (or set `CASA=/path/to/casa`).

---

### 2) Prepare a master config

Make a working copy of the example and fill in your paths (MS files, FITS images, and so on):

```bash
cp configs/example_local.yaml config.yaml
```

Edit `config.yaml`:

* `reference.ms_paths` -> corrected MeasurementSets for your reference field
* `tests[].ms_paths` -> corrected MeasurementSets for each test field
* `reference.images` and `tests[].images` -> PB-corrected continuum, MFS, low, and high FITS files
* `extra.xmatch_pairs`, `extra.positions`, `extra.flux` -> wire the files your wrappers need
* `paths.*` -> where outputs will be written (defaults live under `data/` and are auto-created)

You can keep multiple configs (for example, one per dataset) and pass `--config path/to.yaml` to each command.

**`config.yaml` scaffold (excerpt):**

```yaml
project_name: "MeerKAT-correlator-upgrade-imaging-tests"
paths:
  raw: "data/raw"
  interim: "data/interim"
  processed: "data/processed"
  reports: "data/reports"

reference:
  name: "ref_field"
  ms_paths: []
  images: []

tests:
  - name: "test_field_A"
    ms_paths: []
    images: []

extra:
  xmatch_pairs: []
  positions: []
  flux: {}
```

---

### 3) Run individual steps (surgical control)

Each CLI command wraps a helper in `src/meerkat_corr_imaging/steps/...`, which in turn runs the actual script under `scripts/` with the arguments pulled from `config.yaml`.

#### 3.1 Visibility QA (Step 2 / Step 4 repeat)

```bash
python -m meerkat_corr_imaging.cli vis --config config.yaml
```

What happens:

* For each `reference.ms_paths` and each test `ms_paths`, it runs `scripts/vis_amp_analyze.py ...`.
* Outputs: CSVs, plots, and a concise DOCX report, typically under `data/interim/<msbase>/` (the wrapper passes `--outdir`).

Check after running:

* `data/interim/*/perrow_amp_stats.csv`
* `data/interim/*/mean_vs_scan_split.png` and other QA plots
* A small `.docx` in `data/reports/` (depending on your script)

#### 3.2 Calibrate & Image with CASA (Step 3)

```bash
python -m meerkat_corr_imaging.cli cal --config config.yaml
```

What happens:

* If `extra.force_calibrate: true`, it runs `standalone_xxyy_solve.py` once (rare).
* Then it runs `casa -c scripts/tclean_two_bands.py [--scans=...] <all MS>`.
* Products go next to each MS, usually in `<msdir>/images/...` with FITS exported; QA text files are created.

Before running:

* Ensure `casa` is callable: `which casa`. If not, `export CASA=/full/path/to/casa` and rerun.

#### 3.3 Source finding with PyBDSF (Step 5)

```bash
python -m meerkat_corr_imaging.cli src --config config.yaml
```

What happens:

* Collects images from `reference.images` plus each test `images` (and any `extra.images_globs`).
* Runs `scripts/pybdsf_srcfind.py --images ... [--isl ... --pix ... --freq-* ...]`.
* PyBDSF catalogues land near the images or wherever your script writes them (often under `data/processed/...`).

Sanity check:

* Look for generated catalogues (FITS tables), usually `*_gaul.fits` or `*_srl.fits`.

#### 3.4 Cross-matching catalogues (Step 6)

```bash
python -m meerkat_corr_imaging.cli xm --config config.yaml
```

What happens:

* Reads `extra.xmatch_pairs`: each item is `[input1, input2]` or `[input1, input2, output]`.
* Calls `scripts/xmatch_pybdsf.py` for each pair with `--max-error` and friends.
* Without an explicit `output`, the wrapper creates one under `data/processed/Sky-CrossMatches/` with a sensible name.

Verify:

* `data/processed/Sky-CrossMatches/*.fits` exists and has match columns.

#### 3.5 Astrometry (positions) analysis (Step 7a)

```bash
python -m meerkat_corr_imaging.cli pos --config config.yaml
```

What happens:

* Loops over `extra.positions` entries, each with `xmatch_table`, `ref_fits`, `other_fits`, and optional `per_scan_glob`, `otherdatatag`.
* Calls `scripts/positions_analysis.py ...` to generate the five figures plus a DOCX with interpretation.

Outputs:

* Plots and DOCX files under `data/reports/crossmatched-positions/` (or wherever your script writes them).

#### 3.6 Flux analysis (Step 7b)

```bash
python -m meerkat_corr_imaging.cli flux --config config.yaml
```

What happens:

* Reads `extra.flux` with `ref_low_xmatch`, `ref_high_xmatch`, `ref_mfs_xmatch` (required), and optional `scans_glob`, `docx_name`.
* Calls `scripts/flux_analysis.py ...` to produce plots plus a DOCX summary.

Outputs:

* Plots and DOCX files under `data/reports/flux/` (or your configured location).

---

### 4) Run the whole pipeline (hands-off)

If you have filled the config for every step:

```bash
python -m meerkat_corr_imaging.cli all --config config.yaml
```

Order:

1. `vis`
2. `cal`
3. `src`
4. `xm`
5. `pos`
6. `flux`

Each sub-step logs the exact command it runs. Missing inputs cause a polite skip with a message.

---

### 5) Where things go (default)

* Raw archive downloads: `data/raw/`
* Interim QA and CSVs: `data/interim/<msbase>/...`
* Processed products (images, catalogues, cross-matches): `data/processed/...`
  * Cross-matches specifically: `data/processed/Sky-CrossMatches/`
* Reports (DOCX, figures): `data/reports/...`

You can change these in `config.yaml -> paths.*`. Directories are created automatically.

---

### 6) Quick verification checklist

* After `vis`: CSVs and PNGs under `data/interim/*/`
* After `cal`: CASA images in each MS `images/` directory with FITS and QA exports
* After `src`: PyBDSF catalogues (FITS) near images or in `data/processed/`
* After `xm`: matched FITS tables in `data/processed/Sky-CrossMatches/`
* After `pos` and `flux`: DOCX files plus plots under `data/reports/`

---

### 7) Common gotchas (and fixes)

* `casa` not found -> export `CASA=/full/path/to/casa` or add it to your `PATH`; retry `mci cal`
* Permission errors writing under `data/` -> adjust `paths.*` to point at a writable location
* Wrong FITS paths -> update `config.yaml`; wrappers only forward paths
* No outputs appeared -> read the console; wrappers print the exact script command so you can rerun it by hand
* Long CASA runs -> expected; the wrapper streams CASA logs

---

### 8) Commit your config and results?

* Commit `configs/example_local.yaml` (sanitised), but avoid committing personal `config.yaml` files with private paths.
* Do not commit bulky products unless you have Git LFS set up; prefer a tiny demo under `examples/tiny-demo/`.

---

### TL;DR sequence

```bash
# (once) setup
python3 -m venv .venv && source .venv/bin/activate
pip install -e .
cp configs/example_local.yaml config.yaml  # edit paths inside

# run a step (or all)
python -m meerkat_corr_imaging.cli vis  --config config.yaml
python -m meerkat_corr_imaging.cli cal  --config config.yaml
python -m meerkat_corr_imaging.cli src  --config config.yaml
python -m meerkat_corr_imaging.cli xm   --config config.yaml
python -m meerkat_corr_imaging.cli pos  --config config.yaml
python -m meerkat_corr_imaging.cli flux --config config.yaml
# or
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

## License

Add your license of choice (e.g., MIT, BSD-3-Clause) in `LICENSE`.

---

## Citation

If this work contributes to published research, please cite the repository. Add a `CITATION.cff` file if you want formal citation metadata.

```
@misc{MeerKAT_corr_upgrade_imaging_tests,
  author       = {Legodi, Sam and collaborators},
  title        = {MeerKAT Correlator Upgrade — Imaging Tests},
  year         = {2025},
  howpublished = {\url{https://github.com/Sam-Legodi/MeerKAT-correlator-upgrade-imaging-tests}}
}
```

---
