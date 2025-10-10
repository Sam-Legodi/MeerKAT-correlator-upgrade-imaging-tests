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
* [Configuration](#configuration)
* [Usage](#usage)

  * [Step 2 — Visibility QA](#step-2--visibility-qa)
  * [Step 3 — Calibrate & Image (CASA)](#step-3--calibrate--image-casa)
  * [Step 4 — QA for Newly Corrected Fields](#step-4--qa-for-newly-corrected-fields)
  * [Step 5 — Source Finding (PyBDSF)](#step-5--source-finding-pybdsf)
  * [Step 6 — Cross-Match Catalogues](#step-6--crossmatch-catalogues)
  * [Step 7 — Astrometry & Flux Analyses](#step-7--astrometry--flux-analyses)
* [Outputs](#outputs)
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

Create a local virtual environment and install this package in editable mode:

```bash
# Option A: using uv (fast)
uv venv .venv && source .venv/bin/activate
uv pip install -e .
pre-commit install

# Option B: standard tools
python3 -m venv .venv
source .venv/bin/activate
pip install -e .
pre-commit install
```

Enable Git LFS (recommended if storing artifacts here):

```bash
git lfs install
# .gitattributes already includes patterns like:
# *.fits, *.ms*, *.image, *.mms, *.docx
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

## Configuration

Copy an example config and edit paths/parameters:

```bash
cp configs/example_local.yaml config.yaml
```

**`configs/example_local.yaml` (excerpt):**

```yaml
project_name: "MeerKAT-correlator-upgrade-imaging-tests"
out_root: "data"

raw_dir: "data/raw"
interim_dir: "data/interim"
processed_dir: "data/processed"
reports_dir: "data/reports"

reference:
  name: "ref_field"
  ms_paths: []     # corrected MS paths (from archive)
  images: []       # PB-corrected continuum/cubes

tests:
  - name: "test_field_A"
    ms_paths: []
    images: []

pybdsf:
  thresh_isl: 3.5
  thresh_pix: 5.0
  rms_box: [50, 10]

xmatch:
  max_sep_arcsec: 1.0

casa:
  scans: ""     # "1,2,5" or empty for all
  lowband_hz: [8.98e8, 1.00e9]
  highband_hz: [1.46e9, 1.70e9]

report:
  include_interpretation: true
```

---

## Usage

You can run scripts directly from `scripts/` **or** via the optional CLI in `src/meerkat_corr_imaging/cli.py` (exposed as `mci` if configured).

### Step 2 — Visibility QA

Analyze corrected visibilities for reference and tests.

**Direct script:**

```bash
python scripts/vis_amp_analyze.py --config config.yaml
```

**Generates:**

* CSVs:
  `perrow_amp_stats.csv` (MEAN, RMS, FLAG_FRAC, ANT1, ANT2, SCAN, DDID, SPW, POL),
  `mean_vs_scan.csv`, `rms_vs_scan.csv`, `flagging_vs_scan.csv`, `flagging_by_channel.csv`
* Figures:
  `mean_vs_scan_split.png`, `rms_vs_scan_split.png`,
  `mean_per_antenna_split.png`, `rms_per_antenna_split.png`,
  `mean_vs_baseline_split.png`, `rms_vs_baseline_split.png`,
  `flagging_vs_scan_split.png`, `flagging_per_antenna_split.png`,
  `flagging_by_channel.png`
* A concise `.docx` report

### Step 3 — Calibrate & Image (CASA)

For any fields not already corrected in the archive:

1. **Calibration** (`standalone_xxyy_solve.py`): delay → bandpass → gaincal (amp & phase) → flux bootstrap → leakage (D-terms) → applycal.
2. **Imaging** (`tclean_two_bands.py`):

   * lowband: 8.98e8–1.00e9 Hz
   * highband: 1.46e9–1.70e9 Hz
   * If bands not covered, fallback to lower/upper **25%** of available span.
   * Produces per-band images and optional per-scan images; FITS and QA text emitted alongside.

**CASA 6 usage (example):**

```bash
casa --nologger --log2term -c scripts/tclean_two_bands.py [--scans=1,2,3] <ms1> <ms2> ...
```

> The `tclean_two_bands.py` docstring in `scripts/` describes run modes, scan selection, and QA-only behavior.

### Step 4 — QA for Newly Corrected Fields

Re-run Step 2 on any newly corrected visibilities produced by Step 3.

### Step 5 — Source Finding (PyBDSF)

Run PyBDSF on:

* Per-scan MFS images
* All-scans full-band MFS
* Low-band (lower quartile) MFS
* High-band (upper quartile) MFS

**Command:**

```bash
python scripts/pybdsf_srcfind.py --config config.yaml
```

### Step 6 — Cross-Match Catalogues

Cross-match each test catalogue with the corresponding reference catalogue.

**Command:**

```bash
python scripts/xmatch_pybdsf.py --config config.yaml
```

Outputs one matched catalogue per pair.

### Step 7 — Astrometry & Flux Analyses

**Positions** (`2025_positions_analysis.py`):

* Reads a cross-match FITS table (reference × other)
* Computes ΔRA/ΔDec (east+, north+) using spherical geometry
* Produces five figures (pos_var CMC1×CMC2)
* Optionally overlays per-scan cross-match tables on the same axes
* Draws enclosing ellipses (center = mean, width/height = 2×max deviation)
* Legends outside axes to avoid obscuring markers
* Exports a DOCX report with physical interpretations

**CLI example:**

```bash
python scripts/2025_positions_analysis.py \
  --xmatch-table /path/to/refXother.fits \
  --ref-fits /path/to/reference_image.fits \
  --other-fits /path/to/other_image.fits \
  --per-scan-glob "/path/to/per_scan/refXscan*.fits"
```

**Flux** (`2025_Flux_analysis.ipynb`):

* Compares matched **L-band** and **UHF-band** PyBDSF catalogues against a reference continuum image
* Quantifies flux consistency, builds diagnostic plots, and compiles a DOCX summary

---

## Outputs

* `data/interim/` — calibration products, QA CSVs, intermediary tables
* `data/processed/` — images, PyBDSF catalogues, cross-matches
* `data/reports/` — DOCX reports and figures from Steps 2, 3 (QA), and 7

Each run should record a `run_id` (e.g., `YYYYMMDD_HHMMSS`) and write a `metadata.json` with tool versions, git commit, and config snapshot.

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

### Quick Commands (optional)

```bash
# Environment
uv venv .venv && source .venv/bin/activate
uv pip install -e .
pre-commit install

# Run selected stages (examples):
python scripts/vis_amp_analyze.py --config config.yaml
casa --nologger --log2term -c scripts/tclean_two_bands.py <ms...>
python scripts/pybdsf_srcfind.py --config config.yaml
python scripts/xmatch_pybdsf.py --config config.yaml
python scripts/2025_positions_analysis.py --xmatch-table ... --ref-fits ... --other-fits ...
# Flux analysis (script version; no notebook needed)
# Option A: via project CLI wrapper (recommended)
python -m meerkat_corr_imaging.cli flux --config config.yaml

# Option B: call the script directly
python scripts/flux_analysis.py --config config.yaml

# Outputs:
# - Figures and tables under data/reports/flux/ (or your configured reports_dir)
# - DOCX summary: data/reports/flux/flux_summary.docx

```

This README is intended to be a single, clear starting point for you and your colleagues to rerun analyses or extend the workflow.
