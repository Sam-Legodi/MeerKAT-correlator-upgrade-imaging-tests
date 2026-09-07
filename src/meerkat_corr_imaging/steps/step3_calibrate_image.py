from __future__ import annotations
import os
import json
from pathlib import Path
from typing import Iterable, List
from ..audit import (
    raise_for_failures,
    record_failure,
    record_skip,
    register_inputs,
    run_logged_command,
)
from ..config import Config

def _run_cmd(
    cmd: List[str],
    env=None,
    *,
    inputs: Iterable[str] = (),
):
    return run_logged_command(cmd, prefix="[CAL/IMAGING]", inputs=inputs, env=env)

def run(cfg: Config):
    """
    Calls your CASA scripts:
      - packaged standalone_xxyy_solve.py (rarely needed; force_calibrate)
      - packaged tclean_two_bands.py for each MS in reference+tests
    """
    casa_bin = os.environ.get("CASA", "casa")  # allow override via CASA env var
    package_dir = Path(__file__).resolve().parents[1]
    tclean_script = package_dir / "tclean_two_bands.py"
    solve_script = package_dir / "standalone_xxyy_solve.py"
    tclean_env = os.environ.copy()
    tclean_env["MCI_CASA_SCRIPT_DIR"] = str(package_dir)
    tclean_env["MCI_TCLEAN_EXCLUDE_FIELDS"] = json.dumps(cfg.casa.exclude_fields + cfg.casa.imaging_exclude_fields)
    tclean_env["MCI_TCLEAN_FIELD"] = cfg.casa.field
    tclean_env["MCI_TCLEAN_DATACOL"] = cfg.casa.datacolumn
    tclean_env["MCI_TCLEAN_IMAGE_SCANS"] = "1" if cfg.casa.image_scans else "0"
    lowband_hz = cfg.casa.lowband_hz or cfg.frequency_ranges.lowband_hz
    highband_hz = cfg.casa.highband_hz or cfg.frequency_ranges.highband_hz
    tclean_env["MCI_TCLEAN_LOWBAND_HZ"] = ",".join(str(v) for v in lowband_hz)
    tclean_env["MCI_TCLEAN_HIGHBAND_HZ"] = ",".join(str(v) for v in highband_hz)

    all_ms = list(cfg.reference.ms_paths)
    for t in cfg.tests:
        all_ms += t.ms_paths
    # Resolve against the launch directory, matching the other pipeline steps.
    all_ms = list(dict.fromkeys(str(Path(ms).expanduser().resolve()) for ms in all_ms))
    register_inputs(all_ms)
    missing = [ms for ms in all_ms if not Path(ms).is_dir()]
    if missing:
        raise FileNotFoundError("Configured MeasurementSet directories not found: " + ", ".join(missing))
    if cfg.extra.get("force_calibrate", False) and not all_ms:
        raise ValueError("force_calibrate requires reference.ms_paths or tests[].ms_paths")

    # imaging for all MS (reference + tests)
    if cfg.casa.imaging_enabled and not tclean_script.is_file():
        raise FileNotFoundError(f"Missing {tclean_script}")

    scans_arg = []
    if cfg.casa.scans.strip():
        scans_arg = [f"--scans={cfg.casa.scans.strip()}"]

    if not all_ms:
        message = "No MS paths in config; skipping imaging."
        print(f"[CAL/IMAGING] {message}")
        record_skip(message)
        return

    failures: list[BaseException] = []
    for ms in all_ms:
        try:
            if cfg.extra.get("force_calibrate", False):
                if not solve_script.is_file():
                    raise FileNotFoundError(f"Missing {solve_script}")
                solve_env = os.environ.copy()
                solve_env["MCI_CASA_SCRIPT_DIR"] = str(package_dir)
                solve_env["MCI_CAL_MSFILE"] = ms
                solve_env["MCI_CAL_EXCLUDE_FIELDS"] = json.dumps(cfg.casa.exclude_fields + cfg.casa.calibration_exclude_fields)
                solve_env["MCI_CAL_REFANT"] = cfg.casa.refant
                solve_env["MCI_CAL_FLUX_FIELD"] = cfg.casa.flux_field
                _run_cmd(
                    [casa_bin, "--nologger", "--log2term", "-c", str(solve_script)],
                    env=solve_env,
                    inputs=[ms],
                )
            if not cfg.casa.imaging_enabled:
                print(f"[CAL/IMAGING] Imaging disabled for {ms}")
                continue
            # Do not let an inherited standalone override replace YAML inputs.
            tclean_env["MCI_TCLEAN_MSFILE"] = ms
            _run_cmd(
                [
                    casa_bin,
                    "--nologger",
                    "--log2term",
                    "-c",
                    str(tclean_script),
                    *scans_arg,
                    ms,
                ],
                env=tclean_env,
                inputs=[ms],
            )
        except Exception as exc:
            record_failure(ms, exc)
            failures.append(exc)

    raise_for_failures("calibration/imaging", failures)
