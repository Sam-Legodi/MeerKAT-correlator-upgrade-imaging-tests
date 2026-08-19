from __future__ import annotations
import os
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
    audit_inputs = (["standalone_xxyy_solve"] if cfg.extra.get("force_calibrate", False) else []) + all_ms
    register_inputs(audit_inputs)

    # optional rare calibration step (explicit flag in config)
    if cfg.extra.get("force_calibrate", False):
        if not solve_script.exists():
            raise FileNotFoundError(f"Missing {solve_script}")
        # User maintains inputs inside the CASA script; we just call CASA here.
        _run_cmd(
            [casa_bin, "--nologger", "--log2term", "-c", str(solve_script)],
            inputs=["standalone_xxyy_solve"],
        )

    # imaging for all MS (reference + tests)
    if not tclean_script.is_file():
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
