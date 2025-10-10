from __future__ import annotations
import os, subprocess, shlex
from pathlib import Path
from typing import List
from ..config import Config

def _run_cmd(cmd: List[str], env=None):
    print("[CAL/IMAGING] ->", " ".join(shlex.quote(c) for c in cmd))
    subprocess.run(cmd, check=True, env=env)

def run(cfg: Config):
    """
    Calls your CASA scripts:
      - scripts/standalone_xxyy_solve.py (rarely needed; if cfg.extra.get('force_calibrate'))
      - scripts/tclean_two_bands.py for each MS in reference+tests
    """
    casa_bin = os.environ.get("CASA", "casa")  # allow override via CASA env var
    tclean_script = Path("scripts/tclean_two_bands.py")
    solve_script  = Path("scripts/standalone_xxyy_solve.py")

    # optional rare calibration step (explicit flag in config)
    if cfg.extra.get("force_calibrate", False):
        if not solve_script.exists():
            raise FileNotFoundError(f"Missing {solve_script}")
        # user maintains inputs inside the CASA script; we just call CASA here
        _run_cmd([casa_bin, "--nologger", "--log2term", "-c", str(solve_script)])

    # imaging for all MS (reference + tests)
    if not tclean_script.exists():
        raise FileNotFoundError(f"Missing {tclean_script}")

    scans_arg = []
    if cfg.casa.scans.strip():
        scans_arg = [f"--scans={cfg.casa.scans.strip()}"]

    all_ms = list(cfg.reference.ms_paths)
    for t in cfg.tests:
        all_ms += t.ms_paths

    if not all_ms:
        print("[CAL/IMAGING] No MS paths in config; skipping.")
        return

    _run_cmd([casa_bin, "--nologger", "--log2term", "-c", str(tclean_script), *scans_arg, *all_ms])
