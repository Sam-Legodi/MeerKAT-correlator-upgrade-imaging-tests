from __future__ import annotations
import subprocess, shlex
from pathlib import Path
from typing import Iterable, Dict, Any
from ..config import Config

def _run(cmd: Iterable[str]):
    print("[POS] ->", " ".join(shlex.quote(c) for c in cmd))
    subprocess.run(list(cmd), check=True)

def run(cfg: Config):
    """
    Calls scripts/positions_analysis.py for each analysis entry specified in config.extra.positions.
    Each entry must supply:
      xmatch_table, ref_fits, other_fits
    Optional: per_scan_glob, otherdatatag
    """
    script = Path("scripts/positions_analysis.py")
    if not script.exists():
        raise FileNotFoundError(f"Missing {script}")

    analyses = cfg.extra.get("positions", [])
    if not analyses:
        print("[POS] No positions analyses defined; skipping.")
        return

    for a in analyses:
        cmd = ["python", str(script),
               "--xmatch-table", a["xmatch_table"],
               "--ref-fits", a["ref_fits"],
               "--other-fits", a["other_fits"]]
        if a.get("per_scan_glob"):
            cmd += ["--per-scan-glob", a["per_scan_glob"]]
        if a.get("otherdatatag"):
            cmd += ["--otherdatatag", a["otherdatatag"]]
        _run(cmd)
