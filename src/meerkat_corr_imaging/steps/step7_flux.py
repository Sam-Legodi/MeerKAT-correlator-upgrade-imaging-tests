from __future__ import annotations
import subprocess, shlex
from pathlib import Path
from typing import Iterable
from ..config import Config

def _run(cmd: Iterable[str]):
    print("[FLUX] ->", " ".join(shlex.quote(c) for c in cmd))
    subprocess.run(list(cmd), check=True)

def run(cfg: Config):
    """
    Calls scripts/flux_analysis.py using crossmatch products specified in config.extra.flux.
    Required keys: ref_low_xmatch, ref_high_xmatch, ref_mfs_xmatch
    Optional: scans_glob, docx_name
    """
    script = Path("scripts/flux_analysis.py")
    if not script.exists():
        raise FileNotFoundError(f"Missing {script}")

    fx = cfg.extra.get("flux", {})
    required = ["ref_low_xmatch", "ref_high_xmatch", "ref_mfs_xmatch"]
    if not all(k in fx for k in required):
        print("[FLUX] Missing required keys in config.extra.flux; skipping.\n  Need:", required)
        return

    cmd = ["python", str(script),
           "--ref-low-xmatch", fx["ref_low_xmatch"],
           "--ref-high-xmatch", fx["ref_high_xmatch"],
           "--ref-mfs-xmatch", fx["ref_mfs_xmatch"]]
    if fx.get("scans_glob"):
        cmd += ["--scans-glob", fx["scans_glob"]]
    if fx.get("docx_name"):
        cmd += ["--docx-name", fx["docx_name"]]
    _run(cmd)
