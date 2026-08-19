from __future__ import annotations
import sys
from typing import Iterable
from ..audit import record_skip, register_inputs, run_logged_command
from ..config import Config

def _run(cmd: Iterable[str], *, inputs: Iterable[str] = ()):
    return run_logged_command(cmd, prefix="[FLUX]", inputs=inputs)

def run(cfg: Config):
    """
    Calls the packaged flux analysis using configured crossmatch products.
    Required keys: ref_low_xmatch, ref_high_xmatch, ref_mfs_xmatch
    Optional: scans_glob, docx_name
    """
    fx = cfg.extra.get("flux", {})
    required = ["ref_low_xmatch", "ref_high_xmatch", "ref_mfs_xmatch"]
    if not all(k in fx for k in required):
        message = f"Missing required keys in config.extra.flux; skipping. Need: {required}"
        print(f"[FLUX] {message}")
        record_skip(message)
        return

    inputs = [fx[key] for key in required]
    register_inputs(inputs)

    cmd = [sys.executable, "-m", "meerkat_corr_imaging.flux_analysis",
           "--ref-low-xmatch", fx["ref_low_xmatch"],
           "--ref-high-xmatch", fx["ref_high_xmatch"],
           "--ref-mfs-xmatch", fx["ref_mfs_xmatch"]]
    if fx.get("scans_glob"):
        cmd += ["--scans-glob", fx["scans_glob"]]
    if fx.get("docx_name"):
        cmd += ["--docx-name", fx["docx_name"]]
    _run(cmd, inputs=inputs)
