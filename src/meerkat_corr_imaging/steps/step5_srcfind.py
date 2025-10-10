from __future__ import annotations
import subprocess, shlex
from pathlib import Path
from typing import Iterable, List
from ..config import Config

def _run(cmd: Iterable[str]):
    print("[PYBDSF] ->", " ".join(shlex.quote(c) for c in cmd))
    subprocess.run(list(cmd), check=True)

def run(cfg: Config):
    """
    Walk images from config (reference + tests) and call scripts/pybdsf_srcfind.py
    Users may also point this script at directories/globs inside the config via extra.images_globs.
    """
    script = Path("scripts/pybdsf_srcfind.py")
    if not script.exists():
        raise FileNotFoundError(f"Missing {script}")

    images: List[str] = []
    images += cfg.reference.images
    for t in cfg.tests:
        images += t.images
    # optionally allow arbitrary globs in config.extra
    images += cfg.extra.get("images_globs", [])

    if not images:
        print("[PYBDSF] No images given; nothing to do.")
        return

    cmd = ["python", str(script), "--images", *images]
    if cfg.pybdsf.thresh_isl is not None:
        cmd += ["--isl", str(cfg.pybdsf.thresh_isl)]
    if cfg.pybdsf.thresh_pix is not None:
        cmd += ["--pix", str(cfg.pybdsf.thresh_pix)]
    if cfg.pybdsf.freq_hz is not None:
        cmd += ["--freq-hz", str(cfg.pybdsf.freq_hz)]
    elif cfg.pybdsf.freq_mhz is not None:
        cmd += ["--freq-mhz", str(cfg.pybdsf.freq_mhz)]
    if cfg.pybdsf.base_prefix:
        cmd += ["--base-prefix", cfg.pybdsf.base_prefix]

    _run(cmd)
