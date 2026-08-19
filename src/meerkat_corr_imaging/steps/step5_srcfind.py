from __future__ import annotations
import sys
from typing import Iterable, List
from ..audit import record_skip, register_inputs, run_logged_command
from ..config import Config
from .step4_low_high_slice import source_images as low_high_source_images

def _run(cmd: Iterable[str], *, inputs: Iterable[str] = ()):
    return run_logged_command(cmd, prefix="[PYBDSF]", inputs=inputs)

def run(cfg: Config):
    """
    Walk configured images and call the packaged PyBDSF launcher.
    Users may also point this script at directories/globs inside the config via extra.images_globs.
    """
    images: List[str] = []
    images += cfg.reference.images
    for t in cfg.tests:
        images += t.images
    register_inputs(images)
    images += low_high_source_images(cfg)
    # optionally allow arbitrary globs in config.extra
    images += cfg.extra.get("images_globs", [])

    # Preserve configuration order while avoiding duplicate processing.
    images = list(dict.fromkeys(images))
    register_inputs(images)

    if not images:
        message = "No images given; nothing to do."
        print(f"[PYBDSF] {message}")
        record_skip(message)
        return

    cmd = [
        sys.executable,
        "-m",
        "meerkat_corr_imaging.pybdsf_srcfind",
        "--images",
        *images,
    ]
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

    _run(cmd, inputs=images)
