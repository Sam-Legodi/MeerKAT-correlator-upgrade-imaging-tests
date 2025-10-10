from __future__ import annotations
import subprocess, shlex
from pathlib import Path
from typing import Iterable
from ..config import Config

def _run_cmd(cmd: Iterable[str]):
    print("[VIS-QA] ->", " ".join(shlex.quote(c) for c in cmd))
    subprocess.run(list(cmd), check=True)

def run(cfg: Config):
    """
    For reference + each test target, call scripts/vis_amp_analyze.py
    with MS path(s) from the master config. Outputs are written by the script.
    """
    script = Path("scripts/vis_amp_analyze.py")
    if not script.exists():
        raise FileNotFoundError(f"Missing {script}. Put your real script under scripts/")

    # reference first
    for ms in cfg.reference.ms_paths:
        outdir = Path(cfg.paths.interim_dir) / Path(ms).with_suffix("").name
        cmd = [
            "python", str(script),
            "--ms", ms,
            "--outdir", str(outdir),
        ]
        _run_cmd(cmd)

    # tests (optionally pass reference ms to compare)
    ref_ms0 = cfg.reference.ms_paths[0] if cfg.reference.ms_paths else None
    for t in cfg.tests:
        for ms in t.ms_paths:
            outdir = Path(cfg.paths.interim_dir) / Path(ms).with_suffix("").name
            cmd = ["python", str(script), "--ms", ms, "--outdir", str(outdir)]
            if ref_ms0:
                cmd += ["--ms-ref", ref_ms0]
            _run_cmd(cmd)
