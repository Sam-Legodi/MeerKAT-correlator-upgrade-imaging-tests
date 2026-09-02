from __future__ import annotations
import sys
from pathlib import Path
from typing import Iterable
from ..audit import (
    raise_for_failures,
    record_failure,
    record_skip,
    register_inputs,
    run_logged_command,
)
from ..config import Config

def _run_cmd(cmd: Iterable[str], *, inputs: Iterable[str] = ()):
    return run_logged_command(cmd, prefix="[VIS-QA]", inputs=inputs)

def run(cfg: Config):
    """
    For reference + each test target, call the packaged visibility analyser.
    with MS path(s) from the master config. Outputs are written by the script.
    """
    reference_ms = list(cfg.reference.ms_paths)
    test_ms = [ms for target in cfg.tests for ms in target.ms_paths]
    all_ms = reference_ms + test_ms
    register_inputs(all_ms)
    if not all_ms:
        message = "No MS paths in config; skipping."
        print(f"[VIS-QA] {message}")
        record_skip(message)
        return

    failures: list[BaseException] = []
    reference_results: list[Path] = []

    # reference first
    for ms in reference_ms:
        outdir = Path(cfg.paths.interim_dir) / Path(ms).with_suffix("").name
        cmd = [
            sys.executable, "-m", "meerkat_corr_imaging.vis_amp_analyze",
            "--ms", ms,
            "--outdir", str(outdir),
            "--exact-outdir",
        ]
        try:
            _run_cmd(cmd, inputs=[ms])
            reference_results.append(outdir)
        except Exception as exc:
            record_failure(ms, exc)
            failures.append(exc)

    # Tests reuse the first reference result instead of reading that MS again.
    ref_results0 = reference_results[0] if reference_results else None
    for t in cfg.tests:
        for ms in t.ms_paths:
            outdir = Path(cfg.paths.interim_dir) / Path(ms).with_suffix("").name
            cmd = [
                sys.executable,
                "-m",
                "meerkat_corr_imaging.vis_amp_analyze",
                "--ms",
                ms,
                "--outdir",
                str(outdir),
                "--exact-outdir",
            ]
            if ref_results0:
                cmd += ["--ms-ref-results", str(ref_results0)]
            try:
                _run_cmd(cmd, inputs=[ms])
            except Exception as exc:
                record_failure(ms, exc)
                failures.append(exc)

    raise_for_failures("visibility QA", failures)
