from __future__ import annotations
import sys
from typing import Any, Iterable
from ..audit import (
    raise_for_failures,
    record_failure,
    record_skip,
    register_inputs,
    run_logged_command,
)
from ..config import Config

def _run(cmd: Iterable[str], *, inputs: Iterable[str] = ()):
    return run_logged_command(cmd, prefix="[FLUX]", inputs=inputs)

def _configured_analyses(raw: Any) -> list[dict[str, Any]]:
    if not raw:
        return []
    if isinstance(raw, dict):
        return [raw]
    if isinstance(raw, list):
        return raw
    raise ValueError("config.extra.flux must be a mapping or a list of mappings")

def run(cfg: Config):
    """
    Calls the packaged flux analysis using one or more configured crossmatch sets.
    Required keys: ref_low_xmatch, ref_high_xmatch, ref_mfs_xmatch
    Optional: scans_glob, docx_name
    """
    required = ["ref_low_xmatch", "ref_high_xmatch", "ref_mfs_xmatch"]
    analyses = _configured_analyses(cfg.extra.get("flux"))
    if not analyses:
        message = f"No flux analyses defined; skipping. Each entry needs: {required}"
        print(f"[FLUX] {message}")
        record_skip(message)
        return

    labels = [
        " + ".join(str(analysis.get(key, f"<missing {key}>")) for key in required)
        for analysis in analyses
    ]
    register_inputs(labels)
    failures: list[BaseException] = []

    for label, analysis in zip(labels, analyses):
        try:
            missing = [key for key in required if key not in analysis]
            if missing:
                raise ValueError(f"Flux analysis is missing required keys: {missing}")
            cmd = [sys.executable, "-m", "meerkat_corr_imaging.flux_analysis",
                   "--ref-low-xmatch", analysis["ref_low_xmatch"],
                   "--ref-high-xmatch", analysis["ref_high_xmatch"],
                   "--ref-mfs-xmatch", analysis["ref_mfs_xmatch"]]
            if analysis.get("scans_glob"):
                cmd += ["--scans-glob", analysis["scans_glob"]]
            if analysis.get("docx_name"):
                cmd += ["--docx-name", analysis["docx_name"]]
            _run(cmd, inputs=[label])
        except Exception as exc:
            record_failure(label, exc)
            failures.append(exc)

    raise_for_failures("flux analysis", failures)
