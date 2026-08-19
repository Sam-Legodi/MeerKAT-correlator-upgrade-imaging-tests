from __future__ import annotations
import sys
from typing import Iterable
from ..audit import (
    raise_for_failures,
    record_failure,
    record_skip,
    register_inputs,
    run_logged_command,
)
from ..config import Config

def _run(cmd: Iterable[str], *, inputs: Iterable[str] = ()):
    return run_logged_command(cmd, prefix="[POS]", inputs=inputs)

def run(cfg: Config):
    """
    Calls the packaged position analysis for each configured entry.
    Each entry must supply:
      xmatch_table, ref_fits, other_fits
    Optional: per_scan_glob, otherdatatag
    """
    analyses = cfg.extra.get("positions", [])
    if not analyses:
        message = "No positions analyses defined; skipping."
        print(f"[POS] {message}")
        record_skip(message)
        return

    labels = [
        f"{analysis.get('xmatch_table', '<missing xmatch_table>')} + "
        f"{analysis.get('ref_fits', '<missing ref_fits>')} + "
        f"{analysis.get('other_fits', '<missing other_fits>')}"
        for analysis in analyses
    ]
    register_inputs(labels)
    failures: list[BaseException] = []

    for label, analysis in zip(labels, analyses):
        try:
            cmd = [sys.executable, "-m", "meerkat_corr_imaging.positions_analysis",
                   "--xmatch-table", analysis["xmatch_table"],
                   "--ref-fits", analysis["ref_fits"],
                   "--other-fits", analysis["other_fits"]]
            if analysis.get("per_scan_glob"):
                cmd += ["--per-scan-glob", analysis["per_scan_glob"]]
            if analysis.get("otherdatatag"):
                cmd += ["--otherdatatag", analysis["otherdatatag"]]
            _run(cmd, inputs=[label])
        except Exception as exc:
            record_failure(label, exc)
            failures.append(exc)

    raise_for_failures("position analysis", failures)
