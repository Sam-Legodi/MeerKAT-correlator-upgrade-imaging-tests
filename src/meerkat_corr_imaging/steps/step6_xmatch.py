from __future__ import annotations
import sys
from pathlib import Path
from typing import Iterable, List
from ..audit import (
    raise_for_failures,
    record_failure,
    record_success,
    record_skip,
    register_inputs,
    run_logged_command,
)
from ..config import Config
from ..survey_xmatch import execute_survey_job, resolve_survey_jobs

def _run(cmd: Iterable[str], *, inputs: Iterable[str] = ()):
    return run_logged_command(cmd, prefix="[XMATCH]", inputs=inputs)

def run(cfg: Config):
    """
    Read input pairs from config.extra.xmatch_pairs, each entry:
      - three-item list: [input1, input2, output]
      - or two-item list: [input1, input2] -> auto-make output path under cfg.paths.sky_xmatches_dir
    Then call the packaged catalogue matcher per pair.
    """
    pairs: List[List[str]] = cfg.extra.get("xmatch_pairs", [])
    has_surveys = cfg.extra.get("survey_xmatch_jobs") or cfg.extra.get("survey_xmatches")
    if not pairs and not has_surveys:
        message = "No cross-match pairs defined; skipping."
        print(f"[XMATCH] {message}")
        record_skip(message)
        return

    outroot = Path(cfg.paths.sky_xmatches_dir)
    outroot.mkdir(parents=True, exist_ok=True)

    tasks: list[tuple[str, str, str, str]] = []
    failures: list[BaseException] = []
    for entry in pairs:
        if len(entry) < 2:
            label = repr(entry)
            error = ValueError(f"xmatch pair needs at least two items: {entry}")
            register_inputs([label])
            record_failure(label, error)
            failures.append(error)
            continue
        input1, input2 = entry[0], entry[1]
        if len(entry) >= 3:
            output = entry[2]
        else:
            stem = f"{Path(input1).stem}_X_{Path(input2).stem}.fits"
            output = str(outroot / stem)
        tasks.append((f"{input1} + {input2}", input1, input2, output))

    register_inputs(label for label, _, _, _ in tasks)

    for label, input1, input2, output in tasks:
        cmd = [sys.executable, "-m", "meerkat_corr_imaging.xmatch_pybdsf",
               input1, input2, output,
               "--max-error", f"{cfg.xmatch.max_sep_arcsec} arcsec",
               "--ra-col-1", cfg.xmatch.ra_col_1, "--dec-col-1", cfg.xmatch.dec_col_1,
               "--ra-col-2", cfg.xmatch.ra_col_2, "--dec-col-2", cfg.xmatch.dec_col_2,
               "--coord-frame", cfg.xmatch.coord_frame]
        if cfg.xmatch.one_to_many:
            cmd.append("--one-to-many")
        try:
            _run(cmd, inputs=[label])
        except Exception as exc:
            record_failure(label, exc)
            failures.append(exc)

    try:
        jobs = resolve_survey_jobs(cfg)
    except Exception as exc:
        label = "survey job configuration"
        register_inputs([label])
        record_failure(label, exc)
        failures.append(exc)
        jobs = []
    for job in jobs:
        label = f"{job['input1']['path']} + {job['input2']['path']} :: {job['output']}"
        register_inputs([label])
        try:
            execute_survey_job(job)
            record_success(label)
        except Exception as exc:
            print(f"FAILED {label}: {exc}")
            record_failure(label, exc)
            failures.append(exc)

    raise_for_failures("cross-matching", failures)
