from __future__ import annotations
import subprocess, shlex
from pathlib import Path
from typing import Iterable, List, Tuple
from ..config import Config

def _run(cmd: Iterable[str]):
    print("[XMATCH] ->", " ".join(shlex.quote(c) for c in cmd))
    subprocess.run(list(cmd), check=True)

def run(cfg: Config):
    """
    Read input pairs from config.extra.xmatch_pairs, each entry:
      - three-item list: [input1, input2, output]
      - or two-item list: [input1, input2] -> auto-make output path under cfg.paths.sky_xmatches_dir
    Then call scripts/xmatch_pybdsf.py per pair.
    """
    script = Path("scripts/xmatch_pybdsf.py")
    if not script.exists():
        raise FileNotFoundError(f"Missing {script}")

    pairs: List[List[str]] = cfg.extra.get("xmatch_pairs", [])
    outroot = Path(cfg.paths.sky_xmatches_dir)
    outroot.mkdir(parents=True, exist_ok=True)

    for entry in pairs:
        if len(entry) < 2:
            raise ValueError(f"xmatch pair needs at least two items: {entry}")
        input1, input2 = entry[0], entry[1]
        if len(entry) >= 3:
            output = entry[2]
        else:
            stem = f"{Path(input1).stem}_X_{Path(input2).stem}.fits"
            output = str(outroot / stem)

        cmd = ["python", str(script), input1, input2, output,
               "--max-error", str(cfg.xmatch.max_sep_arcsec),
               "--ra-col-1", cfg.xmatch.ra_col_1, "--dec-col-1", cfg.xmatch.dec_col_1,
               "--ra-col-2", cfg.xmatch.ra_col_2, "--dec-col-2", cfg.xmatch.dec_col_2,
               "--coord-frame", cfg.xmatch.coord_frame]
        if cfg.xmatch.one_to_many:
            cmd.append("--one-to-many")
        _run(cmd)
