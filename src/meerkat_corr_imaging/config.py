from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Dict, Any
import yaml

@dataclass
class Target:
    name: str
    ms_paths: List[str] = field(default_factory=list)
    images: List[str] = field(default_factory=list)

@dataclass
class CasaCfg:
    scans: str = ""
    lowband_hz: List[float] = field(default_factory=lambda: [8.98e8, 1.00e9])
    highband_hz: List[float] = field(default_factory=lambda: [1.46e9, 1.70e9])

@dataclass
class PyBDSFCfg:
    thresh_isl: float = 3.0
    thresh_pix: Optional[float] = None
    freq_mhz: Optional[float] = None
    freq_hz: Optional[float] = None
    pblimit: Optional[float] = None
    base_prefix: Optional[str] = None

@dataclass
class XMatchCfg:
    max_sep_arcsec: float = 1.0
    ra_col_1: str = "RA"
    dec_col_1: str = "DEC"
    ra_col_2: str = "RA"
    dec_col_2: str = "DEC"
    coord_frame: str = "icrs"
    one_to_many: bool = False

@dataclass
class PathsCfg:
    raw_dir: str = "data/raw"
    interim_dir: str = "data/interim"
    processed_dir: str = "data/processed"
    reports_dir: str = "data/reports"
    sky_xmatches_dir: str = "data/processed/Sky-CrossMatches"

@dataclass
class Config:
    project_name: str
    out_root: str = "data"
    paths: PathsCfg = field(default_factory=PathsCfg)
    reference: Target = field(default_factory=lambda: Target(name="reference"))
    tests: List[Target] = field(default_factory=list)
    casa: CasaCfg = field(default_factory=CasaCfg)
    pybdsf: PyBDSFCfg = field(default_factory=PyBDSFCfg)
    xmatch: XMatchCfg = field(default_factory=XMatchCfg)
    extra: Dict[str, Any] = field(default_factory=dict)

def _dict_to_dataclass(d: Dict[str, Any]) -> Config:
    paths = PathsCfg(**d.get("paths", {})) if "paths" in d else PathsCfg(
        raw_dir=d.get("raw_dir", "data/raw"),
        interim_dir=d.get("interim_dir", "data/interim"),
        processed_dir=d.get("processed_dir", "data/processed"),
        reports_dir=d.get("reports_dir", "data/reports"),
        sky_xmatches_dir=d.get("sky_xmatches_dir", "data/processed/Sky-CrossMatches"),
    )
    top_level_extra = {k: v for k, v in d.items() if k not in {
        "project_name","out_root","paths","raw_dir","interim_dir","processed_dir","reports_dir",
        "sky_xmatches_dir","reference","tests","casa","pybdsf","xmatch","extra"}}
    merged_extra = {**(d.get("extra") or {}), **top_level_extra}
    cfg = Config(
        project_name=d["project_name"],
        out_root=d.get("out_root", "data"),
        paths=paths,
        reference=Target(**d.get("reference", {"name": "reference"})),
        tests=[Target(**t) for t in d.get("tests", [])],
        casa=CasaCfg(**d.get("casa", {})),
        pybdsf=PyBDSFCfg(**d.get("pybdsf", {})),
        xmatch=XMatchCfg(**d.get("xmatch", {})),
        extra=merged_extra,
    )
    return cfg

def load_config(path: str | Path) -> Config:
    p = Path(path)
    with p.open("r") as f:
        d = yaml.safe_load(f) or {}
    cfg = _dict_to_dataclass(d)
    # ensure dirs exist
    for dpath in [cfg.paths.raw_dir, cfg.paths.interim_dir, cfg.paths.processed_dir, cfg.paths.reports_dir, cfg.paths.sky_xmatches_dir]:
        Path(dpath).mkdir(parents=True, exist_ok=True)
    return cfg
