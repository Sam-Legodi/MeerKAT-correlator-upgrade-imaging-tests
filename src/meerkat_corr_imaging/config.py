from __future__ import annotations
from dataclasses import dataclass, field
import math
from pathlib import Path
import re
from typing import List, Optional, Dict, Any
import yaml

SLICE_PATH_FIELDS = ("cuboid", "low_image", "high_image")
PATH_VAR_PATTERN = re.compile(r"\$\{([A-Za-z_][A-Za-z0-9_]*)\}")

@dataclass
class Target:
    name: str
    ms_paths: List[str] = field(default_factory=list)
    images: List[str] = field(default_factory=list)
    cuboid: Optional[str] = None
    low_image: Optional[str] = None
    high_image: Optional[str] = None

@dataclass
class FrequencyRangesCfg:
    lowband_hz: List[float] = field(default_factory=lambda: [8.98e8, 1.00e9])
    highband_hz: List[float] = field(default_factory=lambda: [1.46e9, 1.70e9])

@dataclass
class CasaCfg:
    quality_check: bool = False
    quality_min_solution_fraction: float = 0.95
    quality_max_residual: float = 0.1
    stage_timeout_seconds: float = 7200
    timeout_seconds: float = 21600
    shutdown_timeout_seconds: float = 30
    imaging_enabled: bool = True
    exclude_fields: List[str] = field(default_factory=list)
    calibration_exclude_fields: List[str] = field(default_factory=list)
    imaging_exclude_fields: List[str] = field(default_factory=list)
    refant: str = "auto"
    flux_field: str = "J0408-6545"
    scans: str = ""
    lowband_hz: Optional[List[float]] = None
    highband_hz: Optional[List[float]] = None
    field: str = "J2147-8132"
    datacolumn: str = "DATA"
    image_scans: bool = True

@dataclass
class LowHighSliceCfg:
    enabled: bool = False
    overwrite: bool = False
    add_to_source_finding: bool = True
    lowband_hz: Optional[List[float]] = None
    highband_hz: Optional[List[float]] = None

@dataclass
class PyBDSFCfg:
    thresh_isl: float = 3.0
    thresh_pix: Optional[float] = None
    freq_mhz: Optional[float] = None
    freq_hz: Optional[float] = None
    pblimit: Optional[float] = None
    base_prefix: Optional[str] = None
    adaptive_rms_box: bool = True
    overwrite: bool = False

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
    frequency_ranges: FrequencyRangesCfg = field(default_factory=FrequencyRangesCfg)
    casa: CasaCfg = field(default_factory=CasaCfg)
    low_high_slice: LowHighSliceCfg = field(default_factory=LowHighSliceCfg)
    pybdsf: PyBDSFCfg = field(default_factory=PyBDSFCfg)
    xmatch: XMatchCfg = field(default_factory=XMatchCfg)
    extra: Dict[str, Any] = field(default_factory=dict)


def _expand_path_vars(raw: Dict[str, Any]) -> Dict[str, Any]:
    """Expand recursive ``${name}`` variables throughout a raw config."""
    path_vars = raw.get("path_vars") or {}
    if not isinstance(path_vars, dict):
        raise ValueError("path_vars must be a mapping of names to strings")

    definitions: Dict[str, str] = {}
    for name, value in path_vars.items():
        if not isinstance(name, str) or not re.fullmatch(
            r"[A-Za-z_][A-Za-z0-9_]*", name
        ):
            raise ValueError(f"Invalid path_vars name: {name!r}")
        if not isinstance(value, str):
            raise ValueError(f"path_vars.{name} must be a string")
        definitions[name] = value

    resolved: Dict[str, str] = {}
    resolving: List[str] = []

    def resolve(name: str) -> str:
        if name in resolved:
            return resolved[name]
        if name not in definitions:
            raise ValueError(f"Undefined path variable: {name}")
        if name in resolving:
            cycle = " -> ".join([*resolving, name])
            raise ValueError(f"Cyclic path_vars reference: {cycle}")

        resolving.append(name)
        try:
            value = PATH_VAR_PATTERN.sub(
                lambda match: resolve(match.group(1)), definitions[name]
            )
        finally:
            resolving.pop()
        resolved[name] = value
        return value

    for name in definitions:
        resolve(name)

    def expand(value: Any) -> Any:
        if isinstance(value, str):
            return PATH_VAR_PATTERN.sub(lambda match: resolve(match.group(1)), value)
        if isinstance(value, list):
            return [expand(item) for item in value]
        if isinstance(value, tuple):
            return tuple(expand(item) for item in value)
        if isinstance(value, dict):
            return {key: expand(item) for key, item in value.items()}
        return value

    return expand({key: value for key, value in raw.items() if key != "path_vars"})

def _frequency_range(values: Any, label: str) -> List[float]:
    if not isinstance(values, (list, tuple)) or len(values) != 2:
        raise ValueError(f"{label} must contain exactly two frequencies")
    frequencies = [float(value) for value in values]
    if not all(math.isfinite(value) for value in frequencies):
        raise ValueError(f"{label} must contain finite frequencies")
    if frequencies[0] >= frequencies[1]:
        raise ValueError(f"{label} must be increasing")
    return frequencies

def _shared_frequency_ranges(
    raw: Dict[str, Any],
    casa_raw: Dict[str, Any],
    low_high_raw: Dict[str, Any],
) -> FrequencyRangesCfg:
    shared_raw = raw.get("frequency_ranges")
    if shared_raw is not None:
        missing = [
            key for key in ("lowband_hz", "highband_hz") if key not in shared_raw
        ]
        if missing:
            raise ValueError(
                "frequency_ranges requires: " + ", ".join(missing)
            )
        source = shared_raw
    elif all(key in low_high_raw for key in ("lowband_hz", "highband_hz")):
        # Backward compatibility with the original step-specific schema.
        source = low_high_raw
    elif all(key in casa_raw for key in ("lowband_hz", "highband_hz")):
        source = casa_raw
    elif low_high_raw.get("enabled", False):
        raise ValueError(
            "Enabled low_high_slice requires frequency_ranges.lowband_hz and "
            "frequency_ranges.highband_hz"
        )
    else:
        return FrequencyRangesCfg()

    return FrequencyRangesCfg(
        lowband_hz=_frequency_range(source["lowband_hz"], "lowband_hz"),
        highband_hz=_frequency_range(source["highband_hz"], "highband_hz"),
    )

def _step_frequency_override(
    raw: Dict[str, Any],
    key: str,
    shared: List[float],
    section: str,
) -> Optional[List[float]]:
    if key not in raw:
        return None
    values = _frequency_range(raw[key], f"{section}.{key}")
    return None if values == shared else values

def _target(
    raw: Optional[Dict[str, Any]],
    legacy_slice: Optional[Dict[str, Any]] = None,
    default_name: str = "reference",
) -> Target:
    target_raw = dict(raw or {})
    legacy_raw = dict(legacy_slice or {})
    name = target_raw.get("name", default_name)
    legacy_name = legacy_raw.get("name")
    if legacy_name and legacy_name != name:
        raise ValueError(
            f"Legacy low_high_slice target name '{legacy_name}' does not match "
            f"target name '{name}'"
        )

    slice_paths: Dict[str, Optional[str]] = {}
    for key in SLICE_PATH_FIELDS:
        if key in target_raw and key in legacy_raw and target_raw[key] != legacy_raw[key]:
            raise ValueError(f"Conflicting {key} values for target '{name}'")
        slice_paths[key] = target_raw.get(key, legacy_raw.get(key))

    return Target(
        name=name,
        ms_paths=list(target_raw.get("ms_paths") or []),
        images=list(target_raw.get("images") or []),
        **slice_paths,
    )

def _dict_to_dataclass(d: Dict[str, Any]) -> Config:
    d = _expand_path_vars(d)
    paths = PathsCfg(**d.get("paths", {})) if "paths" in d else PathsCfg(
        raw_dir=d.get("raw_dir", "data/raw"),
        interim_dir=d.get("interim_dir", "data/interim"),
        processed_dir=d.get("processed_dir", "data/processed"),
        reports_dir=d.get("reports_dir", "data/reports"),
        sky_xmatches_dir=d.get("sky_xmatches_dir", "data/processed/Sky-CrossMatches"),
    )
    top_level_extra = {k: v for k, v in d.items() if k not in {
        "project_name","out_root","paths","raw_dir","interim_dir","processed_dir","reports_dir",
        "sky_xmatches_dir","reference","tests","frequency_ranges","casa","low_high_slice",
        "pybdsf","xmatch","extra"}}
    merged_extra = {**(d.get("extra") or {}), **top_level_extra}
    casa_raw = dict(d.get("casa") or {})
    low_high_raw = d.get("low_high_slice") or {}
    low_high_enabled = low_high_raw.get("enabled", False)
    frequency_ranges = _shared_frequency_ranges(d, casa_raw, low_high_raw)

    for key in ("exclude_fields", "calibration_exclude_fields", "imaging_exclude_fields"):
        values = casa_raw.get(key, [])
        if not isinstance(values, list) or any(isinstance(v, bool) or not isinstance(v, (str, int)) for v in values):
            raise ValueError("casa." + key + " must be a list of field names or IDs")
    if not isinstance(casa_raw.get("imaging_enabled", True), bool):
        raise ValueError("casa.imaging_enabled must be a boolean")
    if not isinstance(casa_raw.get("quality_check", False), bool):
        raise ValueError("casa.quality_check must be a boolean")
    for key, default in (("stage_timeout_seconds", 7200), ("timeout_seconds", 21600), ("shutdown_timeout_seconds", 30),
                         ("quality_min_solution_fraction", 0.95), ("quality_max_residual", 0.1)):
        value = float(casa_raw.get(key, default))
        if not math.isfinite(value) or value <= 0 or (key == "quality_min_solution_fraction" and value > 1):
            raise ValueError("Invalid casa." + key)
    casa = CasaCfg(
        quality_check=casa_raw.get("quality_check", False),
        quality_min_solution_fraction=float(casa_raw.get("quality_min_solution_fraction", 0.95)),
        quality_max_residual=float(casa_raw.get("quality_max_residual", 0.1)),
        stage_timeout_seconds=float(casa_raw.get("stage_timeout_seconds", 7200)),
        timeout_seconds=float(casa_raw.get("timeout_seconds", 21600)),
        shutdown_timeout_seconds=float(casa_raw.get("shutdown_timeout_seconds", 30)),
        imaging_enabled=casa_raw.get("imaging_enabled", True),
        exclude_fields=[str(v) for v in casa_raw.get("exclude_fields", [])],
        calibration_exclude_fields=[str(v) for v in casa_raw.get("calibration_exclude_fields", [])],
        imaging_exclude_fields=[str(v) for v in casa_raw.get("imaging_exclude_fields", [])],
        refant=str(casa_raw.get("refant") or "auto"),
        flux_field=str(casa_raw.get("flux_field", "J0408-6545")),
        scans=casa_raw.get("scans", ""),
        lowband_hz=_step_frequency_override(
            casa_raw, "lowband_hz", frequency_ranges.lowband_hz, "casa"
        ),
        highband_hz=_step_frequency_override(
            casa_raw, "highband_hz", frequency_ranges.highband_hz, "casa"
        ),
        field=casa_raw.get("field", "J2147-8132"),
        datacolumn=casa_raw.get("datacolumn", "DATA"),
        image_scans=casa_raw.get("image_scans", True),
    )
    low_high_slice = LowHighSliceCfg(
        enabled=low_high_enabled,
        overwrite=low_high_raw.get("overwrite", False),
        add_to_source_finding=low_high_raw.get("add_to_source_finding", True),
        lowband_hz=_step_frequency_override(
            low_high_raw, "lowband_hz", frequency_ranges.lowband_hz, "low_high_slice"
        ),
        highband_hz=_step_frequency_override(
            low_high_raw, "highband_hz", frequency_ranges.highband_hz, "low_high_slice"
        ),
    )

    legacy_tests: Dict[str, Dict[str, Any]] = {}
    for entry in low_high_raw.get("tests", []):
        name = entry.get("name")
        if not name:
            raise ValueError("Each legacy low_high_slice tests entry requires a name")
        if name in legacy_tests:
            raise ValueError(f"Duplicate legacy low_high_slice test name: {name}")
        legacy_tests[name] = entry

    tests = [
        _target(test_raw, legacy_tests.pop(test_raw.get("name"), None), default_name="")
        for test_raw in d.get("tests", [])
    ]
    if legacy_tests:
        raise ValueError(
            "Legacy low_high_slice test(s) have no matching tests[].name entry: "
            + ", ".join(sorted(legacy_tests))
        )

    cfg = Config(
        project_name=d["project_name"],
        out_root=d.get("out_root", "data"),
        paths=paths,
        reference=_target(
            d.get("reference"), low_high_raw.get("reference"), default_name="reference"
        ),
        tests=tests,
        frequency_ranges=frequency_ranges,
        casa=casa,
        low_high_slice=low_high_slice,
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
