from __future__ import annotations
import sys
import json
import copy
from typing import Any, Iterable
from ..audit import (
    raise_for_failures,
    record_failure,
    record_skip,
    register_inputs,
    run_logged_command,
)
from ..config import Config
from ..thermal_noise import band_from_label
from ..sensitivity_workflow import json_path, configured_products, load, absolute


def _thermal_inputs(cfg: Config, analysis: dict, *, infer: bool) -> dict:
    """Infer only an unambiguous one-reference/one-test MS association."""
    configured = copy.deepcopy(analysis.get("thermal_noise") or {})
    if not infer or len(cfg.tests) != 1:
        return configured
    for product in ("low", "high", "mfs"):
        sides = configured.setdefault(product, {})
        for role, target in (("reference", cfg.reference), ("test", cfg.tests[0])):
            if role in sides:
                continue
            if len(target.ms_paths) == 1:
                spec = dict(ms=target.ms_paths[0], field=cfg.casa.field,
                            scans=cfg.casa.scans, datacolumn=cfg.casa.datacolumn)
            elif not target.ms_paths:
                visibility = (cfg.extra.get('verification_report') or {}).get(role + '_visibility_results')
                if not visibility or analysis.get("allow_approximate") is not True:
                    continue
                spec = dict(visibility_results=visibility, scans=cfg.casa.scans, allow_approximate=True)
            else:
                continue
            band = band_from_label(target.name)
            if band:
                spec['band'] = band
            if product != "mfs":
                spec['frequency_range_hz'] = getattr(cfg.frequency_ranges, product + 'band_hz')
            sides[role] = spec
    return configured

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
    Optional: scans_glob, docx_name, thermal_noise (per-product reference/test inputs).
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
            artifact = json_path(cfg, analysis)
            if artifact:
                if analysis.get("thermal_noise") is not None:
                    raise ValueError("Conflicting sensitivity_json and direct thermal_noise inputs")
                products = configured_products(cfg)
                load(artifact, products)
                # Crossmatch paths must be paired with these explicit image associations.
                for product in ("low", "high", "mfs"):
                    entry = next((p for p in cfg.extra.get("positions", [])
                                  if absolute(p["xmatch_table"]) == absolute(analysis[f"ref_{product}_xmatch"])), None)
                    if entry is None or any(absolute(entry[k]) != absolute(products[product][r]["image"])
                        for k,r in (("ref_fits", "reference"), ("other_fits", "test"))):
                        raise ValueError("Flux crossmatch/image sensitivity association mismatch: " + product)
                cmd += ["--sensitivity-json", artifact, "--sensitivity-products-json", json.dumps(products)]
            else:
                thermal = _thermal_inputs(cfg, analysis, infer=len(analyses) == 1)
                cmd += ["--thermal-noise-json", json.dumps(thermal)]
            for key, value in (cfg.extra.get("uncertainty") or {}).items():
                cmd += ["--" + key.replace("_", "-"), str(value)]
            _run(cmd, inputs=[label])
        except Exception as exc:
            record_failure(label, exc)
            failures.append(exc)

    raise_for_failures("flux analysis", failures)
