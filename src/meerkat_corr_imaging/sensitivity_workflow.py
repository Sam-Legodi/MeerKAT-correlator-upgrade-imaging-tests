"""Reusable, versioned sensitivity artifacts and explicit existing-image handoffs.

Only ``generate`` opens MeasurementSets. Consumers validate the stored selection
and image identities, never remeasure visibilities. File-stat identities are
provenance snapshots, not cryptographic hashes of MS data content.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import tempfile
from datetime import datetime, timezone

from .thermal_noise import sensitivity, FORMULA, ASSUMPTIONS

SCHEMA_VERSION = 1
CALCULATION_VERSION = "natural-stokes-i-exposure-1"
EXPOSURE_DEFINITION = (
    "C=sum(abs(CHAN_WIDTH_Hz)*EXPOSURE_s) over selected cross-correlation "
    "row/channel/parallel-hand cells with clear FLAG/FLAG_ROW, finite selected "
    "data and finite positive WEIGHT_SPECTRUM (or WEIGHT if undefined). "
    "Each stored baseline and each XX/YY or RR/LL hand counted once; "
    "C=N*(N-1)*effective_bandwidth_hz*on_source_integration_s. Units Hz s."
)
ROLES = ("reference", "test")


def absolute(path):
    return str(Path(path).expanduser().resolve())


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                      allow_nan=False).encode()).hexdigest()


def file_identity(path):
    p = Path(path)
    if not p.is_file():
        return None
    stat = p.stat()
    return dict(path=absolute(p), size=stat.st_size, mtime_ns=stat.st_mtime_ns)


def association(spec):
    """Canonical requested association, including the actual manifest content."""
    value = dict(spec)
    for key in ("image", "ms", "manifest", "non_pb_image", "visibility_results"):
        if value.get(key):
            value[key] = absolute(value[key])
    value.pop("catalogue", None)  # catalogue location does not affect sensitivity
    if value.get("manifest"):
        value["manifest_sha256"] = hashlib.sha256(Path(value["manifest"]).read_bytes()).hexdigest()
    if isinstance(value.get("channels"), dict):
        value["channels"] = {str(k): sorted(set(v)) for k, v in value["channels"].items()}
    if isinstance(value.get("scans"), list):
        value["scans"] = sorted(set(value["scans"]))
    return value


def configured_products(cfg):
    products = (cfg.extra.get("sensitivity") or {}).get("products")
    if not isinstance(products, dict) or not products:
        raise ValueError("extra.sensitivity.products must explicitly associate image products with selections")
    for product, sides in products.items():
        if not isinstance(product, str) or not product or set(sides) != set(ROLES):
            raise ValueError("Each sensitivity product requires reference and test mappings")
        if any(not isinstance(sides[r], dict) for r in ROLES):
            raise ValueError("Sensitivity observation inputs must be mappings")
    return products


def json_path(cfg, analysis=None):
    paths = [(cfg.extra.get("sensitivity") or {}).get("output_json"),
             cfg.extra.get("sensitivity_json"), (analysis or {}).get("sensitivity_json")]
    paths = {absolute(p) for p in paths if p}
    if len(paths) > 1:
        raise ValueError("Conflicting sensitivity JSON paths")
    return next(iter(paths), None)


def resolved_spec(spec):
    """Require explicit selections, or the completed paired-imaging contract."""
    for key in ("image", "band", "pb_corrected"):
        if key not in spec or spec[key] is None:
            raise ValueError(f"Explicit {key} association required")
    if type(spec["pb_corrected"]) is not bool:
        raise ValueError("pb_corrected must be a boolean")
    value = dict(spec)
    if spec.get("manifest"):
        manifest = json.loads(Path(spec["manifest"]).read_text())
        if manifest.get("processing_status") != "complete":
            raise ValueError("Imaging manifest must be complete")
        contract = manifest["contract"]
        if absolute(contract["output_paths"]["fits"]) != absolute(spec["image"]):
            raise ValueError("Imaging manifest image association mismatch")
        params = contract["tclean_parameters"]
        derived = dict(ms=absolute(params["vis"]), field=contract["field_id"],
                       scans=contract["scan_ids"], datacolumn=params["datacolumn"], channels={})
        for channel in contract["requested_channels"]:
            derived["channels"].setdefault(str(channel["spw"]), []).append(channel["channel"])
        for key, val in derived.items():
            if key in spec:
                other = absolute(spec[key]) if key == "ms" else spec[key]
                if other != val:
                    raise ValueError("Conflicting explicit and manifest selection: " + key)
        value.update(derived)
        value.pop("manifest", None)
        value["imaging_context"] = {k: params[k] for k in
            ("weighting", "robust", "uvtaper", "gridder", "deconvolver") if k in params}
        value["manifest_ms_identity"] = contract["ms_identity"]
    for key in ("ms", "field", "datacolumn", "scans", "channels"):
        if key not in value or value[key] is None or value[key] == "" or value[key] == [] or value[key] == {}:
            raise ValueError(f"Explicit {key} selection required (use 'all' for all scans/channels)")
    if not isinstance(value["channels"], dict) and value["channels"] != "all":
        raise ValueError("channels must be 'all' or a SPW-to-channel-index mapping")
    if str(value["datacolumn"]).lower() not in ("data", "corrected", "corrected_data"):
        raise ValueError("datacolumn must be DATA or CORRECTED_DATA")
    return value


def validate_pb_header(spec):
    """Check an explicit header if present; never guess from the image name."""
    from astropy.io import fits
    for key, expected in (("image", spec.get("pb_corrected")), ("non_pb_image", False)):
        if spec.get(key) and Path(spec[key]).is_file():
            header = fits.getheader(spec[key])
            if "PBCOR" in header:
                val = header["PBCOR"]
                if isinstance(val, str):
                    val = val.strip().lower() in ("true", "t", "yes", "1")
                if bool(val) != expected:
                    raise ValueError("PBCOR header contradicts explicit image state: " + str(spec[key]))


def generate(cfg):
    from .paired_astrometry import ms_identity
    products = configured_products(cfg)
    from .audit import register_inputs
    register_inputs([str(s.get("manifest") or s.get("ms") or "missing MS") for v in products.values() for s in v.values()])
    output = json_path(cfg)
    if not output:
        raise ValueError("extra.sensitivity.output_json is required")
    payload = dict(schema_version=SCHEMA_VERSION, calculation_version=CALCULATION_VERSION,
                   created_utc=datetime.now(timezone.utc).isoformat(), units="uJy/beam", products={})
    for product, sides in products.items():
        payload["products"][product] = {}
        for role, spec in sides.items():
            # Missing manifests are recorded unavailable, not silently replaced by guessed selections.
            try:
                requested = association(spec)
            except OSError:
                requested = dict(spec)
            result = dict(status="unavailable", approximation_status="unavailable", reason=None,
                theoretical_rms_ujy_beam=None, antenna_count=None, antenna_ids=[],
                effective_bandwidth_hz=None, on_source_integration_s=None,
                parallel_hand_exposure_hz_s=None, ms_identity=None, selection=None,
                image=absolute(spec["image"]) if spec.get("image") else None,
                pb_corrected=spec.get("pb_corrected"), band=spec.get("band"),
                sefd_jy=None, sefd_reference=None, formula=FORMULA, units="uJy/beam",
                assumptions=ASSUMPTIONS, calculation_version=CALCULATION_VERSION,
                exposure_definition=EXPOSURE_DEFINITION, association=requested,
                association_sha256=digest(requested),
                image_identity=file_identity(spec["image"]) if spec.get("image") else None,
                non_pb_image_identity=file_identity(spec["non_pb_image"]) if spec.get("non_pb_image") else None)
            try:
                selected = resolved_spec(spec)
                validate_pb_header(spec)
                if Path(absolute(output)).is_relative_to(Path(absolute(selected["ms"]))):
                    raise ValueError("Sensitivity output must not be inside the MS")
                approximate = bool(selected.get("visibility_results"))
                before = None if approximate else ms_identity(selected["ms"])
                if selected.get("manifest_ms_identity") and selected["manifest_ms_identity"] != before:
                    raise ValueError("MS identity differs from completed imaging manifest")
                inputs = dict(selected)
                if approximate:
                    inputs.pop("ms", None)
                    # Cache must be explicitly associated with this field/all scans.
                    if inputs["scans"] != "all" or inputs["channels"] != "all":
                        raise ValueError("Cached summaries require all scans/channels; frequency cuts alone may be approximated")
                    inputs.pop("scans"); inputs.pop("channels")
                computed = sensitivity(inputs)
                if not approximate and before != ms_identity(selected["ms"]):
                    raise ValueError("MS changed during sensitivity calculation")
                result.update(computed)
                result.update(ms_identity=before or dict(path=absolute(selected["ms"]),
                    unavailable_reason="MS not inspected; explicitly opted-in cached approximation"),
                    approximation_status="approximate" if approximate else "exact_selected_exposure",
                    selection={k: selected[k] for k in ("ms", "field", "datacolumn", "scans", "channels")},
                    imaging_context=selected.get("imaging_context", {}), reason=None)
            except (ValueError, OSError, ImportError, RuntimeError, KeyError) as exc:
                result["reason"] = str(exc)
            payload["products"][product][role] = result
            print(f"[SENSITIVITY] {product}/{role}: {result['status']} ({result['approximation_status']})" +
                  (f"; {result['reason']}" if result["reason"] else ""))
    validate(payload)
    atomic_write(output, payload)
    print(f"[SENSITIVITY] Wrote {output}")
    return payload


def validate(payload):
    if not isinstance(payload, dict):
        raise ValueError("Sensitivity payload must be a mapping")
    if type(payload.get("schema_version")) is not int or payload.get("schema_version") != SCHEMA_VERSION or payload.get("calculation_version") != CALCULATION_VERSION:
        raise ValueError("Incompatible sensitivity schema/calculation version")
    if payload.get("units") != "uJy/beam" or not isinstance(payload.get("products"), dict) or not payload["products"]:
        raise ValueError("Invalid sensitivity units/products")
    required = ("status", "approximation_status", "reason", "theoretical_rms_ujy_beam", "antenna_count",
        "antenna_ids", "effective_bandwidth_hz", "on_source_integration_s", "parallel_hand_exposure_hz_s",
        "ms_identity", "selection", "image", "pb_corrected", "band", "sefd_jy", "sefd_reference",
        "formula", "units", "calculation_version", "assumptions", "exposure_definition",
        "association", "association_sha256", "image_identity", "non_pb_image_identity")
    for sides in payload["products"].values():
        if not isinstance(sides, dict):
            raise ValueError("Sensitivity product must be a mapping")
        if set(sides) != set(ROLES):
            raise ValueError("Sensitivity products require reference/test")
        for result in sides.values():
            if not isinstance(result, dict):
                raise ValueError("Sensitivity result must be a mapping")
            if any(k not in result for k in required):
                raise ValueError("Incomplete sensitivity result schema")
            if result["association_sha256"] != digest(result["association"]):
                raise ValueError("Sensitivity association digest mismatch")
            if result["exposure_definition"] != EXPOSURE_DEFINITION or result["formula"] != FORMULA or result["units"] != "uJy/beam" or result["calculation_version"] != CALCULATION_VERSION:
                raise ValueError("Incompatible sensitivity formula/exposure/units")
            if result["status"] == "available":
                if type(result["pb_corrected"]) is not bool or not result["selection"] or not result["ms_identity"] or not result["sefd_reference"]:
                    raise ValueError("Available sensitivity lacks association/provenance")
                if result["approximation_status"] not in ("exact_selected_exposure", "approximate"):
                    raise ValueError("Invalid approximation status")
                n, bw, t = (result[k] for k in ("antenna_count", "effective_bandwidth_hz", "on_source_integration_s"))
                if isinstance(n, bool) or not isinstance(n, int) or n < 2 or not all(
                    isinstance(v, (int, float)) and math.isfinite(v) and v > 0 for v in
                    (bw, t, result["sefd_jy"], result["theoretical_rms_ujy_beam"], result["parallel_hand_exposure_hz_s"])):
                    raise ValueError("Invalid sensitivity quantities")
                if len(set(result["antenna_ids"])) != n:
                    raise ValueError("Participating antenna IDs do not match count")
                if not math.isclose(n*(n-1)*bw*t, result["parallel_hand_exposure_hz_s"], rel_tol=1e-10):
                    raise ValueError("Inconsistent sensitivity calculation/exposure")
            elif result["status"] != "unavailable" or not result["reason"] or result["theoretical_rms_ujy_beam"] is not None:
                raise ValueError("Invalid unavailable sensitivity result")
    return payload


def atomic_write(path, payload):
    p = Path(path); p.parent.mkdir(parents=True, exist_ok=True)
    temp = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", dir=p.parent, prefix=p.name+".", suffix=".tmp", delete=False) as f:
            temp = f.name
            json.dump(payload, f, indent=2, allow_nan=False); f.write("\n")
            f.flush(); os.fsync(f.fileno())
        os.replace(temp, p)
    finally:
        if temp and os.path.exists(temp):
            os.unlink(temp)


def load(path, expected_products):
    """Validate explicit associations and FITS identities without any MS access."""
    payload = validate(json.loads(Path(path).read_text()))
    if not expected_products:
        raise ValueError("Explicit sensitivity product associations required downstream")
    if set(expected_products) != set(payload["products"]):
        raise ValueError("Sensitivity product association mismatch")
    for product, sides in expected_products.items():
        for role, spec in sides.items():
            result = payload["products"][product][role]
            try:
                expected = association(spec)
            except OSError:
                expected = dict(spec)
            if digest(expected) != result["association_sha256"]:
                raise ValueError(f"Sensitivity image/selection association mismatch: {product}/{role}")
            if result["status"] == "available":
                selected = resolved_spec(spec)
                if result["selection"] != {k:selected[k] for k in ("ms", "field", "datacolumn", "scans", "channels")}:
                    raise ValueError("Stored sensitivity selection mismatch")
                if absolute(result["ms_identity"]["path"]) != absolute(selected["ms"]) or result["image"] != absolute(spec["image"]) or result["pb_corrected"] != spec["pb_corrected"] or result["band"] != str(spec["band"]).upper():
                    raise ValueError("Stored sensitivity image/MS/band association mismatch")
            for key, identity_key in (("image", "image_identity"), ("non_pb_image", "non_pb_image_identity")):
                old = result[identity_key]
                if old and old != file_identity(spec[key]):
                    raise ValueError("Image changed since sensitivity snapshot: " + spec[key])
            validate_pb_header(spec)
    return payload


def exposure_ratio(sides):
    reference, test = (sides[r] for r in ROLES)
    if any(r["status"] != "available" for r in (reference, test)):
        return None
    if reference["band"] != test["band"] or reference["sefd_jy"] != test["sefd_jy"]:
        raise ValueError("Exposure scaling requires equal band/SEFD assumptions")
    if reference["approximation_status"] != test["approximation_status"]:
        raise ValueError("Do not mix exact MS and cached approximate exposures")
    return reference["parallel_hand_exposure_hz_s"] / test["parallel_hand_exposure_hz_s"]


def noise_comparison(sides, reference_rms, test_rms):
    ratio = exposure_ratio(sides)
    def multiply(value):
        return 1e6*value if value is not None and math.isfinite(value) else None
    ref, gpu = multiply(reference_rms), multiply(test_rms)
    expected = ref*math.sqrt(ratio) if ref is not None and ratio is not None else None
    rt, gt = (sides[r]["theoretical_rms_ujy_beam"] for r in ROLES)
    def divide(a, b):
        return a/b if a is not None and b is not None and b > 0 else None
    pb = any(sides[r]["pb_corrected"] is True for r in ROLES)
    return dict(measured_gpu_ujy_beam=gpu, expected_gpu_ujy_beam=expected,
        thermal_gpu_ujy_beam=gt, thermal_reference_ujy_beam=rt,
        measured_gpu_over_thermal_gpu=divide(gpu,gt), expected_gpu_over_thermal_gpu=divide(expected,gt),
        measured_reference_over_thermal_reference=divide(ref,rt), measured_reference_ujy_beam=ref,
        exposure_ratio_reference_over_test=ratio, exposure_definition=EXPOSURE_DEFINITION,
        measurement_region="phase-centre annulus 0.25–0.50 deg; 1.4826 * median absolute deviation",
        qualified_pb_comparison=pb,
        qualification=("PB-corrected annular noise is not equivalent to uncorrected on-axis thermal sensitivity; no beam correction applied."
                       if pb else "Uncorrected annular image noise versus on-axis natural-weight thermal expectation; weighting/taper and image systematics may differ."),
        sensitivity_inputs=sides, classification="diagnostic only; no thermal Pass/Concern threshold")


def wire_noise_pipeline(cfg):
    """Wire supplied FITS products; never schedule calibration/imaging/slicing."""
    from .image_pipeline import _catalogue_for_image, _slug
    if len(cfg.tests) != 1:
        raise ValueError("noise-report currently requires exactly one reference/GPU pair")
    products = configured_products(cfg)
    if len({_slug(p) for p in products}) != len(products):
        raise ValueError("Image product names must have unique filename slugs")
    if any(cfg.extra.get(k) for k in ("xmatch_pairs", "positions", "flux")):
        raise ValueError("noise-report derives xmatch/positions/flux from sensitivity.products; remove conflicting handoffs")
    pairs, positions, tables = [], [], {}
    cfg.reference.images = []; cfg.tests[0].images = []
    for product, sides in products.items():
        images = {r: absolute(sides[r]["image"]) for r in ROLES}
        if sides["reference"].get("pb_corrected") != sides["test"].get("pb_corrected"):
            raise ValueError("Paired image PB correction states must match")
        output = str(Path(cfg.paths.sky_xmatches_dir) /
                     f"{_slug(cfg.reference.name)}_X_{_slug(cfg.tests[0].name)}_{_slug(product)}.fits")
        pairs.append([_catalogue_for_image(cfg, images[r]) for r in ROLES]+[output])
        positions.append(dict(xmatch_table=output, ref_fits=images["reference"], other_fits=images["test"],
                              otherdatatag=cfg.tests[0].name+"_"+product, product=product))
        tables[product] = output
        cfg.reference.images.append(images["reference"]); cfg.tests[0].images.append(images["test"])
    cfg.low_high_slice.enabled = False
    cfg.extra.update(xmatch_pairs=pairs, positions=positions, images_globs=[])
    cfg.extra["flux"] = ([{f"ref_{p}_xmatch":tables[p] for p in ("low", "high", "mfs")}]
                         if {"low", "high", "mfs"} <= set(tables) else [])
    return positions
