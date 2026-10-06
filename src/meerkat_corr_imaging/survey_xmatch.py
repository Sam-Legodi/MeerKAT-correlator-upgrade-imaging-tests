"""Explicit local survey jobs; positional associations only, without flux analysis."""
from __future__ import annotations

from copy import deepcopy
from dataclasses import asdict, dataclass, field
import hashlib
import json
import math
from pathlib import Path
import re
from typing import Any

import astropy.units as u
from astropy.io import fits

from .xmatch_pybdsf import (
    _one_to_many, cross_match_catalogues, parse_angle, prepare_catalog,
    prepare_output_directory, write_output,
)


# Only inspected schemas are built in. Missing newer RACS releases require an
# explicit user profile after their actual source tables have been inspected.
PROFILES = {
    "pybdsf_gaul": {
        "survey": "MeerKAT", "release": "PyBDSF gaul", "grain": "component",
        "ra_col": "RA", "dec_col": "DEC", "frame": "icrs",
        "coordinate_unit": "deg", "source_id": "Source_id", "component_id": "Gaus_id",
        "required_columns": ["Source_id", "Gaus_id"],
        "qualification": "Pipeline source-cat filenames contain Gaussian components (catalog_type=gaul).",
    },
    "pybdsf_source": {
        "survey": "MeerKAT", "release": "PyBDSF", "grain": "source",
        "ra_col": "RA", "dec_col": "DEC", "frame": "icrs",
        "coordinate_unit": "deg", "source_id": "Source_id",
        "required_columns": ["Source_id"],
    },
    "sumss_v21r": {
        "survey": "SUMSS", "release": "V2.1r", "grain": "source",
        "ra_col": "_RAJ2000", "dec_col": "_DEJ2000", "frame": "fk5",
        "coordinate_unit": "deg", "frequency_hz": 843e6,
        "required_columns": ["St", "e_St", "Sp", "e_Sp"],
        "qualification": "Sp/e_Sp are mJy/beam; VizieR labels them mJy. St/e_St are mJy.",
        "documentation": "https://cdsarc.cds.unistra.fr/viz-bin/ReadMe/VIII/81B?format=html",
    },
    "racs_low_dr1_gausscut": {
        "survey": "RACS-low", "release": "DR1 J/other/PASA/38.58/gausscut",
        "grain": "component", "ra_col": "RAJ2000", "dec_col": "DEJ2000",
        "frame": "icrs", "coordinate_unit": "deg", "frequency_hz": 887.5e6,
        "component_id": "GID", "source_id": "ID", "required_columns": ["GID", "ID"],
        "qualification": "Gaussian rows retained; parent Ftot is repeated and must not be summed.",
    },
}


def absolute(path: str | Path) -> str:
    return str(Path(path).expanduser().resolve(strict=False))


@dataclass(frozen=True)
class CatalogueSpec:
    path: str
    survey: str
    release: str
    grain: str
    ra_col: str
    dec_col: str
    frame: str
    coordinate_unit: str
    format: str | None = None
    hdu: int | str | None = None
    table_id: int | str | None = None
    frequency_hz: float | None = None
    image: str | None = None
    source_id: str | None = None
    component_id: str | None = None
    required_columns: list[str] = field(default_factory=list)
    qualification: str = ""
    documentation: str = ""

    @classmethod
    def from_config(cls, raw: dict, profiles: dict | None = None) -> CatalogueSpec:
        raw = deepcopy(raw)
        name = raw.pop("profile", None)
        available = {**PROFILES, **(profiles or {})}
        if name is not None and name not in available:
            raise ValueError(f"Unknown catalogue profile: {name}")
        values = {**deepcopy(available.get(name, {})), **raw}
        values["path"] = absolute(values["path"])
        if values.get("image"):
            values["image"] = absolute(values["image"])
        spec = cls(**values)
        if spec.grain not in {"source", "component"}:
            raise ValueError("Catalogue grain must be source or component")
        if not spec.survey or not spec.release:
            raise ValueError("Catalogue survey and release must be explicit")
        if spec.grain == "component" and not spec.component_id:
            raise ValueError("Component catalogue requires component_id")
        if spec.frequency_hz is not None:
            frequency_hz = float(spec.frequency_hz)
            if not math.isfinite(frequency_hz) or frequency_hz <= 0:
                raise ValueError("Catalogue frequency_hz must be positive and finite")
            values["frequency_hz"] = frequency_hz
            spec = cls(**values)
        # INITIAL data retain the preliminary qualification even in custom profiles.
        if ("INITIAL" in (spec.release + Path(spec.path).name).upper()
                and "Preliminary INITIAL measurements." not in spec.qualification):
            values["qualification"] = "Preliminary INITIAL measurements. " + spec.qualification
            spec = cls(**values)
        return spec

    def prepare(self, prefix):
        required = list(self.required_columns)
        required.extend(c for c in (self.source_id, self.component_id) if c)
        catalog = prepare_catalog(
            self.path, self.ra_col, self.dec_col, self.frame, prefix,
            format=self.format, hdu=self.hdu, table_id=self.table_id,
            coordinate_unit=self.coordinate_unit, required_columns=required,
            strict_table_selection=True,
        )
        if self.grain == "source" and self.source_id:
            ids = catalog.table[self.source_id]
            ids = ids.compressed() if hasattr(ids, "compressed") else ids
            if len(set(ids)) != len(ids):
                raise ValueError(f"Repeated source IDs in source-grain input {self.path}; "
                                 "declare component grain and component_id for Gaussian rows")
        return catalog


@dataclass(frozen=True)
class SurveyJob:
    input1: CatalogueSpec
    input2: CatalogueSpec
    output: str
    radius: str
    association_mode: str = "all-candidates"
    coverage: str = "unknown"
    coverage_evidence: str = ""

    @classmethod
    def from_config(cls, raw: dict, profiles: dict | None = None) -> SurveyJob:
        parse_angle(raw["radius"])
        mode = raw.get("association_mode", "all-candidates")
        if mode not in {"all-candidates", "one-to-one"}:
            raise ValueError(f"Invalid association_mode: {mode}")
        coverage = raw.get("coverage", "unknown")
        evidence = raw.get("coverage_evidence", "")
        if coverage not in {"unknown", "covered", "absent"}:
            raise ValueError(f"Invalid coverage: {coverage}")
        if coverage != "unknown" and not evidence:
            raise ValueError("Coverage claims require independent coverage_evidence")
        output = Path(absolute(raw["output"]))
        if output.parent.name != "Sky-CrossMatches":
            output = output.parent / "Sky-CrossMatches" / output.name
        return cls(CatalogueSpec.from_config(raw["input1"], profiles),
                   CatalogueSpec.from_config(raw["input2"], profiles), str(output),
                   raw["radius"], mode, coverage, evidence)

    def signature(self) -> str:
        values = asdict(self)
        values.pop("output")
        # Canonicalize equivalent angular quantities for deduplication.
        values["radius"] = float(parse_angle(self.radius).arcsec)
        return json.dumps(values, sort_keys=True)


def resolve_survey_jobs(cfg, products=None) -> list[dict[str, Any]]:
    """Expand selected image products for both standalone xm and image wiring.

    Explicit jobs take priority over equivalent generated jobs. Path resolution
    collapses the commissioning directory symlink before deduplication.
    """
    profiles = cfg.extra.get("catalogue_profiles", {})
    jobs = deepcopy(cfg.extra.get("survey_xmatch_jobs") or [])
    templates = cfg.extra.get("survey_xmatches") or []
    if templates and products is None:
        from .image_pipeline import _image_products
        products = _image_products(cfg)
    for template in templates:
        if not template.get("enabled", True):
            continue
        selected_names = template.get("targets")
        bands = template.get("bands", ["mfs"])
        for product in products or []:
            if product.band not in bands:
                continue
            if selected_names is not None and product.target_name not in selected_names:
                continue
            raw = {key: deepcopy(template[key]) for key in
                   ("radius", "association_mode", "coverage", "coverage_evidence") if key in template}
            raw["input1"] = {"profile": template.get("input_profile", "pybdsf_gaul"), "path": product.catalogue,
                             "image": product.image}
            raw["input2"] = deepcopy(template["catalogue"])
            # Include settings in the output name so separate radii/modes/releases
            # cannot overwrite one another, even for similarly named products.
            digest_job = SurveyJob.from_config({**raw, "output": "unused.fits"}, profiles)
            digest = hashlib.sha256(digest_job.signature().encode()).hexdigest()[:12]
            label = f"{product.target_name}_{product.band}_X_{template.get('name', 'survey')}"
            slug = re.sub(r"[^A-Za-z0-9._-]+", "_", label)[:120]
            raw["output"] = str(Path(cfg.paths.sky_xmatches_dir) / f"{slug}_{digest}.fits")
            jobs.append(raw)
    unique = {}
    outputs = {}
    for raw in jobs:
        if not raw.get("enabled", True):
            continue
        normalized = {key: value for key, value in raw.items() if key != "enabled"}
        job = SurveyJob.from_config(normalized, profiles)
        signature = job.signature()
        if signature in unique:
            continue
        if job.output in outputs and outputs[job.output] != signature:
            raise ValueError(f"Different survey jobs request the same output: {job.output}")
        outputs[job.output] = signature
        unique[signature] = asdict(job)
    return list(unique.values())


def frequency_provenance(spec: CatalogueSpec, table) -> dict:
    """Report the FITS observing frequency without inventing effective frequencies."""
    result = {"frequency_hz": spec.frequency_hz,
              "origin": "configured" if spec.frequency_hz is not None else "unknown"}
    image = spec.image
    if not image and table.meta.get("INIMAGE"):
        candidate = Path(str(table.meta["INIMAGE"]))
        if not candidate.is_absolute():
            candidate = Path(spec.path).parent / candidate
        if candidate.is_file():
            image = str(candidate)
    if not image:
        return result
    # An explicitly requested image must exist and be readable.
    header = fits.getheader(image)
    result["image"] = absolute(image)
    result["configured_frequency_hz"] = spec.frequency_hz
    for axis in range(1, int(header.get("WCSAXES", header.get("NAXIS", 0))) + 1):
        if str(header.get(f"CTYPE{axis}", "")).upper().startswith("FREQ"):
            key = f"CRVAL{axis}"
            if key in header:
                hz = (float(header[key]) * u.Unit(header.get(f"CUNIT{axis}", "Hz"))).to_value(u.Hz)
                if math.isfinite(hz) and hz > 0:
                    result.update(frequency_hz=hz, origin=f"FITS {key} spectral reference")
                    break
    else:
        for key in ("FREQ", "CFREQ"):
            hz = float(header.get(key, float("nan")))
            if math.isfinite(hz) and hz > 0:
                result.update(frequency_hz=hz, origin=f"FITS {key}")
                break
    result["qualification"] = "Image header frequency; an effective weighted MFS frequency is not inferred."
    return result


def execute_survey_job(raw: dict) -> dict:
    job = SurveyJob.from_config(raw)
    output, _, _ = prepare_output_directory(job.output)
    diagnostics_path = Path(output).with_suffix(".diagnostics.json")
    if Path(output).exists() or diagnostics_path.exists():
        raise FileExistsError(f"Survey output already exists; choose a new output: {output}")
    cat1, cat2 = job.input1.prepare("T1"), job.input2.prepare("T2")
    table, matches = cross_match_catalogues(
        cat1, cat2, parse_angle(job.radius), job.association_mode == "all-candidates", enrich=False)
    frequencies = [frequency_provenance(spec, cat.table)
                   for spec, cat in ((job.input1, cat1), (job.input2, cat2))]
    diagnostics = {"job": asdict(job), "frequencies": frequencies,
                   "associations": len(matches), "catalogues": [],
                   "coverage": {"status": job.coverage, "evidence": job.coverage_evidence,
                                "qualification": "Catalogue row absence alone does not establish survey noncoverage."},
                   "association_qualification": "Positional candidates; secure physical counterparts are not inferred."}
    candidates = matches if job.association_mode == "all-candidates" else _one_to_many(
        cat1.coords, cat2.coords, cat1.valid_indices, cat2.valid_indices, parse_angle(job.radius))
    n1, n2 = {}, {}
    for a, b, _ in candidates:
        n1[a] = n1.get(a, 0) + 1
        n2[b] = n2.get(b, 0) + 1
    for side, spec, cat, freq in ((1, job.input1, cat1, frequencies[0]),
                                  (2, job.input2, cat2, frequencies[1])):
        counts = {key.lower(): table.meta[f"HIERARCH XM{side}_{key}"]
                  for key in ("INPUT", "VALID", "INVALID", "MATCHED", "UNMATCHED", "AMBIGUOUS")}
        counts["grain"] = spec.grain
        if spec.source_id:
            parents = cat.table[spec.source_id]
            counts["distinct_source_ids"] = len(set(parents.compressed() if hasattr(parents, "compressed") else parents))
            valid_parents = parents[cat.valid_indices]
            valid_parents = set(valid_parents.compressed() if hasattr(valid_parents, "compressed") else valid_parents)
            matched_parents = {parents[m[side - 1]] for m in matches
                               if not bool(getattr(parents[m[side - 1]], "mask", False))}
            # Ambiguity includes all candidates, even if one-to-one discards them.
            ambiguous_parents = {parents[m[side - 1]] for m in candidates
                                 if (n1[m[0]] > 1 or n2[m[1]] > 1)
                                 and not bool(getattr(parents[m[side - 1]], "mask", False))}
            counts.update(valid_source_ids=len(valid_parents), matched_source_ids=len(matched_parents),
                          unmatched_source_ids=len(valid_parents - matched_parents),
                          ambiguous_source_ids=len(ambiguous_parents))
            for key in ("valid_source_ids", "matched_source_ids", "unmatched_source_ids", "ambiguous_source_ids"):
                table.meta[f"HIERARCH XM{side}_{key.upper()}"] = counts[key]
        diagnostics["catalogues"].append(counts)
        # Full paths/settings also reside in HISTORY and the inspectable JSON.
        table.meta.setdefault("HISTORY", [])
        if isinstance(table.meta["HISTORY"], str):
            table.meta["HISTORY"] = [table.meta["HISTORY"]]
        provenance = json.dumps({"catalogue": asdict(spec), "frequency": freq}, sort_keys=True)
        table.meta["HISTORY"].extend(provenance[i:i+70] for i in range(0, len(provenance), 70))
        for key, value in {"SURVEY": spec.survey, "RELEASE": spec.release,
                           "GRAIN": spec.grain, "FRAME": spec.frame,
                           "CATALOG": Path(spec.path).name,
                           "FREQ_HZ": freq["frequency_hz"] if freq["frequency_hz"] is not None else "unknown"}.items():
            table.meta[f"HIERARCH XM{side}_{key}"] = value
    table.meta["HIERARCH XM_COVERAGE"] = job.coverage
    table.meta["HIERARCH XM_ASSOC_MODE"] = job.association_mode
    table.meta["HIERARCH XM_ASSOC_QUAL"] = "positional candidates"
    diagnostics["coverage"]["unmatched_input1"] = diagnostics["catalogues"][0]["unmatched"]
    diagnostics["coverage"]["unmatched_interpretation"] = {
        "unknown": "Unmatched rows have uncertain survey coverage; do not diagnose MeerKAT astrometry.",
        "absent": "Independent evidence reports absent coverage; counterpart search is not informative.",
        "covered": "No selected association within radius; sensitivity, morphology and ambiguity remain relevant.",
    }[job.coverage]
    write_output(table, output, overwrite=False)
    with diagnostics_path.open("x", encoding="utf-8") as handle:
        json.dump(diagnostics, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
    print(f"[XMATCH SURVEY] {json.dumps(diagnostics['catalogues'], sort_keys=True)}")
    print(f"[XMATCH SURVEY] coverage={job.coverage}; catalogue absence is not evidence of noncoverage.")
    return diagnostics
