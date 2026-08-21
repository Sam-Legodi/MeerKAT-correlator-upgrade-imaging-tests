from __future__ import annotations

from dataclasses import dataclass
import hashlib
from pathlib import Path
import re
from typing import Any, Iterable

from .config import Config, Target
from .pybdsf_srcfind import catalogue_paths, compute_base_name
from .steps.step4_low_high_slice import resolved_targets


@dataclass(frozen=True)
class ImageProduct:
    """An image and the deterministic PyBDSF catalogue produced from it."""

    role: str
    target_name: str
    band: str
    image: str
    catalogue: str


@dataclass(frozen=True)
class XMatchProduct:
    """A reference/test image pair and its cross-match table."""

    target_name: str
    band: str
    reference: ImageProduct
    test: ImageProduct
    table: str


@dataclass(frozen=True)
class ImagePipelinePlan:
    """Inputs handed from source finding through the analysis steps."""

    products: tuple[ImageProduct, ...]
    matches: tuple[XMatchProduct, ...]
    xmatch_pairs: tuple[tuple[str, ...], ...]
    positions: tuple[dict[str, Any], ...]
    flux: tuple[dict[str, Any], ...]


def _absolute(path: str | Path) -> str:
    return str(Path(path).expanduser().resolve(strict=False))


def _slug(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("._-") or "target"
    if len(cleaned) <= 64:
        return cleaned
    digest = hashlib.sha1(cleaned.encode("utf-8")).hexdigest()[:8]
    return f"{cleaned[:55]}-{digest}"


def _catalogue_for_image(cfg: Config, image: str) -> str:
    image_path = _absolute(image)
    base = compute_base_name(Path(image_path).stem, cfg.pybdsf.base_prefix)
    fits_catalogue, _ = catalogue_paths(image_path, base)
    return _absolute(fits_catalogue)


def _product(
    cfg: Config,
    *,
    role: str,
    target_name: str,
    band: str,
    image: str,
) -> ImageProduct:
    image_path = _absolute(image)
    return ImageProduct(
        role=role,
        target_name=target_name,
        band=band,
        image=image_path,
        catalogue=_catalogue_for_image(cfg, image_path),
    )


def _configured_products(cfg: Config, role: str, target: Target) -> list[ImageProduct]:
    """Resolve full-band inputs conservatively when more than one is configured."""

    images = list(dict.fromkeys(target.images))
    if len(images) == 1:
        bands = ["mfs"]
    else:
        # Multiple arbitrary images cannot safely be identified as one MFS product.
        # Pair them by configuration order, but exclude them from the three-band
        # flux handoff, which specifically requires a unique ``mfs`` product.
        bands = [f"image_{index}" for index in range(1, len(images) + 1)]
    return [
        _product(
            cfg,
            role=role,
            target_name=target.name or role,
            band=band,
            image=image,
        )
        for band, image in zip(bands, images)
    ]


def _image_products(cfg: Config) -> list[ImageProduct]:
    products = _configured_products(cfg, "reference", cfg.reference)
    for target in cfg.tests:
        products.extend(_configured_products(cfg, "test", target))

    if cfg.low_high_slice.enabled and cfg.low_high_slice.add_to_source_finding:
        for target in resolved_targets(cfg):
            products.extend(
                (
                    _product(
                        cfg,
                        role=target.role,
                        target_name=target.name,
                        band="low",
                        image=str(target.low_image),
                    ),
                    _product(
                        cfg,
                        role=target.role,
                        target_name=target.name,
                        band="high",
                        image=str(target.high_image),
                    ),
                )
            )

    unique: dict[tuple[str, str, str], ImageProduct] = {}
    for product in products:
        unique.setdefault((product.role, product.target_name, product.band), product)
    return list(unique.values())


def _pair_key(input1: str, input2: str) -> tuple[str, str]:
    return _absolute(input1), _absolute(input2)


def _canonical_xmatch_output(cfg: Config, entry: Iterable[str]) -> str:
    values = list(entry)
    input1, input2 = values[:2]
    if len(values) >= 3:
        requested = Path(values[2]).expanduser().resolve(strict=False)
    else:
        stem = f"{Path(input1).stem}_X_{Path(input2).stem}.fits"
        requested = Path(cfg.paths.sky_xmatches_dir).expanduser().resolve(strict=False) / stem

    if requested.parent.name != "Sky-CrossMatches":
        requested = requested.parent / "Sky-CrossMatches" / requested.name
    return str(requested)


def _configured_flux_analyses(raw: Any) -> list[dict[str, Any]]:
    if not raw:
        return []
    if isinstance(raw, dict):
        return [dict(raw)]
    if isinstance(raw, list):
        return [dict(entry) for entry in raw]
    raise ValueError("config.extra.flux must be a mapping or a list of mappings")


def build_image_pipeline_plan(cfg: Config) -> ImagePipelinePlan:
    """Build deterministic downstream inputs for the image-only pipeline.

    Explicit ``extra.xmatch_pairs``, ``extra.positions`` and ``extra.flux``
    entries are retained. Generated entries fill only reference/test product
    pairs that are not already wired explicitly.
    """

    products = _image_products(cfg)
    references = {product.band: product for product in products if product.role == "reference"}

    tests_by_name: dict[str, dict[str, ImageProduct]] = {}
    for product in products:
        if product.role != "test":
            continue
        tests_by_name.setdefault(product.target_name, {})[product.band] = product

    configured_pairs = [list(entry) for entry in (cfg.extra.get("xmatch_pairs") or [])]
    configured_pair_outputs: dict[tuple[str, str], str] = {}
    for entry in configured_pairs:
        if len(entry) >= 2:
            configured_pair_outputs[_pair_key(entry[0], entry[1])] = _canonical_xmatch_output(
                cfg, entry
            )

    inferred_pairs: list[list[str]] = []
    matches: list[XMatchProduct] = []
    for target in cfg.tests:
        test_products = tests_by_name.get(target.name, {})
        for band, reference in references.items():
            test = test_products.get(band)
            if test is None:
                continue
            key = _pair_key(reference.catalogue, test.catalogue)
            output = configured_pair_outputs.get(key)
            if output is None:
                output_name = (
                    f"{_slug(cfg.reference.name)}_X_{_slug(target.name)}_{_slug(band)}.fits"
                )
                output = _absolute(Path(cfg.paths.sky_xmatches_dir) / output_name)
                inferred_pairs.append([reference.catalogue, test.catalogue, output])
            matches.append(
                XMatchProduct(
                    target_name=target.name,
                    band=band,
                    reference=reference,
                    test=test,
                    table=output,
                )
            )

    all_pairs = configured_pairs + inferred_pairs

    configured_positions = [
        dict(entry) for entry in (cfg.extra.get("positions") or [])
    ]
    positioned_tables = {
        _absolute(entry["xmatch_table"])
        for entry in configured_positions
        if entry.get("xmatch_table")
    }
    inferred_positions: list[dict[str, Any]] = []
    for match in matches:
        if _absolute(match.table) in positioned_tables:
            continue
        inferred_positions.append(
            {
                "xmatch_table": match.table,
                "ref_fits": match.reference.image,
                "other_fits": match.test.image,
                "otherdatatag": f"{match.target_name}_{match.band}",
            }
        )

    configured_flux = _configured_flux_analyses(cfg.extra.get("flux"))
    configured_flux_keys = {
        (
            _absolute(entry["ref_low_xmatch"]),
            _absolute(entry["ref_high_xmatch"]),
            _absolute(entry["ref_mfs_xmatch"]),
        )
        for entry in configured_flux
        if all(
            key in entry
            for key in ("ref_low_xmatch", "ref_high_xmatch", "ref_mfs_xmatch")
        )
    }
    matches_by_test: dict[str, dict[str, XMatchProduct]] = {}
    for match in matches:
        matches_by_test.setdefault(match.target_name, {})[match.band] = match

    inferred_flux: list[dict[str, Any]] = []
    for target in cfg.tests:
        by_band = matches_by_test.get(target.name, {})
        if not all(band in by_band for band in ("low", "high", "mfs")):
            continue
        analysis = {
            "ref_low_xmatch": by_band["low"].table,
            "ref_high_xmatch": by_band["high"].table,
            "ref_mfs_xmatch": by_band["mfs"].table,
        }
        key = tuple(_absolute(analysis[name]) for name in analysis)
        if key not in configured_flux_keys:
            inferred_flux.append(analysis)

    return ImagePipelinePlan(
        products=tuple(products),
        matches=tuple(matches),
        xmatch_pairs=tuple(tuple(str(value) for value in entry[:3]) for entry in all_pairs),
        positions=tuple(configured_positions + inferred_positions),
        flux=tuple(configured_flux + inferred_flux),
    )


def wire_image_pipeline(cfg: Config) -> ImagePipelinePlan:
    """Add the planned downstream handoffs to the in-memory configuration."""

    plan = build_image_pipeline_plan(cfg)
    cfg.extra = dict(cfg.extra)
    cfg.extra["xmatch_pairs"] = [list(entry) for entry in plan.xmatch_pairs]
    cfg.extra["positions"] = [dict(entry) for entry in plan.positions]
    cfg.extra["flux"] = [dict(entry) for entry in plan.flux]
    return plan
