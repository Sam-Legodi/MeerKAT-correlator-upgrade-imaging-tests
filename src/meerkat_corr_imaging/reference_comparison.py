"""Statistical comparisons with a paired CMC1 reference observation."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2

from .uncertainty import UncertaintyConfig, interval


COMPARISON_CONFIDENCE = 0.95


def comparison_config(config: UncertaintyConfig | None = None) -> UncertaintyConfig:
    config = config or UncertaintyConfig()
    return UncertaintyConfig(config.bootstrap_samples, COMPARISON_CONFIDENCE, config.random_seed)


def parity_decision(info: dict | None, reference: float = 1.0) -> dict:
    """Two-sided concern when the 95% interval excludes reference parity."""
    low, high = (info or {}).get("ci_low"), (info or {}).get("ci_high")
    if low is None or high is None or not np.isfinite([low, high]).all():
        return {"status": "Not assessed", "reason": "95% comparison interval unavailable"}
    return {"status": "Concern" if low > reference or high < reference else "Pass",
            "ci_low": float(low), "ci_high": float(high), "reference": reference,
            "confidence_level": COMPARISON_CONFIDENCE}


def degradation_decision(info: dict | None, reference: float = 0.0) -> dict:
    """One-sided concern only for a supported increase above CMC1."""
    low, high = (info or {}).get("ci_low"), (info or {}).get("ci_high")
    if low is None or high is None or not np.isfinite([low, high]).all():
        return {"status": "Not assessed", "reason": "95% comparison interval unavailable"}
    return {"status": "Concern" if low > reference else "Pass",
            "ci_low": float(low), "ci_high": float(high), "reference": reference,
            "confidence_level": COMPARISON_CONFIDENCE}


def position_reference_test(separation_p95: float | None, covariance: np.ndarray,
                            translation: np.ndarray | None, translation_covariance: np.ndarray | None,
                            config: UncertaintyConfig | None = None) -> dict:
    """Compare radial p95 with catalogue-noise null and test rigid translation."""
    cfg = comparison_config(config)
    cov = np.asarray(covariance, float)
    valid = np.isfinite(cov).all(axis=(1, 2)) & (np.linalg.eigvalsh(np.nan_to_num(cov)).min(axis=1) >= -1e-12)
    cov = cov[valid]
    if separation_p95 is None or not np.isfinite(separation_p95) or len(cov) < 3:
        return {"status": "Not assessed", "reason": "insufficient finite source position covariance"}
    rng = np.random.default_rng(cfg.random_seed)
    eigenvalue, eigenvector = np.linalg.eigh(cov)
    root = eigenvector * np.sqrt(np.maximum(eigenvalue, 0))[:, None, :]
    draws = []
    for start in range(0, cfg.bootstrap_samples, 128):
        n = min(128, cfg.bootstrap_samples - start)
        offsets = np.einsum("nij,bnj->bni", root, rng.normal(size=(n, len(cov), 2)))
        draws.extend(np.percentile(np.linalg.norm(offsets, axis=2), 95, axis=1))
    null = interval(draws, float(np.median(draws)), cfg, method="catalogue covariance noise-only Monte Carlo", n=len(cov))
    chi_square = None
    translation_concern = False
    if translation is not None and translation_covariance is not None:
        vector = np.asarray(translation, float)
        matrix = np.asarray(translation_covariance, float)
        if vector.shape == (2,) and matrix.shape == (2, 2) and np.isfinite(vector).all() and np.isfinite(matrix).all() and np.linalg.eigvalsh(matrix).min() > 0:
            chi_square = float(vector @ np.linalg.solve(matrix, vector))
            translation_concern = chi_square > float(chi2.ppf(cfg.confidence_level, 2))
    p95_concern = separation_p95 > null["ci_high"]
    return {"status": "Concern" if p95_concern or translation_concern else "Pass",
            "observed_p95_arcsec": float(separation_p95), "noise_only_p95": null,
            "translation_chi2": chi_square, "translation_chi2_limit": float(chi2.ppf(cfg.confidence_level, 2)),
            "p95_concern": bool(p95_concern), "translation_concern": bool(translation_concern),
            "confidence_level": cfg.confidence_level}


def visibility_exposure(directory: str | Path, polarisation: str = "XX") -> dict:
    """Effective exposure proxy: unflagged cross samples × channel width × dump time."""
    directory = Path(directory)
    flags = pd.read_csv(directory / "flagging_vs_scan.csv")
    selected = flags[(flags.CLASS == "cross") & (flags.POL == polarisation)]
    unflagged = float((selected.TOTAL - selected.FLAGGED).sum())
    frequency = np.sort(pd.read_csv(directory / "rfi_free_channel_mask.csv", usecols=["FREQ_HZ"]).FREQ_HZ.dropna().unique())
    time = np.sort(pd.read_csv(directory / "perrow_amp_stats.csv", usecols=["TIME"]).TIME.dropna().unique())
    frequency_step = np.diff(frequency)
    time_step = np.diff(time)
    frequency_step = frequency_step[frequency_step > 0]
    time_step = time_step[time_step > 0]
    if not unflagged > 0 or not len(frequency_step) or not len(time_step):
        raise ValueError("effective exposure requires cross samples, frequency spacing and dump spacing")
    width = float(np.median(frequency_step))
    dump = float(np.median(time_step))
    return {"polarisation_proxy": polarisation,
            "unflagged_cross_samples": unflagged, "channel_width_hz": width,
            "dump_seconds": dump, "effective_exposure": unflagged * width * dump,
            "assumption": "XX exposure represents the imaged parallel hands; equal effective visibility weights and comparable imaging weighting"}
