"""Pipeline wrapper for the consolidated imaging verification report."""

from __future__ import annotations

from ..config import Config
from ..verification_report import build_verification_report


def run(cfg: Config):
    return build_verification_report(cfg)

