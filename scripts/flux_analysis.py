#!/usr/bin/env python3
"""Compatibility entry point for the packaged flux analysis."""

from meerkat_corr_imaging.flux_analysis import main


if __name__ == "__main__":
    raise SystemExit(main())
