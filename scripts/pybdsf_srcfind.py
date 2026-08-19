#!/usr/bin/env python3
"""Compatibility entry point for the packaged PyBDSF launcher."""

from meerkat_corr_imaging.pybdsf_srcfind import main


if __name__ == "__main__":
    raise SystemExit(main())
