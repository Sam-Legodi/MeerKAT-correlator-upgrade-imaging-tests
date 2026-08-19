#!/usr/bin/env python
"""Compatibility entry point for the packaged CASA calibration script."""

from __future__ import print_function

import os

import meerkat_corr_imaging


_module_path = os.path.join(
    os.path.dirname(meerkat_corr_imaging.__file__),
    "standalone_xxyy_solve.py",
)
with open(_module_path, "rb") as _handle:
    _code = compile(_handle.read(), _module_path, "exec")
globals()["__file__"] = _module_path
exec(_code, globals(), globals())
