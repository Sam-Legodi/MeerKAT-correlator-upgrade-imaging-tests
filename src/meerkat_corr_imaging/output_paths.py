"""Naming helpers for pipeline output artifacts."""

from __future__ import annotations

import os
from os import PathLike


def draft_docx_path(path: str | PathLike[str]) -> str:
    """Return *path* with ``draft_`` prefixed to its DOCX basename once."""
    path_string = os.fspath(path)
    directory, basename = os.path.split(path_string)
    if basename.startswith("draft_"):
        return path_string
    return os.path.join(directory, f"draft_{basename}")
