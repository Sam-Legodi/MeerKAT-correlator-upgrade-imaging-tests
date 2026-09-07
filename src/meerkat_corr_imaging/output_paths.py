"""Naming, registration and PDF export for pipeline reports."""
from __future__ import annotations

import json
import os
from pathlib import Path
import re


def draft_docx_path(path, observation=None) -> str:
    """Prefix drafts once and append an observation label when available."""
    path = Path(path)
    if observation is None:
        observation = next((p.name for p in path.parents
                            if re.match(r"^\d{10}(?:_|$)", p.name)), None)
    observation = observation or os.environ.get("MCI_REPORT_OBSERVATION")
    stem = path.stem if path.stem.startswith("draft_") else "draft_" + path.stem
    if observation:
        label = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(observation))
        if not stem.endswith("_" + label):
            stem += "_" + label
    return str(path.with_name(stem + path.suffix))


def save_report(document, path):
    """Save a DOCX and register it across pipeline subprocesses."""
    path = Path(path).resolve()
    document.save(path)
    manifest = os.environ.get("MCI_REPORT_MANIFEST")
    if manifest:
        with open(manifest, "a", encoding="utf-8") as stream:
            stream.write(json.dumps(str(path)) + "\n")
    else:
        export_pdf(path)


def export_pdf(path):
    """Generate a PDF using the pipeline's existing Python dependencies."""
    from .report_pdf import render_docx_pdf

    return render_docx_pdf(Path(path))
