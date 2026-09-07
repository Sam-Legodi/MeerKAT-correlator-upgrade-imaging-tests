"""Completion records for reference visibility QA, scoped to its output directory."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

MARKER = 'reference_qa_complete.json'


def _identity(ms):
    source = Path(ms).resolve()
    if not source.exists():
        return None
    # Include table data and metadata, including MMS sub-tables, without reading
    # visibility arrays. Casacore lock files change simply from opening a table.
    files = sorted(p for p in source.rglob('*') if p.is_file() and p.name != 'table.lock')
    signature = [(str(p.relative_to(source)), p.stat().st_size, p.stat().st_mtime_ns)
                 for p in files]
    analyser = Path(__file__).with_name('vis_amp_analyze.py').read_bytes()
    return {'ms': str(source), 'files': signature,
            'analyser': hashlib.sha256(analyser).hexdigest()}


def reference_is_complete(ms, output):
    output = Path(output)
    try:
        record = json.loads((output / MARKER).read_text())
        identity = _identity(ms)
        # JSON normalizes signature tuples to lists.
        if identity is None or record['identity'] != json.loads(json.dumps(identity)):
            return False
        if not record['outputs']:
            return False
        return all((output / name).is_file() and (output / name).stat().st_size == size
                   and size > 0 for name, size in record['outputs'].items())
    except (OSError, ValueError, KeyError, TypeError):
        return False


def record_reference_completion(ms, output):
    output = Path(output)
    identity = _identity(ms)
    required = ['scan_averaged_amp_stats.csv', 'acceptance_summary.csv',
                'amplitude_by_baseline.csv', 'flagging_by_channel.csv']
    reports = list(output.glob('draft_vis_amp_summary*.docx'))
    if identity is None or not reports or not all((output / name).is_file() for name in required):
        return
    outputs = {p.name: p.stat().st_size for p in output.iterdir()
               if p.suffix in {'.csv', '.png', '.docx'} and p.is_file()}
    (output / MARKER).write_text(json.dumps({'identity': identity, 'outputs': outputs}, indent=2))
