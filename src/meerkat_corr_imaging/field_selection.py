"""CASA 5/6 compatible exact field-name/ID exclusion helpers."""
from __future__ import print_function


def select_fields(names, excluded, included=None, allow_empty=False):
    """Return explicit IDs; never return CASA's empty-string/all-fields sentinel."""
    def resolve(tokens):
        ids = set()
        for token in tokens:
            token = str(token).strip()
            matches = [i for i, name in enumerate(names) if name == token]
            if token.isdigit() and int(token) < len(names):
                matches = [int(token)]
            if not matches:
                raise ValueError("Unknown field name or ID: " + token)
            ids.update(matches)
        return ids
    selected = set(range(len(names))) if included is None else resolve(included)
    selected -= resolve(excluded)
    if not selected and not allow_empty:
        raise ValueError("Field selection excludes every selected field")
    return ','.join(str(i) for i in sorted(selected)) if selected else None


def ms_field_selection(vis, excluded, included=None, allow_empty=False):
    try:
        from casatools import msmetadata
    except ImportError:
        from taskinit import msmdtool as msmetadata
    tool = msmetadata()
    try:
        tool.open(vis)
        return select_fields(list(tool.fieldnames()), excluded, included, allow_empty)
    finally:
        tool.done()
