# Copyright (C) 2022 Inter-University Institute for Data Intensive Astronomy
# Adapted from legacy_scripts/calc_refant.py; see processMeerKAT.py for license details.
"""Reference antenna selection, importable in CASA 5 (Python 2) and CASA 6."""
from __future__ import print_function

import json
import numpy as np


def get_ref_ant(visname, fluxfield, stats_path=None, msmd_factory=None, table_factory=None):
    """Rank antennas by flagged fraction on flux-calibrator cross-correlations.

    Include both baseline ends and row flags. Ties select the lowest antenna ID.
    Return the CASA antenna name and IDs with more than 80 percent flagged data.
    Tool factories allow reuse with CASA 5 and testing without CASA installed.
    """
    if msmd_factory is None or table_factory is None:
        try:
            from casatools import msmetadata, table
        except ImportError:
            from taskinit import msmdtool as msmetadata, tbtool as table
        msmd_factory = msmd_factory or msmetadata
        table_factory = table_factory or table
    msmd = msmd_factory()
    tb = table_factory()
    try:
        msmd.open(visname)
        field_ids = []
        for selector in str(fluxfield).split(','):
            selector = selector.strip()
            if selector.isdigit():
                ids = [int(selector)]
                if ids[0] >= msmd.nfields():
                    raise ValueError("Unknown flux field ID: " + selector)
            else:
                ids = list(msmd.fieldsforname(selector)) if selector else []
            if not ids:
                raise ValueError("Flux field not found: " + selector)
            field_ids.extend(int(value) for value in ids)
        field_ids = sorted(set(field_ids))
        tb.open(visname)
        field_query = 'FIELD_ID IN [{}] AND ANTENNA1!=ANTENNA2'.format(
            ','.join(str(value) for value in field_ids))
        subset = tb.query(field_query)
        try:
            antennas = sorted(set(int(a) for a in subset.getcol('ANTENNA1')) |
                              set(int(a) for a in subset.getcol('ANTENNA2')))
        finally:
            subset.close()
        stats = []
        for ant in antennas:
            rows = tb.query('({}) AND (ANTENNA1=={} OR ANTENNA2=={})'.format(
                field_query, ant, ant))
            total = flagged = 0
            try:
                # Row-wise reads support variable channel counts between SPWs.
                for row in range(rows.nrows()):
                    flags = np.asarray(rows.getcell('FLAG', row), dtype=bool)
                    total += flags.size
                    flagged += flags.size if rows.getcell('FLAG_ROW', row) else np.count_nonzero(flags)
            finally:
                rows.close()
            fraction = float(flagged) / total if total else 1.0
            stats.append({'antenna_id': ant, 'flagged_fraction': fraction,
                          'samples': int(total)})
        usable = [entry for entry in stats if entry['samples'] and entry['flagged_fraction'] < 1.0]
        if not usable:
            raise ValueError("No unflagged flux-calibrator cross-correlations in " + visname)
        best = min(usable, key=lambda entry: (entry['flagged_fraction'], entry['antenna_id']))
        refant = str(msmd.antennanames(best['antenna_id'])[0])
        badants = [entry['antenna_id'] for entry in stats if entry['flagged_fraction'] > 0.8]
        if stats_path:
            with open(stats_path, 'w') as handle:
                json.dump({'ms': visname, 'field_ids': field_ids, 'refant': refant,
                           'badants': badants, 'antennas': stats}, handle, indent=2)
        return refant, badants
    finally:
        # Release both tools even when metadata lookup or table reads fail.
        try:
            tb.done()
        finally:
            msmd.done()
