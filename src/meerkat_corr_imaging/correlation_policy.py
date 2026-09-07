"""CASA 5/6 compatible calibration policy for linear parallel-hand inputs."""
from __future__ import print_function
import os


def parallel_linear_only(correlation_rows):
    """MS Stokes enumeration: XX=9, YY=12; reject absent metadata."""
    rows = [list(row) for row in correlation_rows]
    if not rows or any(not row for row in rows):
        raise ValueError("Missing CORR_TYPE metadata; cannot choose calibration policy")
    return all(int(corr) in (9, 12) for row in rows for corr in row)


def ms_parallel_linear_only(vis, table_factory=None):
    if table_factory is None:
        try:
            from casatools import table as table_factory
        except ImportError:
            from taskinit import tbtool as table_factory
    tool = table_factory()
    try:
        tool.open(os.path.join(vis, 'POLARIZATION'))
        return parallel_linear_only([tool.getcell('CORR_TYPE', row)
                                     for row in range(tool.nrows())])
    finally:
        tool.done()
