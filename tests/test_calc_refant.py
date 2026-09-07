import numpy as np
import pytest
from meerkat_corr_imaging.calc_refant import get_ref_ant


class Metadata:
    def open(self, path): pass
    def fieldsforname(self, name): return [0] if name == 'flux' else []
    def nfields(self): return 1
    def antennanames(self, ant): return ['m%03d' % ant]
    def done(self): self.closed = True


class Rows:
    def __init__(self, rows): self.rows = rows
    def getcol(self, name): return np.array([row[name] for row in self.rows])
    def nrows(self): return len(self.rows)
    def getcell(self, name, row): return self.rows[row][name]
    def close(self): pass


class Table:
    def __init__(self, rows): self.rows = rows
    def open(self, path): pass
    def query(self, query):
        assert 'ANTENNA1!=ANTENNA2' in query
        rows = self.rows
        if ' OR ' in query:
            import re
            ant = int(re.search(r'ANTENNA1==(\d+)', query).group(1))
            assert 'ANTENNA2==%d' % ant in query
            rows = [r for r in rows if ant in (r['ANTENNA1'], r['ANTENNA2'])]
        return Rows(rows)
    def done(self): self.closed = True


def row(a, b, flags, row_flag=False):
    return dict(ANTENNA1=a, ANTENNA2=b, FLAG=np.array(flags), FLAG_ROW=row_flag)


def test_actual_id_both_baseline_ends_and_row_flags(tmp_path):
    import json
    md = Metadata()
    tb = Table([row(2, 7, [False], True), row(7, 11, [False, False])])
    output = tmp_path / 'stats.json'
    refant, bad = get_ref_ant('test.ms', 'flux', str(output), lambda: md, lambda: tb)
    assert refant == 'm011'  # only appears as ANTENNA2; not list index 2
    assert bad == [2]
    assert json.loads(output.read_text())['refant'] == refant
    assert md.closed and tb.closed


def test_all_flagged_fails_and_closes_tools():
    md = Metadata()
    tb = Table([row(2, 7, [True])])
    with pytest.raises(ValueError, match='No unflagged'):
        get_ref_ant('test.ms', '0', msmd_factory=lambda: md, table_factory=lambda: tb)
    assert md.closed and tb.closed


def test_missing_field_fails_and_closes_tools():
    md, tb = Metadata(), Table([])
    with pytest.raises(ValueError, match='Flux field not found'):
        get_ref_ant('test.ms', 'missing', msmd_factory=lambda: md, table_factory=lambda: tb)
    assert md.closed and tb.closed


def test_tie_uses_lowest_antenna_id():
    assert get_ref_ant('test.ms', 'flux', msmd_factory=Metadata,
                       table_factory=lambda: Table([row(7, 11, [False])])) == ('m007', [])
