import numpy as np
import pytest
from meerkat_corr_imaging import calibration_checks


class Rows:
    def __init__(self, values, flagged=False):
        self.values = np.asarray(values)
        self.flagged = flagged
        self.closed = False
    def nrows(self): return 1
    def query(self, *args):
        pytest.fail('Nested queries break CASA 5 temporary tables')
    def getcell(self, column, row):
        return {'CORRECTED_DATA': self.values,
                'FLAG': np.zeros(self.values.shape, dtype=bool),
                'FLAG_ROW': self.flagged}[column]
    def close(self): self.closed = True


class Table:
    def __init__(self, rows):
        self.rows = rows
        self.queries = []
        self.closed = False
    def open(self, vis): assert vis == 'input.ms'
    def colnames(self): return ['CORRECTED_DATA']
    def query(self, expression):
        self.queries.append(expression)
        assert 'AND ANTENNA1!=ANTENNA2' in expression
        return self.rows[len(self.queries)-1]
    def done(self): self.closed = True


def test_each_field_queries_original_ms(monkeypatch):
    rows = [Rows([[1+2j, 3]]), Rows([[4, 5, 6]])]
    table = Table(rows)
    monkeypatch.setattr(calibration_checks, 'table_tool', lambda: table)
    assert calibration_checks.inspect_corrected('input.ms', '0,2') == {'0': 2, '2': 3}
    assert table.queries == ['FIELD_ID==0 AND ANTENNA1!=ANTENNA2',
                             'FIELD_ID==2 AND ANTENNA1!=ANTENNA2']
    assert table.closed and all(row.closed for row in rows)


@pytest.mark.parametrize('values,flagged,message', [([np.nan], False, 'Nonfinite'),
                                                  ([1], True, 'No usable')])
def test_failure_still_closes_tables(monkeypatch, values, flagged, message):
    rows = Rows(values, flagged)
    table = Table([rows])
    monkeypatch.setattr(calibration_checks, 'table_tool', lambda: table)
    with pytest.raises(RuntimeError, match=message):
        calibration_checks.inspect_corrected('input.ms', '0')
    assert table.closed and rows.closed
