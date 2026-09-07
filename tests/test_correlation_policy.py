import ast
from pathlib import Path
import pytest
from meerkat_corr_imaging.correlation_policy import parallel_linear_only, ms_parallel_linear_only


@pytest.mark.parametrize('rows,expected', [([[9]], True), ([[12]], True),
    ([[9,12]], True), ([[9], [12]], True), ([[9,10,11,12]], False),
    ([[5,8]], False), ([[9,12], [9,10,11,12]], False)])
def test_policy(rows, expected):
    assert parallel_linear_only(rows) is expected


@pytest.mark.parametrize('rows', [[], [[]]])
def test_missing_metadata(rows):
    with pytest.raises(ValueError):
        parallel_linear_only(rows)


def test_table_cleanup():
    class Table:
        def open(self, path): assert path.endswith('/POLARIZATION')
        def nrows(self): return 1
        def getcell(self, column, row):
            assert column == 'CORR_TYPE'
            return [9, 12]
        def done(self): self.closed = True
    table = Table()
    assert ms_parallel_linear_only('input.ms', lambda: table)
    assert table.closed


@pytest.mark.parametrize('parallel', [True, False])
def test_solver_task_policy(parallel):
    import meerkat_corr_imaging.correlation_policy as module
    source = Path(module.__file__).with_name('standalone_xxyy_solve.py').read_text()
    # Exercise the actual task-call section with recording CASA stubs.
    source = source.split('# setjy (exact behavior from your setjy.py)', 1)[1]
    tree = ast.parse(source)
    names = {n.id for n in ast.walk(tree) if isinstance(n, ast.Name)}
    ns = dict.fromkeys(names, 'value')
    calls = []
    for task in ('gaincal', 'bandpass', 'applycal', 'polcal', 'fluxscale'):
        ns[task] = lambda _task=task, **kw: calls.append((_task, kw))
    ns.update(len=len, apply_parang=not parallel, do_leakage=not parallel, do_split=False,
              interp_all=['linear']*4, _log=lambda *a: None,
              do_setjy=lambda **kw: None, Fields=lambda x: x,
              out=lambda *a: 'output', datetime=__import__('datetime').datetime)
    exec(compile(tree, 'solver-tasks', 'exec'), ns)
    for task, kw in calls:
        if task in ('gaincal', 'bandpass', 'applycal'):
            assert kw['parang'] is (False if kw.get('gaintype') == 'K' else not parallel)
    assert any(task == 'polcal' for task, _ in calls) is (not parallel)
    apply = next(kw for task, kw in calls if task == 'applycal')
    assert len(apply['gaintable']) == (4 if parallel else 5)
