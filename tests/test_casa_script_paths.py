"""CASA 5 executes scripts without the usual Python module globals."""
import ast
import os
from pathlib import Path
from types import SimpleNamespace

import pytest
from meerkat_corr_imaging import config


@pytest.mark.parametrize('script', ['standalone_xxyy_solve.py', 'tclean_two_bands.py'])
@pytest.mark.parametrize('mode', ['environment', 'argv', 'file'])
def test_helpers_without_casa_imports(script, mode, monkeypatch, tmp_path):
    package = Path(config.__file__).parent
    tree = ast.parse((package / script).read_text())
    function = next(node for node in tree.body
                    if isinstance(node, ast.FunctionDef) and node.name == '_casa_helper_path')
    namespace = {'os': os, 'sys': SimpleNamespace(argv=['casa'])}
    monkeypatch.delenv('MCI_CASA_SCRIPT_DIR', raising=False)
    monkeypatch.chdir(tmp_path)
    if mode == 'environment':
        monkeypatch.setenv('MCI_CASA_SCRIPT_DIR', str(package))
    elif mode == 'argv':
        namespace['sys'].argv += ['--nologger', '-c', str(package / script)]
    else:
        namespace['__file__'] = str(package / script)
    exec(compile(ast.Module(body=[function], type_ignores=[]), script, 'exec'), namespace)
    assert namespace['_casa_helper_path']('field_selection.py') == str(package / 'field_selection.py')
    with pytest.raises(IOError, match='CASA helper not found'):
        namespace['_casa_helper_path']('missing.py')
