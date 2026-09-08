"""Check the CASA model against the supplied MeerCals polynomial."""
import ast
import logging
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from meerkat_corr_imaging import config


@pytest.fixture
def model():
    path = Path(config.__file__).with_name('standalone_xxyy_solve.py')
    tree = ast.parse(path.read_text())
    functions = [node for node in tree.body if isinstance(node, ast.FunctionDef)]
    namespace = {'np': np, 'logger': logging.getLogger(__name__)}
    exec(compile(ast.Module(body=functions, type_ignores=[]), str(path), 'exec'), namespace)
    return namespace


class Metadata:
    def __init__(self, frequencies):
        self.frequencies = frequencies

    def nspw(self):
        return len(self.frequencies)

    def chanfreqs(self, spw, unit):
        assert unit == 'Hz'
        return self.frequencies[spw]

    def open(self, vis):
        pass

    def namesforfields(self, field):
        return ['J0408-6545'] if field == 2 else ['other']

    def fieldnames(self):
        return ['other', 'J0408-6545']

    def done(self):
        pass


@pytest.mark.parametrize('limits', [(0.5, 1.01), (0.856, 1.70), (2.0, 3.5)])
@pytest.mark.parametrize('count', [1, 2, 200])
def test_spectrum_matches_supplied_model(model, limits, count):
    frequencies = np.linspace(*limits, count) * 1e9
    reference, flux, s0, s1, s2 = model['j0408_flux_model_from_msmd'](Metadata([frequencies]))
    # Evaluate across the entire band, including when metadata has only one channel.
    nu = np.linspace(*limits, 200) * 1e9
    x = np.log10(nu / 1e6)
    expected = 10 ** (-0.9790 + 3.3662*x - 1.1216*x**2 + 0.0861*x**3)
    ratio = nu / reference
    log_ratio = np.log10(ratio)
    actual = flux * ratio ** (s0 + s1*log_ratio + s2*log_ratio**2)
    np.testing.assert_allclose(actual, expected, rtol=1e-12)


@pytest.mark.parametrize('field,manual', [('J0408-6545', True), ('0408-6545', True),
                                        ('2', True), (' 2 ', True), ('other', False)])
def test_setjy_dispatch(model, field, manual):
    calls = []
    model.update(msmd=Metadata([[0.6e9], [2e9, 3e9]]),
                 delmod=lambda **kwargs: None, setjy=lambda **kwargs: calls.append(kwargs))
    model['do_setjy']('test.ms', '1', SimpleNamespace(fluxfield=field), 'Perley-Butler 2017')
    call, = calls
    assert call['field'] == field.strip()
    assert call['spw'] == '1'
    assert call['scalebychan'] is True
    assert call['standard'] == ('manual' if manual else 'Perley-Butler 2017')
    if manual:
        assert call['reffreq'] == '2500000000.000000 Hz'
        assert call['fluxdensity'][1:] == [0.0, 0.0, 0.0]
        assert call['spix'][2:] == [0.0861, 0.0]


def test_reject_nonpositive_frequency(model):
    with pytest.raises(ValueError, match='positive'):
        model['j0408_flux_model_from_msmd'](Metadata([[0, 1e9]]))
