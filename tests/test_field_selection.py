import pytest
from meerkat_corr_imaging.field_selection import select_fields
from meerkat_corr_imaging.config import _dict_to_dataclass
from meerkat_corr_imaging.cli import _normalize_cli_args


def test_all_fields_and_exclusions():
    names = ['flux', 'gain', 'target']
    assert select_fields(names, []) == '0,1,2'
    assert select_fields(names, ['gain', 2]) == '0'
    assert select_fields(names, ['gain'], ['gain', 'target']) == '2'
    assert select_fields(names, ['target'], ['target'], True) is None
    with pytest.raises(ValueError, match='every'):
        select_fields(names, names)
    with pytest.raises(ValueError, match='Unknown'):
        select_fields(names, ['typo'])


def test_config_controls():
    cfg = _dict_to_dataclass({'project_name': 'test', 'casa': {
        'imaging_enabled': False, 'exclude_fields': [1, 'foo']}})
    assert not cfg.casa.imaging_enabled
    assert cfg.casa.exclude_fields == ['1', 'foo']
    with pytest.raises(ValueError):
        _dict_to_dataclass({'project_name': 'test', 'casa': {'exclude_fields': 'foo'}})


def test_cli_flag_normalization():
    assert _normalize_cli_args(['cal', '--no-imaging', '--config', 'x']) == [
        '--config', 'x', '--no-imaging', 'cal']
