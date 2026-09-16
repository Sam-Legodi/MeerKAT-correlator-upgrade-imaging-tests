#!/usr/bin/env python3
"""Read-only preflight. Does not run pipeline stages or create reports/directories."""
import argparse
from pathlib import Path
import sys
import yaml
from meerkat_corr_imaging.config import _dict_to_dataclass

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('kind', choices=['imaging', 'visibility'])
p.add_argument('configs', nargs='+')
a = p.parse_args()
failures = []
for filename in a.configs:
    try:
        cfg = _dict_to_dataclass(yaml.safe_load(Path(filename).read_text()))
        if a.kind == 'imaging':
            from astropy.io import fits
            from astropy.table import Table
            for image in cfg.reference.images + [v for t in cfg.tests for v in t.images]:
                if not Path(image).is_file():
                    raise ValueError('missing image: ' + image)
                h = fits.getheader(image)
                if not all(k in h for k in ['BMAJ', 'BMIN', 'CTYPE1', 'CTYPE2']):
                    raise ValueError('beam/WCS metadata missing: ' + image)
            for pair in cfg.extra['xmatch_pairs']:
                for cat in pair[:2]:
                    t = Table.read(cat)
                    required = {'RA','DEC','Peak_flux','Total_flux','E_Peak_flux'}
                    if not required.issubset(t.colnames):
                        raise ValueError(f'missing required columns {required-set(t.colnames)}: {cat}')
                    optional = {'E_RA','E_DEC','E_Total_flux'}-set(t.colnames)
                    if optional:
                        print(f'WARNING measurement errors unavailable {optional}: {cat}')
        else:
            meta = cfg.extra['refresh_metadata']
            ms_paths = cfg.reference.ms_paths + [ms for t in cfg.tests for ms in t.ms_paths]
            if any('REPLACE' in ms for ms in ms_paths):
                raise ValueError('replace reference_ms/test_ms placeholders first')
            if any('REPLACE' in meta[k] for k in ['corrected_data_provenance','spectral_setup']):
                raise ValueError('record SDP correction provenance and spectral setup in refresh_metadata')
            from casacore.tables import table
            import numpy as np
            for ms in ms_paths:
                if not Path(ms).is_dir():
                    raise ValueError('missing MS: '+ms)
                with table(ms+'::FIELD', readonly=True, ack=False) as ft:
                    names = list(ft.getcol('NAME'))
                ids = [i for i,name in enumerate(names) if name == meta['expected_field']]
                if not ids:
                    raise ValueError(f"expected {meta['expected_field']}; FIELD names are {names}: {ms}")
                with table(ms, readonly=True, ack=False) as mt:
                    if not mt.nrows(): raise ValueError('empty MS: '+ms)
                    for start in range(0, mt.nrows(), 100000):
                        used = mt.getcol('FIELD_ID', startrow=start, nrow=min(100000,mt.nrows()-start))
                        if not np.isin(used,ids).all():
                            raise ValueError('contains other fields; wrapper analyzes ALL rows: '+ms)
                    chosen = next((c for c in ['CORRECTED_DATA','DATA','MODEL_DATA'] if c in mt.colnames()),None)
                    if chosen not in ['CORRECTED_DATA','DATA'] or not mt.iscelldefined(chosen,0):
                        raise ValueError('no populated observed data column: '+ms)
                    required = {'FLAG','TIME','SCAN_NUMBER','ANTENNA1','ANTENNA2','DATA_DESC_ID'}
                    if not required.issubset(mt.colnames()): raise ValueError('required MS columns missing: '+ms)
                    print(f'  {ms}: {mt.nrows()} rows, field={meta["expected_field"]}, selected column={chosen}')
                with table(ms+'::SPECTRAL_WINDOW',readonly=True,ack=False) as spw:
                    for i in range(spw.nrows()):
                        freq=spw.getcell('CHAN_FREQ',i);width=spw.getcell('CHAN_WIDTH',i)
                        print(f'    SPW {i}: {len(freq)} channels, {min(freq)/1e6:.3f}–{max(freq)/1e6:.3f} MHz, median width={np.median(np.abs(width)):.1f} Hz')
        print('OK', filename)
    except Exception as exc:
        failures.append(filename)
        print('NEEDS INPUT',filename,':',exc)
sys.exit(bool(failures))
