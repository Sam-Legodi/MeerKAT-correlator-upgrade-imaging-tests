#!/usr/bin/env python3
"""Validate explicit local survey inputs in a fresh directory, without imaging.

From the repository/environment: PYTHONPATH=src python scripts/validate_local_survey_matches.py
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import tempfile

from astropy.coordinates import SkyCoord

from meerkat_corr_imaging.audit import run_step_with_audit
from meerkat_corr_imaging.config import load_config
from meerkat_corr_imaging.steps.step6_xmatch import run
from meerkat_corr_imaging.survey_xmatch import CatalogueSpec, resolve_survey_jobs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', default='configs/survey_xmatch/local_mfs.yaml')
    parser.add_argument('--output-root', help='New directory; defaults to a fresh temporary directory')
    parser.add_argument('--racs-wide', help='Also validate this inspected RACS DR1 gausscut export')
    args = parser.parse_args()
    cfg = load_config(args.config)
    if cfg.extra.get('xmatch_pairs'):
        parser.error('Use a survey-only configuration (extra.xmatch_pairs must be empty)')
    if args.output_root:
        root = Path(args.output_root).expanduser().resolve()
        root.mkdir(parents=True, exist_ok=False)
    else:
        root = Path(tempfile.mkdtemp(prefix='mci-local-survey-validation-')).resolve()
    cfg.paths.sky_xmatches_dir = str(root / 'Sky-CrossMatches')
    cfg.paths.reports_dir = str(root / 'reports')
    if args.racs_wide:
        cfg.extra.setdefault('survey_xmatches', []).append({
            'name': 'RACS_low_DR1_wide', 'radius': '10 arcsec',
            'catalogue': {'profile': 'racs_low_dr1_gausscut', 'path': args.racs_wide, 'hdu': 1},
            'association_mode': 'all-candidates', 'coverage': 'unknown'})
    center = SkyCoord('21h47m23.62s', '-81d32m08.6s')
    inputs = {}
    for raw in resolve_survey_jobs(cfg):
        for key in ('input1', 'input2'):
            spec = CatalogueSpec.from_config(raw[key])
            if spec.path in inputs:
                continue
            cat = spec.prepare('T')
            separations = center.separation(cat.coords).deg
            digest = hashlib.sha256()
            with Path(spec.path).open('rb') as handle:
                for chunk in iter(lambda: handle.read(1024 * 1024), b''):
                    digest.update(chunk)
            inputs[spec.path] = {
                'specification': raw[key], 'input_rows': len(cat.table),
                'valid_rows': len(cat.valid_indices), 'sha256': digest.hexdigest(),
                'nearest_to_J2147_deg': float(min(separations)) if len(separations) else None,
                'furthest_from_J2147_deg': float(max(separations)) if len(separations) else None}
    run_step_with_audit('xm_local_validation', cfg.paths.reports_dir, lambda: run(cfg))
    reports = [json.loads(p.read_text()) for p in
               sorted((root / 'Sky-CrossMatches').glob('*.diagnostics.json'))]
    summary = {'config': str(Path(args.config).resolve()), 'inputs': inputs, 'jobs': reports}
    (root / 'validation.json').write_text(json.dumps(summary, indent=2, sort_keys=True) + '\n')
    print(f'Validation: {root / "validation.json"}')
    for report in reports:
        print(Path(report['job']['input1']['path']).name,
              Path(report['job']['input2']['path']).name, report['associations'],
              report['catalogues'][0])


if __name__ == '__main__':
    main()
