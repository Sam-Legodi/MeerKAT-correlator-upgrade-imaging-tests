from meerkat_corr_imaging.visibility_cache import record_reference_completion, reference_is_complete


def test_reference_cache_requires_complete_unchanged_outputs_and_source(tmp_path):
    ms = tmp_path / 'ref.ms'
    ms.mkdir()
    source = ms / 'table.dat'
    source.write_bytes(b'visibility data')
    output = tmp_path / 'result'
    output.mkdir()
    assert not reference_is_complete(ms, output)
    names = ['scan_averaged_amp_stats.csv', 'acceptance_summary.csv',
             'amplitude_by_baseline.csv', 'flagging_by_channel.csv',
             'draft_vis_amp_summary_epoch.docx']
    for name in names:
        (output / name).write_bytes(b'completed output')
    record_reference_completion(ms, output)
    assert reference_is_complete(ms, output)
    (ms / 'table.lock').write_bytes(b'lock')
    assert reference_is_complete(ms, output)
    (output / names[0]).unlink()
    assert not reference_is_complete(ms, output)
    (output / names[0]).write_bytes(b'completed output')
    source.write_bytes(b'changed visibilities')
    assert not reference_is_complete(ms, output)


def test_step_reuses_reference_but_runs_test(tmp_path, monkeypatch):
    from meerkat_corr_imaging.steps import step2_vis_analysis as step
    from meerkat_corr_imaging.config import Config, PathsCfg, Target
    cfg = Config(project_name='cache', paths=PathsCfg(interim_dir=str(tmp_path)),
                 reference=Target(name='ref', ms_paths=['ref.ms']),
                 tests=[Target(name='test', ms_paths=['test.ms'])])
    monkeypatch.setattr(step, 'reference_is_complete', lambda *a: True)
    calls = []
    monkeypatch.setattr(step, '_run_cmd', lambda cmd, **kw: calls.append(cmd))
    step.run(cfg)
    assert len(calls) == 1
    assert calls[0][calls[0].index('--ms') + 1] == 'test.ms'
    assert calls[0][calls[0].index('--ms-ref-results') + 1] == str(tmp_path / 'ref')
