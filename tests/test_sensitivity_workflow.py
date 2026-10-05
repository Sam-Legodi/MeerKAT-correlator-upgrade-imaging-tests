"""Remote workflow contract, exposure scaling and downstream no-MS guarantees."""
import json
from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits
from astropy.table import Table
from docx import Document

from test_thermal_noise import ms  # real casacore fixture; skips only if dependency absent
from meerkat_corr_imaging.config import Config, Target, PathsCfg, load_config
from meerkat_corr_imaging import sensitivity_workflow as sw


def config(tmp_path, ms, products=("mfs", "low", "high"), images=False):
    cfg = Config(project_name="noise-test", reference=Target("CMC1_L"), tests=[Target("GPU_L")],
        paths=PathsCfg(reports_dir=str(tmp_path/"reports"), interim_dir=str(tmp_path/"interim"),
                       sky_xmatches_dir=str(tmp_path/"Sky-CrossMatches")))
    cfg.extra = dict(sensitivity=dict(output_json=str(tmp_path/"sensitivity.json"), products={}),
                     uncertainty=dict(bootstrap_samples=20))
    for product in products:
        sides = {}
        for role in sw.ROLES:
            image = tmp_path/f"{role}_{product}.fits"
            if images:
                header = fits.Header(dict(CTYPE1="RA---TAN", CTYPE2="DEC--TAN", CRVAL1=1., CRVAL2=-30.,
                    CRPIX1=101., CRPIX2=101., CDELT1=-.01, CDELT2=.01,
                    BMAJ=.003, BMIN=.002, BUNIT="Jy/beam", PBCOR=False))
                data = np.random.default_rng(42).normal(0, 1e-5, (201,201)).astype("float32")
                fits.writeto(image, data, header)
            sides[role] = dict(image=str(image), ms=str(ms), field="target", datacolumn="DATA",
                scans="all", channels="all" if product=="mfs" else {"0":[2 if product=="low" else 3]},
                band="L", pb_corrected=False, imaging_context=dict(weighting="briggs", robust=-.5))
        cfg.extra["sensitivity"]["products"][product] = sides
    return cfg


def test_roundtrip_explicit_all_and_no_images(ms,tmp_path):
    cfg = config(tmp_path, ms)
    payload = sw.generate(cfg)
    loaded = sw.load(sw.json_path(cfg), sw.configured_products(cfg))
    assert loaded == payload
    result = loaded["products"]["mfs"]["test"]
    assert result["status"] == "available"
    assert result["scans"] == [1,2]
    assert result["channel_selection"] == {"0":[0,1,2,3]}
    assert result["parallel_hand_exposure_hz_s"] == 240e6
    assert result["image_identity"] is None
    assert not list(tmp_path.glob("sensitivity.json.*.tmp"))
    # Input JSON can be consumed after the MS is gone: no stat/reopen/recalculation.
    import shutil
    shutil.rmtree(ms)
    assert sw.load(sw.json_path(cfg), sw.configured_products(cfg)) == payload


@pytest.mark.parametrize("change", ["version", "units", "missing", "exposure", "digest", "resolved"])
def test_schema_and_stored_selection_rejections(ms,tmp_path,change):
    cfg=config(tmp_path,ms); payload=sw.generate(cfg)
    result=payload["products"]["mfs"]["test"]
    if change=="version": payload["schema_version"]=999
    elif change=="units": payload["units"]="Jy/beam"
    elif change=="missing": result.pop("antenna_ids")
    elif change=="exposure": result["parallel_hand_exposure_hz_s"]*=2
    elif change=="digest": result["association"]["scans"]=[1]
    else: result["selection"]["scans"]=[2]
    sw.atomic_write(sw.json_path(cfg),payload)
    with pytest.raises(ValueError): sw.load(sw.json_path(cfg),sw.configured_products(cfg))


@pytest.mark.parametrize("key,value", [("image","different.fits"), ("scans",[1]), ("channels",{"0":[3]}),
                                       ("pb_corrected",True), ("band","S4")])
def test_requested_association_rejections(ms,tmp_path,key,value):
    cfg=config(tmp_path,ms);sw.generate(cfg)
    cfg.extra["sensitivity"]["products"]["mfs"]["test"][key]=value
    with pytest.raises(ValueError, match="association mismatch"):
        sw.load(sw.json_path(cfg),sw.configured_products(cfg))


def test_missing_selection_unavailable_and_atomic_failure(ms,tmp_path,monkeypatch):
    cfg=config(tmp_path,ms)
    cfg.extra["sensitivity"]["products"]["mfs"]["test"].pop("scans")
    result=sw.generate(cfg)["products"]["mfs"]["test"]
    assert result["status"]=="unavailable" and "scans" in result["reason"]
    old=Path(sw.json_path(cfg)).read_text()
    monkeypatch.setattr(sw.os,"replace",lambda *args: (_ for _ in ()).throw(OSError("write failure")))
    with pytest.raises(OSError): sw.generate(cfg)
    assert Path(sw.json_path(cfg)).read_text()==old
    assert not list(tmp_path.glob("sensitivity.json.*.tmp"))


def test_completed_manifest_identity_and_image_association(ms,tmp_path):
    from meerkat_corr_imaging.paired_astrometry import ms_identity
    cfg=config(tmp_path,ms,products=("mfs",))
    sides=sw.configured_products(cfg)["mfs"]
    spec=sides["test"]
    manifest=tmp_path/"image.manifest.json"
    contract=dict(output_paths=dict(fits=spec["image"]),ms_identity=ms_identity(ms),field_id=0,
        scan_ids=[2],requested_channels=[dict(spw=0,channel=3)],
        tclean_parameters=dict(vis=str(ms),datacolumn="data",weighting="natural"))
    manifest.write_text(json.dumps(dict(processing_status="complete",contract=contract)))
    sides["test"]={k:spec[k] for k in ("image","band","pb_corrected")}
    sides["test"]["manifest"]=str(manifest)
    result=sw.generate(cfg)["products"]["mfs"]["test"]
    assert result["status"]=="available" and result["scans"]==[2]
    assert result["effective_bandwidth_hz"]==1e6
    assert result["imaging_context"]==dict(weighting="natural")
    # A subsequent MS change must invalidate the completed image/MS provenance.
    (ms/"new_provenance_file").write_text("changed")
    result=sw.generate(cfg)["products"]["mfs"]["test"]
    assert result["status"]=="unavailable" and "identity" in result["reason"]
    contract["output_paths"]["fits"]="wrong.fits"
    manifest.write_text(json.dumps(dict(processing_status="complete",contract=contract)))
    assert "image association" in sw.generate(cfg)["products"]["mfs"]["test"]["reason"]


def test_cache_requires_opt_in_and_does_not_mix_exact_exposure(ms,tmp_path):
    import pandas as pd
    cache=tmp_path/"cache";cache.mkdir()
    pd.DataFrame([dict(TIME=t,ANT1=0,ANT2=1,CLASS="cross",SCAN=1,SPW=0,POL=p,FLAG_FRAC=.5)
        for t in (100.,110.) for p in ("XX","YY")]).to_csv(cache/"perrow_amp_stats.csv",index=False)
    pd.DataFrame([dict(CLASS="cross",SPW=0,CHAN=c,FREQ_HZ=1.3e9+c*1e6,POL=p,FLAG_FRAC=.5)
        for c in range(4) for p in ("XX","YY")]).to_csv(cache/"rfi_free_channel_mask.csv",index=False)
    cfg=config(tmp_path,ms,products=("mfs",));sides=sw.configured_products(cfg)["mfs"]
    sides["test"]["visibility_results"]=str(cache)
    result=sw.generate(cfg)["products"]["mfs"]
    assert "allow_approximate" in result["test"]["reason"]
    sides["test"]["allow_approximate"]=True
    result=sw.generate(cfg)["products"]["mfs"]
    assert result["test"]["status"]=="available"
    assert result["test"]["approximation_status"]=="approximate"
    with pytest.raises(ValueError,match="mix exact"):sw.exposure_ratio(result)


def test_product_exposure_scaling(ms,tmp_path):
    cfg=config(tmp_path,ms)
    sw.configured_products(cfg)["low"]["test"]["scans"]=[2]
    result=sw.generate(cfg)["products"]
    assert sw.exposure_ratio(result["mfs"])==1
    assert sw.exposure_ratio(result["low"])==2
    assert sw.exposure_ratio(result["high"])==1
    comparison=sw.noise_comparison(result["low"],10e-6,12e-6)
    assert comparison["expected_gpu_ujy_beam"]==pytest.approx(10*np.sqrt(2))
    assert comparison["measured_gpu_ujy_beam"]==12
    assert comparison["expected_gpu_over_thermal_gpu"]==pytest.approx(
        comparison["expected_gpu_ujy_beam"]/comparison["thermal_gpu_ujy_beam"])


def write_crossmatches(cfg):
    sw.wire_noise_pipeline(cfg)
    table=Table(dict(RA_1=[1.,1.1,1.2],DEC_1=[-30.]*3,RA_2=[1.,1.1,1.2],DEC_2=[-30.]*3,
        Peak_flux_1=[.1,.2,.3],Peak_flux_2=[.11,.21,.31],Total_flux_1=[.1,.2,.3],Total_flux_2=[.11,.21,.31]))
    for p in cfg.extra["positions"]:
        path=Path(p["xmatch_table"]);path.parent.mkdir(parents=True,exist_ok=True);table.write(path)


def forbid_ms(monkeypatch):
    def fail(*args,**kwargs): raise AssertionError("Downstream attempted MS sensitivity")
    monkeypatch.setattr("meerkat_corr_imaging.thermal_noise.ms_sensitivity",fail)
    monkeypatch.setattr("meerkat_corr_imaging.paired_astrometry.ms_identity",fail)
    monkeypatch.setattr(sw,"sensitivity",fail)


def test_flux_consumes_json_without_ms_and_rejects_conflict(ms,tmp_path,monkeypatch):
    from meerkat_corr_imaging.flux_analysis import compare_fluxes_across_band_and_scans
    from meerkat_corr_imaging.uncertainty import UncertaintyConfig
    cfg=config(tmp_path,ms);sw.generate(cfg);write_crossmatches(cfg);forbid_ms(monkeypatch)
    monkeypatch.setenv("MCI_REPORT_MANIFEST",str(tmp_path/"reports.jsonl"))
    kwargs=dict(sensitivity_json=sw.json_path(cfg),sensitivity_products=sw.configured_products(cfg),
                uncertainty_config=UncertaintyConfig(20))
    paths=[cfg.extra["flux"][0][f"ref_{p}_xmatch"] for p in ("low","high","mfs")]
    result=compare_fluxes_across_band_and_scans(*paths,**kwargs)
    assert json.loads(Path(result["metrics_path"]).read_text())["thermal_noise"]["mfs"]["test"]["status"]=="available"
    with pytest.raises(ValueError,match="Conflicting"):
        compare_fluxes_across_band_and_scans(*paths,thermal_noise={},**kwargs)


def test_report_product_ratios_pb_nonpb_and_no_ms(ms,tmp_path,monkeypatch):
    from meerkat_corr_imaging import verification_report as report
    cfg=config(tmp_path,ms,images=True);products=sw.configured_products(cfg)
    products["low"]["test"]["scans"]=[2]
    # MFS PB images; corresponding uncorrected pair with known doubled noise.
    for role in sw.ROLES:
        spec=products["mfs"][role];spec["pb_corrected"]=True
        nonpb=tmp_path/f"{role}_nonpb.fits"
        with fits.open(spec["image"]) as hdus:
            fits.writeto(nonpb,hdus[0].data,hdus[0].header)
        fits.setval(spec["image"],"PBCOR",value=True)
        spec["non_pb_image"]=str(nonpb)
    sw.generate(cfg);write_crossmatches(cfg);forbid_ms(monkeypatch)
    monkeypatch.setenv("MCI_REPORT_MANIFEST",str(tmp_path/"reports.jsonl"))
    # Ensure the old XX proxy cannot be called when precomputed JSON is present.
    monkeypatch.setattr(report,"visibility_exposure",lambda *a:pytest.fail("XX proxy called"))
    path=report.build_verification_report(cfg)
    data=json.loads(path.with_name(path.stem+"_metrics.json").read_text())
    byband={v["band"]:v for v in data["bands"]}
    assert byband["Low band"]["noise_comparison"]["exposure_ratio_reference_over_test"]==2
    assert byband["MFS"]["noise_comparison"]["qualified_pb_comparison"] is True
    assert byband["MFS"]["noise_comparison"]["non_pb_comparison"]["qualified_pb_comparison"] is False
    assert byband["MFS"]["noise_comparison"]["non_pb_comparison"]["status"]=="available"
    assert byband["High band"]["noise_comparison"]["qualified_pb_comparison"] is False
    assert "XX" not in byband["MFS"]["decisions"]["rms"]["assumption"].split("Product-specific")[0]
    doc=Document(path)
    assert any("Natural-weighting thermal-noise" in p.text for p in doc.paragraphs)
    assert any("Measured GPU" in c.text for t in doc.tables for row in t.rows for c in row.cells)


def test_pb_contradiction_and_changed_image(ms,tmp_path):
    cfg=config(tmp_path,ms,images=True);sw.generate(cfg)
    spec=sw.configured_products(cfg)["mfs"]["test"]
    fits.setval(spec["image"],"PBCOR",value=True)
    with pytest.raises(ValueError,match="Image changed"):sw.load(sw.json_path(cfg),sw.configured_products(cfg))
    result=sw.generate(cfg)["products"]["mfs"]["test"]
    assert result["status"]=="unavailable" and "contradicts" in result["reason"]


@pytest.mark.parametrize("reuse",[False,True])
def test_cli_single_launch_runs_only_required_steps(ms,tmp_path,monkeypatch,reuse):
    from meerkat_corr_imaging import cli
    cfg=config(tmp_path,ms)
    if reuse:sw.generate(cfg)
    events=[]
    monkeypatch.setenv("TMUX","test")
    monkeypatch.setattr(cli,"load_config",lambda p:cfg)
    monkeypatch.setattr(cli,"run_step_with_audit",lambda n,d,runner:runner())
    monkeypatch.setattr(cli,"STEP_RUNNERS",{p:(p,lambda c,p=p:events.append(p))
        for p in ("sensitivity","src","xm","pos","flux","report")})
    cli.main(["noise-report","--config","unused.yaml"]+(["--reuse-sensitivity"] if reuse else []))
    assert events==([] if reuse else ["sensitivity"])+["src","xm","pos","flux","report"]
    assert not cfg.low_high_slice.enabled
    assert {p["product"] for p in cfg.extra["positions"]}=={"mfs","low","high"}


def test_json_precedence_and_remote_templates(tmp_path):
    cfg=Config("test",extra=dict(sensitivity=dict(output_json=str(tmp_path/"a.json")),
                                sensitivity_json=str(tmp_path/"b.json")))
    with pytest.raises(ValueError,match="Conflicting"):sw.json_path(cfg)
    root=Path(__file__).parents[1]
    for name in ("l_band","s4"):
        import yaml
        template=yaml.safe_load((root/f"configs/noise_report/{name}.yaml").read_text())
        template["path_vars"]["output_root"]=str(tmp_path/name)
        path=tmp_path/f"{name}.yaml";path.write_text(yaml.safe_dump(template))
        cfg=load_config(path)
        assert set(sw.configured_products(cfg))=={"mfs","low","high"}
        assert sw.json_path(cfg).startswith(str(tmp_path))
        assert sw.configured_products(cfg)["mfs"]["test"]["ms"].startswith("/srv/")
        assert not cfg.casa.imaging_enabled


def test_ms_missing_baselines_and_weight_fallback(ms):
    from casacore.tables import table
    from meerkat_corr_imaging.thermal_noise import ms_sensitivity
    with table(str(ms),readonly=False,ack=False) as tb:
        # Remove 0-2 and 1-2 from the second dump. Union N stays 3, only 4 baseline dumps contribute.
        tb.removerows([4,5])
        tb.removecols("WEIGHT_SPECTRUM")
        weights=tb.getcol("WEIGHT");weights[:,0]=0;tb.putcol("WEIGHT",weights)
    result=ms_sensitivity(dict(ms=str(ms),field="target",channels={"0":[2,3]}))
    assert result["antenna_ids"]==[0,1,2]
    assert result["on_source_integration_s"]==20
    # Four usable baseline dumps * one YY hand * 2 MHz * 10 s.
    assert result["parallel_hand_exposure_hz_s"]==80e6
    assert result["effective_bandwidth_hz"]==pytest.approx(2e6/3)


def test_standalone_cli_does_not_wire_or_find_sources(ms,tmp_path,monkeypatch):
    from meerkat_corr_imaging import cli
    cfg=config(tmp_path,ms)
    monkeypatch.setenv("TMUX","test")
    monkeypatch.setattr(cli,"load_config",lambda p:cfg)
    monkeypatch.setattr(cli,"wire_noise_pipeline",lambda c:pytest.fail("Image wiring called"))
    cli.main(["--config","unused.yaml","sensitivity"])
    assert Path(sw.json_path(cfg)).is_file()
    assert not cfg.reference.images


def test_flux_step_json_precedence_and_no_ms(ms,tmp_path,monkeypatch):
    from meerkat_corr_imaging.steps import step7_flux
    cfg=config(tmp_path,ms);sw.generate(cfg);write_crossmatches(cfg);forbid_ms(monkeypatch)
    commands=[]
    monkeypatch.setattr(step7_flux,"_run",lambda cmd,**kw:commands.append(cmd))
    monkeypatch.setattr(step7_flux,"_thermal_inputs",lambda *a,**kw:pytest.fail("Direct calculation inferred"))
    step7_flux.run(cfg)
    assert "--sensitivity-json" in commands[0] and "--thermal-noise-json" not in commands[0]
    cfg.extra["flux"][0]["thermal_noise"]={}
    with pytest.raises(Exception,match="flux analysis"):step7_flux.run(cfg)


def test_old_cached_inference_is_opt_in(tmp_path):
    from meerkat_corr_imaging.steps.step7_flux import _thermal_inputs
    cfg=Config("cache",reference=Target("CMC1_S4"),tests=[Target("GPU_S4")],
               extra=dict(verification_report=dict(reference_visibility_results="refcache",test_visibility_results="gpucache")))
    assert not _thermal_inputs(cfg,{},infer=True)["mfs"]
    result=_thermal_inputs(cfg,dict(allow_approximate=True),infer=True)
    assert result["mfs"]["test"]["allow_approximate"] is True


def test_sensitivity_audit_counts_products_not_ms_access(ms, tmp_path, monkeypatch):
    from meerkat_corr_imaging.audit import run_step_with_audit
    cfg = config(tmp_path, ms)
    run_step_with_audit('sensitivity_metadata', tmp_path/'audit', lambda: sw.generate(cfg))
    log = next((tmp_path/'audit'/'pipeline_audits').glob('*.log')).read_text()
    assert 'Inputs: 6 succeeded, 0 failed, 0 not run' in log
    assert 'mfs/test:' in log
    cfg.extra['sensitivity']['products']['mfs']['test'].pop('scans')
    run_step_with_audit('sensitivity_metadata', tmp_path/'unavailable', lambda: sw.generate(cfg))
    log = next((tmp_path/'unavailable'/'pipeline_audits').glob('*.log')).read_text()
    assert 'Inputs: 5 succeeded, 1 failed, 0 not run' in log
