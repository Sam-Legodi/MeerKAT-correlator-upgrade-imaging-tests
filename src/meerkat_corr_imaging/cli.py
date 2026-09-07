import argparse
import sys
import os
import re
import shlex
import hashlib
import json
import tempfile
import traceback
from pathlib import Path
from .output_paths import export_pdf
from .audit import run_step_with_audit
from .config import load_config
from .image_pipeline import wire_image_pipeline
from .steps import (
    step2_vis_analysis,
    step3_calibrate_image,
    step4_low_high_slice,
    step5_srcfind,
    step6_xmatch,
    step7_positions,
    step7_flux,
    step8_verification_report,
)


STEP_RUNNERS = {
    "vis": ("step2_visibility_qa", step2_vis_analysis.run),
    "cal": ("step3_calibration_imaging", step3_calibrate_image.run),
    "low_high_slice": ("step4_low_high_slice", step4_low_high_slice.run),
    "src": ("step5_source_finding", step5_srcfind.run),
    "xm": ("step6_cross_matching", step6_xmatch.run),
    "pos": ("step7a_positions", step7_positions.run),
    "flux": ("step7b_flux", step7_flux.run),
    "report": ("step8_verification_report", step8_verification_report.run),
}
ALL_STEPS = tuple(STEP_RUNNERS)
IMAGE_STEPS = ("low_high_slice", "src", "xm", "pos", "flux", "report")
IMAGE_COMMANDS = {"images", "image", "image_all", "all_images"}


def _normalize_cli_args(argv):
    """
    Allow `--config ...` to appear before or after the subcommand by
    rewriting the argv list so argparse always sees it first.
    """
    cfg_tokens = []
    rest = []
    skip_next = False
    for idx, token in enumerate(argv):
        if skip_next:
            skip_next = False
            continue
        if token == "--config":
            if idx + 1 >= len(argv):
                raise SystemExit("error: --config requires a value")
            cfg_tokens = ["--config", argv[idx + 1]]
            skip_next = True
        elif token.startswith("--config="):
            cfg_tokens = [token]
        else:
            rest.append(token)
    return cfg_tokens + rest

def _confirm_session(args, argv):
    if os.environ.get("TMUX") or os.environ.get("STY"):
        return
    label = re.sub(r"[^A-Za-z0-9]", "", Path(args.config).stem)[:5] or "mci"
    digest = hashlib.sha256((str(Path(args.config).resolve()) + args.cmd).encode()).hexdigest()[:4]
    session = label + "-" + digest
    command = shlex.join([sys.executable, "-m", "meerkat_corr_imaging.cli", *argv])
    print("WARNING: This pipeline is running outside tmux or screen.")
    print(f"Start a session: tmux new -s {session}")
    print(f"Or: screen -S {session}")
    print(f"Then rerun: {command}")
    try:
        answer = input("Continue outside tmux/screen? [y/N]: ")
    except (EOFError, OSError):
        answer = ""
    if answer.strip().lower() not in {"y", "yes"}:
        raise SystemExit("Pipeline cancelled; restart inside tmux or screen.")


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    norm_argv = _normalize_cli_args(list(argv))
    p = argparse.ArgumentParser(prog="mci", description="MeerKAT correlator imaging orchestration")
    p.add_argument("--config", required=True, help="Path to master YAML config")
    sub = p.add_subparsers(dest="cmd", required=True)

    sub.add_parser("vis", help="Run visibility QA (step 2)")
    sub.add_parser("cal", help="Run CASA calibration/imaging (step 3)")
    sub.add_parser("low_high_slice", help="Extract low/high planes from MFImage cuboids (step 4)")
    sub.add_parser("src", help="Run PyBDSF source finding (step 5)")
    sub.add_parser("xm",  help="Run cross-matching (step 6)")
    sub.add_parser("pos", help="Run astrometric analysis (step 7a)")
    sub.add_parser("flux",help="Run flux analysis (step 7b)")
    sub.add_parser("report", help="Build consolidated image-domain verification report (step 8)")
    sub.add_parser(
        "all",
        help="Run vis -> cal -> low_high_slice -> src -> xm -> pos -> flux -> report",
    )
    sub.add_parser(
        "images",
        aliases=["image", "image_all", "all_images"],
        help=(
            "Run low_high_slice -> src -> xm -> pos -> flux -> report and pass generated "
            "image products to downstream steps"
        ),
    )

    args = p.parse_args(norm_argv)
    _confirm_session(args, list(argv))
    cfg = load_config(args.config)

    if args.cmd in IMAGE_COMMANDS:
        plan = wire_image_pipeline(cfg)
        print(
            "[IMAGE PIPELINE] Wired "
            f"{len(plan.products)} image product(s), "
            f"{len(plan.xmatch_pairs)} cross-match pair(s), "
            f"{len(plan.positions)} position analysis input(s), and "
            f"{len(plan.flux)} flux analysis input(s)."
        )
        requested_steps = IMAGE_STEPS
    else:
        requested_steps = ALL_STEPS if args.cmd == "all" else (args.cmd,)
    previous = {key: os.environ.get(key) for key in
                ("MCI_REPORT_MANIFEST", "MCI_REPORT_OBSERVATION")}
    failed = False
    with tempfile.TemporaryDirectory(prefix="mci-reports-") as tmp:
        manifest = Path(tmp) / "reports.jsonl"
        os.environ["MCI_REPORT_MANIFEST"] = str(manifest)
        os.environ["MCI_REPORT_OBSERVATION"] = "_vs_".join(t.name for t in cfg.tests) or cfg.reference.name
        try:
            for command_name in requested_steps:
                step_name, runner = STEP_RUNNERS[command_name]
                run_step_with_audit(step_name, cfg.paths.reports_dir,
                                    lambda runner=runner: runner(cfg))
        except Exception:
            traceback.print_exc()
            failed = True
        finally:
            for key, value in previous.items():
                if value is None:
                    os.environ.pop(key, None)
                else:
                    os.environ[key] = value
            reports = sorted(set(json.loads(line) for line in manifest.read_text().splitlines())) if manifest.exists() else []
            produced = []
            for report in reports:
                produced.append(report)
                try:
                    produced.append(str(export_pdf(report)))
                except Exception as exc:
                    print(f"[REPORT] PDF export failed for {report}: {exc}")
                    failed = True
            summary = "[REPORTS] Produced report files:\n" + ("\n".join(produced) or "(none)") + "\n"
            report_dir = Path(cfg.paths.reports_dir).expanduser()
            report_dir.mkdir(parents=True, exist_ok=True)
            (report_dir / "produced_reports.log").write_text(summary, encoding="utf-8")
            print(summary, end="", flush=True)
    if failed:
        raise SystemExit(1)

if __name__ == "__main__":
    main()
