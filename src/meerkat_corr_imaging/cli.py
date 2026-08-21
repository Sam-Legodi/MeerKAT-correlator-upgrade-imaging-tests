import argparse
import sys
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
)


STEP_RUNNERS = {
    "vis": ("step2_visibility_qa", step2_vis_analysis.run),
    "cal": ("step3_calibration_imaging", step3_calibrate_image.run),
    "low_high_slice": ("step4_low_high_slice", step4_low_high_slice.run),
    "src": ("step5_source_finding", step5_srcfind.run),
    "xm": ("step6_cross_matching", step6_xmatch.run),
    "pos": ("step7a_positions", step7_positions.run),
    "flux": ("step7b_flux", step7_flux.run),
}
ALL_STEPS = tuple(STEP_RUNNERS)
IMAGE_STEPS = ("low_high_slice", "src", "xm", "pos", "flux")
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

def main(argv=None):
    argv = argv or sys.argv[1:]
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
    sub.add_parser(
        "all",
        help="Run vis -> cal -> low_high_slice -> src -> xm -> pos -> flux",
    )
    sub.add_parser(
        "images",
        aliases=["image", "image_all", "all_images"],
        help=(
            "Run low_high_slice -> src -> xm -> pos -> flux and pass generated "
            "image products to downstream steps"
        ),
    )

    args = p.parse_args(norm_argv)
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
    for command_name in requested_steps:
        step_name, runner = STEP_RUNNERS[command_name]
        run_step_with_audit(
            step_name,
            cfg.paths.reports_dir,
            lambda runner=runner: runner(cfg),
        )

if __name__ == "__main__":
    main()
