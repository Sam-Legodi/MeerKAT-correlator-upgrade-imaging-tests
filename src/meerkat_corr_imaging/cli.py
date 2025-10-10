import argparse
import sys
from .config import load_config
from .steps import (
    step2_vis_analysis,
    step3_calibrate_image,
    step5_srcfind,
    step6_xmatch,
    step7_positions,
    step7_flux,
)

def main(argv=None):
    argv = argv or sys.argv[1:]
    p = argparse.ArgumentParser(prog="mci", description="MeerKAT correlator imaging orchestration")
    p.add_argument("--config", required=True, help="Path to master YAML config")
    sub = p.add_subparsers(dest="cmd", required=True)

    sub.add_parser("vis", help="Run visibility QA (step 2 / step 4)")
    sub.add_parser("cal", help="Run CASA calibration/imaging (step 3)")
    sub.add_parser("src", help="Run PyBDSF source finding (step 5)")
    sub.add_parser("xm",  help="Run cross-matching (step 6)")
    sub.add_parser("pos", help="Run astrometric analysis (step 7a)")
    sub.add_parser("flux",help="Run flux analysis (step 7b)")
    sub.add_parser("all", help="Run vis -> cal -> src -> xm -> pos -> flux")

    args = p.parse_args(argv)
    cfg = load_config(args.config)

    if args.cmd == "vis":
        step2_vis_analysis.run(cfg)
    elif args.cmd == "cal":
        step3_calibrate_image.run(cfg)
    elif args.cmd == "src":
        step5_srcfind.run(cfg)
    elif args.cmd == "xm":
        step6_xmatch.run(cfg)
    elif args.cmd == "pos":
        step7_positions.run(cfg)
    elif args.cmd == "flux":
        step7_flux.run(cfg)
    elif args.cmd == "all":
        step2_vis_analysis.run(cfg)
        step3_calibrate_image.run(cfg)
        step5_srcfind.run(cfg)
        step6_xmatch.run(cfg)
        step7_positions.run(cfg)
        step7_flux.run(cfg)

if __name__ == "__main__":
    main()
