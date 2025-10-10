# standalone_xxyy_solve.py
# Standalone CASA calibration script: delays -> bandpass -> gains -> flux bootstrapping -> leakage (D-terms) -> applycal
# Runs in CASA 5.x (Py2.7) and CASA 6.x (Py3)
# Usage:
#   casa --log2term --nologger --nogui -c standalone_xxyy_solve.py
# or
#   casa -c standalone_xxyy_solve.py

from __future__ import print_function
import os, sys, shutil, logging
import numpy as np
from datetime import datetime
from time import gmtime

# ----------------------------- CASA COMPAT ---------------------------------
# Do NOT import casatasks in CASA 5; tasks are globals. We only resolve tools.
# msmetadata tool
msmd = None
try:
    # CASA 6-style tools
    from casatools import msmetadata as _msmetadata
    msmd = _msmetadata()
except Exception:
    try:
        # CASA 5-style tools
        from taskinit import msmdtool
        msmd = msmdtool()
    except Exception:
        msmd = None  # will error when used

# ----------------------------- LOGGING -------------------------------------
logging.Formatter.converter = gmtime
logger = logging.getLogger(__name__)
logging.basicConfig(format="%(asctime)-15s %(levelname)s: %(message)s",
                    level=logging.INFO)

def _get_casalog():
    """Retrieve CASA logger without shadowing."""
    # CASA 6
    try:
        from casatools import casalog as tool
        return tool
    except Exception:
        pass
    # CASA 5
    try:
        from taskinit import casalog as tool
        return tool
    except Exception:
        pass
    # Fallback to injected global
    try:
        import __main__
        return getattr(__main__, 'casalog', None)
    except Exception:
        return None

def _log(msg):
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    tool = _get_casalog()
    if tool is not None:
        try:
            tool.post(msg)
            return
        except Exception:
            pass
    print("[{}] {}".format(ts, msg))

def _set_logfile(path):
    tool = _get_casalog()
    if tool is not None:
        try:
            tool.setlogfile(path)
            return
        except Exception:
            pass
    _log("Logging to {} (print fallback)".format(path))

def out(path, outdir, make_outdir):
    return os.path.join(outdir, path) if make_outdir else path

# ----------------------------- INLINE setjy.py -----------------------------
def linfit(xInput, xDataList, yDataList):
    """
    Linear fit helper (exact from setjy.py).
    """
    y_predict = np.poly1d(np.polyfit(xDataList, yDataList, 1))
    yPredict = y_predict(xInput)
    return yPredict

def do_setjy(visname, spw, fields, standard, dopol=False, createmms=True):
    """
    Exact implementation from your setjy.py.
    """
    # Use global CASA task 'delmod'
    delmod(vis=visname)  # clear existing model (prevents exit code 1)

    fluxlist = ["J0408-6545", "0408-6545", ""]
    ismms = createmms

    if msmd is None:
        raise RuntimeError("msmetadata tool unavailable in this CASA build.")
    msmd.open(visname)
    fnames = fields.fluxfield.split(",")
    for fname in fnames:
        if fname.isdigit():
            fname = msmd.namesforfields(int(fname))

    do_manual = False
    for ff in fluxlist:
        if ff in fnames:
            setjyname = ff
            do_manual = True
            break
        else:
            setjyname = fields.fluxfield.split(",")[0]

    if do_manual:
        smodel = [17.066, 0.0, 0.0, 0.0]
        spix = [-1.179]
        reffreq = "1284MHz"

        logger.info("Using manual flux density scale - ")
        logger.info("Flux model: %s ", smodel)
        logger.info("Spix: %s", spix)
        logger.info("Ref freq %s", reffreq)

        setjy(vis=visname, field=setjyname, scalebychan=True, standard="manual",
              fluxdensity=smodel, spix=spix, reffreq=reffreq, ismms=ismms)
    else:
        setjy(vis=visname, field=setjyname, spw=spw, scalebychan=True,
              standard=standard, ismms=ismms)

    fieldnames = msmd.fieldnames()

    if dopol:
        # Check if 3C286 exists in the data
        is3C286 = False
        try:
            calibrator_3C286 = list(set(["3C286", "1328+307", "1331+305", "J1331+3030"]).intersection(set(fieldnames)))[0]
        except IndexError:
            calibrator_3C286 = []

        if len(calibrator_3C286):
            is3C286 = True
            id3C286 = str(msmd.fieldsforname(calibrator_3C286)[0])

        if is3C286:
            logger.info("Detected calibrator name(s):  %s" % calibrator_3C286)
            logger.info("Flux and spectral index taken/calculated from:  https://science.nrao.edu/facilities/vla/docs/manuals/oss/performance/fdscale")
            logger.info("Estimating polarization index and position angle of polarized emission from linear fit based on: Perley & Butler 2013 (https://ui.adsabs.harvard.edu/abs/2013ApJS..204...19P/abstract)")
            # central freq of spw
            spwMeanFreq = msmd.meanfreq(0, unit='GHz')
            freqList = np.array([1.05, 1.45, 1.64, 1.95])
            # fractional linear polarisation
            fracPolList = [0.086, 0.095, 0.099, 0.101]
            polindex = linfit(spwMeanFreq, freqList, fracPolList)
            logger.info("Predicted polindex at frequency %s: %s", spwMeanFreq, polindex)
            # position angle of polarized intensity
            polPositionAngleList = [33, 33, 33, 33]
            polangle = linfit(spwMeanFreq, freqList, polPositionAngleList)
            logger.info("Predicted pol angle at frequency %s: %s", spwMeanFreq, polangle)

            reffreq = "1.45GHz"
            logger.info("Ref freq %s", reffreq)
            setjy(vis=visname,
                field=id3C286,
                scalebychan=True,
                standard="manual",
                fluxdensity=[-14.6, 0.0, 0.0, 0.0],
                #spix=-0.52, # between 1465MHz and 1565MHz
                reffreq=reffreq,
                polindex=[polindex],
                polangle=[polangle],
                rotmeas=0,ismms=ismms)

        # Check if 3C138 exists in the data
        is3C138 = False
        try:
            calibrator_3C138 = list(set(["3C138", "0518+165", "0521+166", "J0521+1638"]).intersection(set(fieldnames)))[0]
        except IndexError:
            calibrator_3C138 = []

        if len(calibrator_3C138):
            is3C138 = True
            id3C138 = str(msmd.fieldsforname(calibrator_3C138)[0])

        if is3C138:
            logger.info("Detected calibrator name(s):  %s" % calibrator_3C138)
            logger.info("Flux and spectral index taken/calculated from:  https://science.nrao.edu/facilities/vla/docs/manuals/oss/performance/fdscale")
            logger.info("Estimating polarization index and position angle of polarized emission from linear fit based on: Perley & Butler 2013 (https://ui.adsabs.harvard.edu/abs/2013ApJS..204...19P/abstract)")
            # central freq of spw
            spwMeanFreq = msmd.meanfreq(0, unit='GHz')
            freqList = np.array([1.05, 1.45, 1.64, 1.95])
            # fractional linear polarisation
            fracPolList = [0.056, 0.075, 0.084, 0.09]
            polindex = linfit(spwMeanFreq, freqList, fracPolList)
            logger.info("Predicted polindex at frequency %s: %s", spwMeanFreq, polindex)
            # position angle of polarized intensity
            polPositionAngleList = [-14, -11, -10, -10]
            polangle = linfit(spwMeanFreq, freqList, polPositionAngleList)
            logger.info("Predicted pol angle at frequency %s: %s", spwMeanFreq, polangle)

            reffreq = "1.45GHz"
            logger.info("Ref freq %s", reffreq)
            setjy(vis=visname,
                field=id3C138,
                scalebychan=True,
                standard="manual",
                fluxdensity=[-8.26, 0.0, 0.0, 0.0],
                #spix=-0.57,  # between 1465MHz and 1565MHz
                reffreq=reffreq,
                polindex=[polindex],
                polangle=[polangle],
                rotmeas=0,ismms=ismms)

    msmd.done()

# ----------------------------- USER INPUTS ---------------------------------

# Minimal shim to satisfy do_setjy(fields=...) without external config modules
class Fields(object):
    def __init__(self, fluxfield):
        self.fluxfield = fluxfield  # comma-separated CASA field selector

# MeasurementSet to calibrate
visname = "data/reference_obs/1692135074_sdp_l0.full.ms"
antennas = ""
refant   = "m060"

# Fields
flux_field   = "J0408-6545"
bp_field     = "J0408-6545"
delay_field  = "J0408-6545"
gain_fields  = "J1619-8418"
target_fields= "J2147-8132"
apply_to_all = False  # applycal to all fields (True) or only listed fields (False)

# Solution intervals
delay_solint  = "inf"
bp_solint     = "inf"
phase_solint  = "int"
amp_solint    = "inf"

# Combine
bp_combine    = "scan"
phase_combine = ""
amp_combine   = "scan"

# Leakage calibration
do_leakage      = True
leakage_field   = flux_field
leakage_poltype = "D"
leakage_solint  = "inf"

# Output
# Derive clean output paths so filenames dont duplicate directories
parentdir = os.path.dirname(visname)
basename  = os.path.basename(visname)

# Prefix used for all caltable filenames (no directories here)
out_root  = basename.replace(".ms", "cal")
# Output directory lives next to the MS
outdir    = os.path.join(parentdir, out_root + "_" + datetime.now().strftime("%Y%m%d_%H%M%S"))
make_outdir = True  # ensure we write into outdir via out()

listobs_txt = False
do_split    = False
split_ms    = "target_split.ms"
split_datacolumn = "corrected"

# Applycal
interp_all  = ["linear","linear","linear","linear"]
calwt       = False
apply_mode  = "calonly"

# Flagging
do_initial_flagging = False
initial_flagging_args = dict(
    vis=visname,
    mode="manual",
    autocorr=True
)

# ----------------------------- PREP & LOGGING ------------------------------

if make_outdir and not os.path.isdir(outdir):
    os.makedirs(outdir)
logfile = os.path.join(outdir if make_outdir else ".", out_root + ".casa.log")
_set_logfile(logfile)

_log("==============================================")
_log(" CASA Calibration Script: delays -> bandpass -> gains -> flux bootstrapping -> leakage (D-terms) -> applycal")
_log(" Start: {}".format(datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
_log(" MS: {}".format(visname))
_log(" Outdir: {}".format(outdir if make_outdir else "(cwd)"))
_log("==============================================")

# Output caltables
ct_delay    = out(out_root + ".K.delays.cal",    outdir, make_outdir)
ct_gp_prebp = out(out_root + ".Gpre.phase.cal",  outdir, make_outdir)
ct_bandpass = out(out_root + ".B.bandpass.cal",  outdir, make_outdir)
ct_gphase   = out(out_root + ".G.phase.cal",     outdir, make_outdir)
ct_gamp     = out(out_root + ".G.amp.cal",       outdir, make_outdir)
ct_flux     = out(out_root + ".fluxscale.cal",   outdir, make_outdir)
ct_leakage  = out(out_root + ".D.leakage.cal",   outdir, make_outdir)

# ----------------------------- PIPELINE ------------------------------------

# Listobs
if listobs_txt:
    listfile = out(out_root + ".listobs.txt", outdir, make_outdir)
    _log("listobs -> {}".format(listfile))
    listobs(vis=visname, listfile=listfile, overwrite=True, verbose=False)

# Flagging
if do_initial_flagging:
    _log("Initial flagging...")
    flagdata(**initial_flagging_args)

# setjy (exact behavior from your setjy.py)
_log("Running do_setjy...")
do_setjy(visname=visname,
         spw="",
         fields=Fields(flux_field),
         standard="Stevens-Reynolds 2016",
         dopol=False,
         createmms=True)

# Delay solve
_log("gaincal: delays (K)...")
gaincal(vis=visname, caltable=ct_delay, field=delay_field,
        solint=delay_solint, refant=refant, gaintype="K", calmode="p",
        minblperant=4, minsnr=3.0)

# Pre-bandpass phase
_log("gaincal: pre-bandpass phase...")
gaincal(vis=visname, caltable=ct_gp_prebp, field=bp_field,
        solint=phase_solint, combine=phase_combine, refant=refant,
        gaintype="G", calmode="p", minblperant=4, minsnr=3.0,
        gaintable=[ct_delay], parang=True)

# Bandpass
_log("bandpass...")
bandpass(vis=visname, caltable=ct_bandpass, field=bp_field,
         solint=bp_solint, combine=bp_combine, refant=refant,
         minblperant=4, minsnr=3.0, gaintable=[ct_delay, ct_gp_prebp],
         parang=True)

# Phase-only gains
_log("gaincal: phase-only...")
gaincal(vis=visname, caltable=ct_gphase, field=gain_fields,
        solint=phase_solint, combine=phase_combine, refant=refant,
        gaintype="G", calmode="p", minblperant=4, minsnr=3.0,
        gaintable=[ct_delay, ct_bandpass], parang=True)

# Amp+phase gains
_log("gaincal: flux_field and gain_field amp+phase...")
# Include BOTH flux_field and gain_fields so fluxscale can find the reference in ct_gamp
gaincal(vis=visname, caltable=ct_gamp, field=",".join([flux_field, gain_fields]),
        solint=amp_solint, combine=amp_combine, refant=refant,
        gaintype="G", calmode="ap", minblperant=4, minsnr=3.0,
        gaintable=[ct_delay, ct_bandpass, ct_gphase], parang=True)

# Fluxscale
_log("fluxscale...")
fluxscale(vis=visname, caltable=ct_gamp, fluxtable=ct_flux,
          reference=flux_field, transfer=gain_fields,
          listfile=out(out_root + ".fluxscale.txt", outdir, make_outdir))

# Leakage (D-terms only)
if do_leakage:
    _log("polcal: leakage-only...")
    polcal(vis=visname, caltable=ct_leakage, field=leakage_field,
           solint=leakage_solint, poltype=leakage_poltype,
           gaintable=[ct_delay, ct_bandpass, ct_gphase, ct_flux])

# Applycal
_log("applycal...")
gaintables = [ct_delay, ct_bandpass, ct_gphase, ct_flux]
if do_leakage:
    gaintables.append(ct_leakage)
interp = (interp_all + ["linear"] * 10)[:len(gaintables)]

if apply_to_all:
    _log("  applying to all fields")
    applycal(vis=visname,
             field=",".join([flux_field, gain_fields, target_fields]),
             gaintable=gaintables,
             interp=interp,
             calwt=calwt,
             parang=True,
             applymode=apply_mode)
else:
    _log("  applying to select fields: {}".format(",".join([flux_field, gain_fields])))
    applycal(vis=visname,
             field=",".join([flux_field, gain_fields]),
             gaintable=gaintables,
             interp=interp,
             calwt=calwt,
             parang=True,
             applymode=apply_mode)

# Optional split of targets
if do_split:
    _log("split -> {}".format(split_ms))
    split(vis=visname,
          outputvis=out(split_ms, outdir, make_outdir),
          datacolumn=split_datacolumn,
          field=target_fields,
          keepflags=False)

_log("==============================================")
_log(" Calibration complete.")
_log(" End: {}".format(datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
_log(" Log: {}".format(logfile))
_log("==============================================")
