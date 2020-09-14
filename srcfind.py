import sys, time, os, string, shutil
from scipy import *
from matplotlib.pyplot import *
from astropy.io import fits
from astropy.wcs import WCS as w
from astropy.io.votable import parse_single_table
from matplotlib import pyplot as plt
from astropy.table import Table,Column

try:
    import bdsf
except ImportError as ioe:
    print(' >>> ', ioe)
    pass

import katdal, csv
import traceback, atpy
from astropy.wcs import WCS
import numpy as np
from scipy import stats
from scipy import optimize
import matplotlib.cm as cm
import astropy.units as u
import scipy.interpolate as si
import matplotlib.colors as plc
import matplotlib.colorbar as plo
import matplotlib.colors as mpc
import matplotlib.pyplot as mpl
import matplotlib.colorbar as mpo
import mpl_toolkits.axisartist as AA
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1 import host_subplot
from astropy.io.votable import parse_single_table
from astropy.cosmology import FlatLambdaCDM
import matplotlib.gridspec as gridspec
from scipy import interpolate
from scipy.optimize import curve_fit
from astropy.stats import poisson_conf_interval
from astropy.coordinates import SkyCoord as sk
from PyAstronomy import pyasl
#from recipes.almapolhelpers import *

#################################################################################

class source:
    def __init__(self,ra,raerr,dec,decerr,flux,dflux,peak,dpeak,Smaj,dmaj,Smin,dmin,Spa,dpa,rms,x,y):
        self.ra = ra
        self.dra = raerr
        self.dec = dec
        self.ddec = decerr
        self.flux = flux
        self.df = dflux   
        self.peak = peak
        self.dp = dpeak 
        self.Smaj = Smaj
        self.dmaj = dmaj
        self.Smin = Smin
        self.dmin = dmin
        self.Spa = Spa
        self.dpa = dpa
        self.rms = rms
        self.x = x
        self.y = y
        
class beams:
    def __init__(self,bmaj,bmin,bpa):
        self.maj = bmaj
        self.min = bmin
        self.pa  = bpa
    def show(self):
        print(("bmaj = %1.3e, bmin = %1.3e, bpa = %1.3e" % (self.maj, self.min, self.pa)))
        
def source_errors(Smaj,Smin,Spa,flux,peak,rms,beam): 
    # calculate source Gaussian fit error using the Condon et al prescription
    rho = sqrt( pi/(8.0*log(2.0)) * Smaj*Smin/(bmin*bmin/4.0) \
                * (peak*peak)/(rms*rms) )
    dA = sqrt(2.0)*(peak/rho)
    dI = sqrt(2.0)*(flux/rho)
    dx = sqrt(1.0/(4.0*log(2.0)) * (Smaj/rho))
    dy = sqrt(1.0/(4.0*log(2.0)) * (Smin/rho))
    dalpha = sqrt(dx*dx*sin(Spa*pi/180.0)*sin(Spa*pi/180.0) + \
                  dy*dy*cos(Spa*pi/180.0)*cos(Spa*pi/180.0))
    ddelta = sqrt(dx*dx*cos(Spa*pi/180.0)*cos(Spa*pi/180.0) +\
                  dy*dy*sin(Spa*pi/180.0)*sin(Spa*pi/180.0)) 
    return(dI,dA,dalpha/3600.0,ddelta/3600.0)


def write_myipac(ipacfilename, ipactable):
    try:
        ipactable.write(ipacfilename,type='ipac')
        print(('file: %12s written out ' %ipacfilename))
    except Exception:
        print(('   File (%s) exists! ' % ipacfilename))
        overwrite = 'y' #raw_input('   overwrite the existing file?   ')
        if (overwrite.lower() == 'y'):
            os.system('rm %s' % ipacfilename)
            print(('   [!] file: %s  removed .. will write a new version now.' % ipacfilename))
            ipactable.write(ipacfilename,type='ipac')
        else:
            print(('   [!] file: %s left as is.' % ipacfilename))
    return


def deconvolve(Smaj,Smin,Spa,dmaj,dmin,dpa,beam):
    # deconvolve beam from source dimensions. Input angles in degrees. Returns source dimensions.
    bmaj = beam.maj; bmin = beam[1]; bpa = beam[2]
    diffpa = (Spa - beam.bpa) * (pi/90.0)  
    rho_c = (Smaj*Smaj - Smin*Smin)*cos(diffpa) - (beam.maj*beam.maj - beam.min*beam.min)
    if(rho_c == 0.0):
        sigic2 = 0.0
        rho_A = 0.0
    else:
        sigic2 = arctan( (Smaj*Smaj - Smin*Smin)*sin(diffpa)/rho_c)
        rho_A = ((beam.maj*beam.maj - beam.min*beam.min) - (Smaj*Smaj - Smin*Smin)*cos(diffpa))/(2.0*cos(sigic2))
    Rpa = sigic2 * (90.0/pi) + bpa
    det = ((Smaj*Smaj + Smin*Smin) - (beam.maj*beam.maj + beam.min*beam.min))/2.0
    Imaj = det - rho_A
    Imin = det + rho_A

    ierr = 0
    if (Imin < 0.0):
        ierr = ierr + 1
        Imin = 0.0
    if (Imaj < 0.0):
        ierr = ierr + 1
        Imaj = 0.0
    Imaj = sqrt(fabs(Imaj))
    Imin = sqrt(fabs(Imin))
    if(Imaj < Imin):           # swap if Imaj < Imin
        temp = Imaj
        Imaj = Imin
        Imin = temp
        Rpa = Rpa + 90.0
    while(Rpa > 180.0):
        Rpa = Rpa - 180.0
    while(Rpa < -180.0):
        Rpa = Rpa + 180.0
    if (Imaj == 0.0):
        Rpa = 0.0
    return(Imaj,Imin,Rpa,ierr)   

def quick(items, int, count):
    qs(items,0,count);

def pbcorr(x,y,cell,FWHM):
    r = fabs(cell)*sqrt(x*x + y*y)
    r2 = r*r; FWHM2 = FWHM*FWHM
    return( exp(-2.7726*r2/FWHM2) )

def evenpoly(X1, X2, C):
    Y = []
    polyorder = 2
    if X2 == None:
        R = X1
    else:
        R = sqrt(X1*X1 + X2*X2)
    
    for i0 in range(len(R)):
        yi = 0.
        for i1 in range(polyorder):
            yi = yi + C*R[i0]**(2*i1)
        Y.append(yi)
    return(array(Y)/max(Y))

def safe_div(x,y):
    try:
        if y == 0.:
            return 0
    
    except ValueError:
        for yi in y:
            if yi == 0.:
                return 0
    return (x/y)

def qs(items, left, right):
    i = left; j=right;
    x = items[(left+right)/2].ra;

    while True:
        while ((items[i].ra < x) and (i < right)):
            i= i+1
        while ((items[j].ra > x) and (j > left) ): 
            j= j-1
        if(i <= j):
            temp = items[i]
            items[i] = items[j]
            items[j] = temp
            i = i+1
            j = j-1
        if not(i <= j):
            break
    if(left < j):
        qs(items,left,j);
    if(i < right):
        qs(items,i,right);

def complex_amp_phase(A):
    try:
        amp  = sqrt(A.real**2 + A.imag**2)
        phse = phase(A)
    except TypeError as te0:
        print(te0)
        amp, phse = [], []
        for i in range(len(A)):
            amp.append( sqrt(A.real**2 + A.imag**2) )
            phse.append( phase(A) )
    amp, phse = np.array(amp), np.array(phse)
    return (amp, phse)

def get_refant(MS_in,fluxfield):
    # must be run in CASA ...
    # returns best refant based on the deviation of it's amplitude from the median of all antenna amplitiudes.
    # plots antenna amplitudes vs. antenna ID and saves plot in same directory as MS file.
    import numpy as np
    import matplotlib.pyplot as plt
    msmd.open(MS_in)
    fluxscans = msmd.scansforfield(int(fluxfield))
    antennas  = msmd.antennasforscan(fluxscans[0])
    numchans  = msmd.nchan(0)
    gainchanstart = int(0.2*numchans)
    gainchanstop  = int(0.8*numchans)
    gainspw       = '0:'+str(gainchanstart)+'~'+str(gainchanstop) # gainspw for gaincal in selfcal
    msmd.done()
    antamp =[]; antrms = []
    for ant in antennas:
        ant = str(ant)
        t = visstat(vis=MS_in,field = fluxfield,antenna = ant, timeaverage=True, timebin='500min',
                    timespan='state,scan',reportingaxes='field')
        item = str(list(t.keys())[0])
        amp = float(t[item]['median'])
        rms = float(t[item]['rms'])
        antamp.append(amp)
        antrms.append(rms)
    antamp = np.array(antamp)
    antrms = np.array(antrms)
    medamp = np.median(antamp)
    medrms = np.median(antrms)
    goodrms= []; goodamp=[]; goodant=[]
    for i in range(len(antamp)):
        if (antamp[i] > medamp):
            goodant.append(int(antennas[i]))
            goodamp.append(antamp[i])
            goodrms.append(antrms[i])
    goodrms = np.array(goodrms)
    
    plt.plot(goodant, goodamp, '.')
    plt.axhline(y = medamp, ls='-.', c='g', label = 'median')
    #plt.axhline(y = medamp+medrms, ls = ':', c='r', alpha=.5, label = 'median$\pm\sigma_{rms}$')
    #plt.axhline(y = medamp-medrms, ls = ':', c='r', alpha=.5)
    plt.grid(True,ls="-.", alpha=0.35)
    plt.xlabel('Antenna ID'); plt.ylabel('Amp')
    plt.legend(loc=0)#, fancybox=True, framealpha=0.35)
    plt.savefig('%s-antennas.png'%MS_in)
    plt.show()
    
    j = np.argmin(goodrms)
    referenceant = "%03d"%int(goodant[j])
    return (referenceant, gainspw)

def check_directory(dir_in):
    # run this to check if output directory exists to avoid file conflicts
    if os.path.isdir(dir_in+"/"):
        dir_new = '%s_new' % dir_in
    elif os.path.isdir(dir_in+"_new/"):
        dir_new = '%s_new2' % dir_in
    else:
        dir_new = dir_in
    return dir_new

def separation(ra1,dec1,ra2,dec2):
    """This function returns the angular separation between two points defined by
    ra1,dec1, and ra2,dec2. 
    Input and output is in degrees """
    from PyAstronomy import pyasl
    pi = 3.1415928
    deg2rad = pi/180.0
    sep = pyasl.getAngDist(float(ra1), float(dec1), float(ra2), float(dec2))
    return(sep)

def get_trgt_gcal_sep(TARGETS_in, GCALS_in, output_dir):
    # chooses TARGET list based on angular distance between targets and gain cal.
    from astropy.coordinates import SkyCoord as sk
    import numpy as np
    import matplotlib.pyplot as plt
    
    if output_dir == None:
        output_dir = 'target_gcal_seps'
        if not (os.path.isdir(output_dir)):
            os.mkdir(output_dir)
    print(('\nTarget - gaincal separation results are in: {}\n'.format(output_dir)))
    outfile = open(output_dir+'/target_gcal_dist.txt','w',0)
    RAs, DECs, Dr, out_TARGET = [],[],[],[]
    GCALS_in   = GCALS_in[0]
    cal_radius = 16. 
    for src in TARGETS_in:
        for gcal in GCALS_in:
            try:
                ra_src0  = src[1:7]; dec_src0 = src[7:]
                ra_gcal0 = gcal[1:7]; dec_gcal0 = gcal[7:]
                ra_src1  = str( ra_src0[0:2]+'h'+ra_src0[2:4]+'m'+ra_src0[4:]+'s' )
                dec_src1 = str( dec_src0[0:3]+'d'+dec_src0[3:5]+'m'+dec_src0[5:]+'s' )
                ra_gcal1 = str( ra_gcal0[0:2]+'h'+ra_gcal0[2:4]+'m'+ra_gcal0[4:]+'s' )
                dec_gcal1= str( dec_gcal0[0:3]+'d'+dec_gcal0[3:5]+'m'+dec_gcal0[5:]+'s' )
                
                D_src   = sk(ra=ra_src1,dec=dec_src1,frame='icrs')
                D_gcal  = sk(ra=ra_gcal1,dec=dec_gcal1,frame='icrs')
                
                ra_src  = D_src.ra.deg; dec_src = D_src.dec.deg
                ra_gcal = D_gcal.ra.deg; dec_gcal = D_gcal.dec.deg
                seprtn  = separation(ra_src, dec_src, ra_gcal, dec_gcal)
                RAs.append(ra_src), DECs.append(dec_src), Dr.append(seprtn)
                if abs(seprtn) <= cal_radius and src != gcal:
                    outfile.write( '*** TARGET-GCAL separation: {} deg (TARGET: {}, GCAL: {})\n'.format(seprtn,src,gcal) )
                    print(('TARGET-GCAL separation: {} deg (TARGET: {}, GCAL: {})'.format(seprtn,src,gcal)))
                    #print 'TARGET RA: {}, DEC: {}\nGCAL RA: {}, DEC: {}'.format(ra_src, dec_src, ra_gcal, dec_gcal)
                else:
                    outfile.write( '... TARGET-GCAL separation: {} deg (TARGET: {}, GCAL: {})\n'.format(seprtn,src,gcal) )
            except ValueError as ve:
                outfile.write( '!!!! {}\n\n'.format(ve))
                ra_src0  = src[1:5]; dec_src0 = src[5:]
                ra_gcal0 = gcal[1:5]; dec_gcal0 = gcal[5:]
                ra_src1  = str( ra_src0[0:2]+'h'+ra_src0[2:4]+'m'+'00s' )
                dec_src1 = str( dec_src0[0:3]+'d'+dec_src0[3:5]+'m'+'00s' )
                ra_gcal1 = str( ra_gcal0[0:2]+'h'+ra_gcal0[2:4]+'m'+'00s' )
                dec_gcal1= str( dec_gcal0[0:3]+'d'+dec_gcal0[3:5]+'m''00s' )
                
                D_src   = sk(ra=ra_src1,dec=dec_src1,frame='icrs')
                D_gcal  = sk(ra=ra_gcal1,dec=dec_gcal1,frame='icrs')
                
                ra_src  = D_src.ra.deg; dec_src = D_src.dec.deg
                ra_gcal = D_gcal.ra.deg; dec_gcal = D_gcal.dec.deg
                seprtn  = separation(ra_src, dec_src, ra_gcal, dec_gcal)
                RAs.append(ra_src), DECs.append(dec_src), Dr.append(seprtn)
                if abs(seprtn) <= 16. and src != gcal:
                    outfile.write( '*** TARGET-GCAL separation: {} deg (TARGET: {}, GCAL: {})\n'.format(seprtn,src,gcal) )
                    print(('TARGET-GCAL separation: {} deg (TARGET: {}, GCAL: {})'.format(seprtn,src,gcal)))
                    out_TARGET.append(src)
                    #print 'TARGET RA: {}, DEC: {}\nGCAL RA: {}, DEC: {}'.format(ra_src, dec_src, ra_gcal, dec_gcal)
                else:
                    outfile.write( '... TARGET-GCAL separation: {} deg (TARGET: {}, GCAL: {})\n'.format(seprtn,src,gcal) )
                pass
    outfile.close()
    Dr = np.array(Dr)
    out_TARGET.append(GCALS_in[0])
    RAs, DECs = np.array(RAs), np.array(DECs)
    #plt.plot(RAs, DECs, 'k*')
    #plt.xlabel('RA / deg'), plt.ylabel('DEC / deg')
    #plt.grid(True,ls="-.", alpha=0.35)
    #plt.savefig(output_dir+'calibrators.png')
    #plt.close()
    return (Dr,out_TARGET)

def get_field_names(MS_in):
    # this needs to run in a CASA session
    import numpy as np
    msmd.open(MS_in)
    fieldnames = msmd.fieldnames()
    msmd.done()
    fx = {}
    for i in range(len(fieldnames)):
        fx[fieldnames[i]] = '%d'%i
    return( fieldnames, fx )

def get_bpcal(fieldnames, Field_dict):
    # explicitly looks for and chooses known cal as primary cal
    
    '''
    J0137+3309 = 3C48
    J0521+1638 = 3C138
    J1331+3030 = 3C286
    '''
    
    bp_cal_list = ["J1939-6342", "J193925-634243","J0137+3309",\
        "J0521+1638", "J1331+3030","J0408-6545", "J0252-7104"]
    
    for fld in fieldnames:
        for cal in bp_cal_list:
            if (fld == cal):
                bp_cal0 = fld
    #print " (!) Primary cal: ", bp_cal0
    bp_calibrator = [bp_cal0]
    bpcal_ID = Field_dict[bp_cal0]
    return (bp_calibrator,bpcal_ID)

def get_field_ID(field_in, Field_dict):
    field_ID = Field_dict[field_in]
    return field_ID

def user_defnd_inputs(prefix0):
    """
    function to Associate target name with field ID in the .MS file and also sets the reference antenna.
    """
    if prefix0.endswith('.ms'):
        ms0    = prefix0
        suffix = prefix0[:-10]
    else:
        suffix = '_sdp_l0_1284.full_pol.ms'             #full pol. data
        ms0    = pwd+'msdir/'+prefix0+suffix
    
    #TARGETS,FDB = get_field_names(ms0)

    msmd.open(ms0)
    TARGETS = msmd.fieldnames()
    FDB     = {ti: str( msmd.fieldsforname(ti)[0] ) for ti in TARGETS}
    bpcalID = msmd.fieldsforintent('CALIBRATE_BANDPASS')[0]
    msmd.done()
    
    if prefix0     == "1532725253":
        REF_ANT = "m033"

    elif prefix0   == "1556645454":
        REF_ANT = "m059"
    
    elif prefix0   == "1555162676":
        REF_ANT = "m043"
    
    elif prefix0   == "1541456887":
        REF_ANT = "m011"
    
    elif prefix0   == "1541095095":
        REF_ANT = "m040"
    
    elif prefix0   == "1541522024":
        REF_ANT = "m041"
    
    elif prefix0   == "1544842869":
        REF_ANT = 'm029'
        
    elif prefix0   == "1541107285":
        REF_ANT = "m010"
    
    elif prefix0   == "1543346311":
        REF_ANT = "m045"
    
    elif prefix0   == "1539238540":
        REF_ANT = "m029"
        
    elif prefix0   == "1538697656":
        REF_ANT = "m022"
        
    elif prefix0   == "1538924761":
        REF_ANT = "m049"

    elif prefix0   == "1538843813":
        REF_ANT = "m010" #run the function "get_refant(ms0,BCAL_ID)[0]" within casa to get this

    elif prefix0   == "1541558432":
        REF_ANT = "m039"

    elif prefix0   == "1527564238":
        REF_ANT = "m033"
    
    else:
        REF_ANT = get_refant(ms0,bpcalID)[0]
    
    return(TARGETS,FDB,REF_ANT)


def cal_inputs(gcal_id, obs):
    pwd        = "/scratch2/slegodi/raw_obs/Mcals/"
    prefix0    = obs[0]                                         # 1st element in 'obs' must be the observation ID
    specres    = obs[1]                                         # 2nd element must be the observation mode; e.g '4K'
    ms0        = pwd+'msdir/%s_sdp_l0.full_1284.full_pol.ms'%prefix0
    stim_code  = pwd+'my_stimela_script.py'
    in_dir     = pwd+'msdir'
    MSNAME     = '%s_sdp_l0.full_1284.full_pol.ms'%prefix0
    gc_0_split = check_directory(prefix0 + ".0GC.ms")
    gc_1_split = check_directory(prefix0 + ".1GC.ms")
    
    TARGETS,FDB,ref_ant = user_defnd_inputs(prefix0)
    #GCALS: gain calibrator (to be looped through the "TARGET" list of sources with 16 deg of each other)
    #       run the function "get_trgt_gcal_sep(TARGETS, GCALS[0], out_dir)" to get groups of these "TARGET" sources
    
    cal_radius = 16.            #max distance between targets and gain cal. in degrees
    GCALS0     = TARGETS
    GCALS      = [[GCALS0[gcal_id]]] # cycle through the TARGETS list to switch between gain calibrators via 'gcal_id'
    TARGET     = get_trgt_gcal_sep(TARGETS, GCALS, None)[1] 
    GCAL_ID    = get_field_ID(GCALS[0][0], FDB)
    BPCAL      = get_bpcal(TARGETS,FDB)[0]
    BCAL_ID    = get_bpcal(TARGETS,FDB)[1]
    out_dir    = pwd+'%s-%s-stimela-output-%s'%(specres,prefix0,GCALS[0][0])    #check_directory(pwd+'stimela-output-%s'%GCALS[0][0])
    primcal_ms = check_directory("primcal-1GC-%s.ms"%GCALS[0][0])
    primcal_im = check_directory("primcal-%s"%GCALS[0][0])
    
    return (obs,pwd,prefix0,specres,ms0,stim_code,in_dir,MSNAME,gc_0_split,gc_1_split,TARGETS,FDB,ref_ant,\
        cal_radius,GCALS0,GCALS,TARGET,GCAL_ID,BPCAL,BCAL_ID,out_dir,primcal_ms,primcal_im)


def find_sources(gcal_in, fits_in,base, pwd, show_fitted):
    # PYBDSF spurce finding in python (not via stimela!)
    #pwd            = "/scratch2/slegodi/raw_obs/Mcals/"
    #inputdir       = pwd+"sourcefind_images/test-images/"
    inputdir       = pwd #pwd+"sourcefind_images/"
    fitsfile       = fits_in #files0[5]   # will need to loop through list of these files!
    if base == None:
        base           = fitsfile[3:].rstrip('-1GC-.fits')
    n_cores        = 16
    snrthreshold   = 10.0                # SNR limit for good sources
    isl_thresh     = 10.0                # SNR limit for source finding island region
    pix_thresh     = 1.5*isl_thresh      # SNR limit for source detection/identification
    fluxthreshold  = 1e-04
    pblimit        = 0.05  #this should be constant at about 5%, hence 0.05 
    
    if gcal_in == None:
        gcal_in = base
        #outdir  = 'SNR-thresh_isl-10-'+str(snrthreshold)+'.out.'+base+'/'
    else:
        base    = fitsfile.rstrip('-1GC-.fits')
        #outdir  = base+'/'
    outdir  = pwd+base+'pybdsf.results'
    if not (os.path.isdir(outdir)):
        os.mkdir(outdir)
        
    #try:
        #if not (os.path.isdir(outdir)):
            #os.mkdir(outdir)
    #except OSError as ose2:
        #print ose2,"\n\n"
        #base0  = os.path.basename(fitsfile)
        #base   = base0.rstrip('-image.fits')
        #outdir = base+"-PYBDSF"
        #os.mkdir(outdir)
    
    print((' > OUTPUT dir: {}'.format(outdir)))

    try:
        img = bdsf.process_image(fitsfile,indir=inputdir, rms_box=None, output_opts=True,output_all=True,\
            shapelet_do=True,solnname=outdir+base,thresh='hard',thresh_isl=isl_thresh, thresh_pix=pix_thresh,\
                psf_vary_do=False,ncores=n_cores, rms_value=fluxthreshold)
    except RuntimeError as re:
        print(re)
        inputdir = pwd+"sourcefind_images/test-images/"
        #os.system("cp *.fits "+inputdir+"*.fits .")
        fitsfile = fits_in
        img      = bdsf.process_image(fitsfile, rms_box=None, output_opts=True,output_all=True,\
                        shapelet_do=True,solnname=outdir+base,thresh='hard',thresh_isl=isl_thresh, thresh_pix=pix_thresh,\
                            psf_vary_do=False,ncores=n_cores, rms_value=fluxthreshold)

    meanimage   = outdir+base+'-meanimage.fits'
    chan0image  = outdir+base+'-chan0image.fits'
    rmsfile     = outdir+base+'-rms.fits'
    resfile     = outdir+base+'-residuals.fits'
    #polimage   = outdir+'base+-polimage.fits'
    modimage    = outdir+base+'-modimage.fits'
    srcfilename = outdir+base+'-source-cat.txt'

    img.write_catalog(outfile = srcfilename,clobber =True,format='ascii',catalog_type ='gaul')
    #img.export_image(outfile =resfile,img_format = 'fits',img_type = 'gaus_resid',clobber =True)
    img.export_image(outfile =rmsfile,img_format = 'fits',img_type = 'rms',clobber =True)
    #img.export_image(outfile =meanimage,img_format = 'fits',img_type = 'mean',clobber =True)
    #img.export_image(outfile =chan0image,img_format = 'fits',img_type = 'ch0',clobber =True)
    img.export_image(outfile =modimage,img_format = 'fits',img_type = 'gaus_model',clobber =True)

    #img.export_image(outfile =polimage,img_format = 'fits',img_type = 'pi',clobber =True)
    #img.export_image(outfile =psf_ratio_image,img_format = 'fits',img_type = 'psf_ratio',clobber =True)

    print('\n\n done exporting bdsf.process_image result files ... Making catalog.')

    hdulist_I      = fits.open(fitsfile)
    w              = WCS(hdulist_I[0].header)
    Idata          = hdulist_I[0].data
    Imax           = Idata.max()
    ###
    table          = Table.read(srcfilename,format = 'ascii')
    src_id         = table['col1']
    ra             = table['col3']#+360   # don't remember why I had to add 360 deg??
    raerr          = table['col4']
    dec            = table['col5']
    decerr         = table['col6']
    peak_flux      = table['col9']
    peak_flux_err  = table['col10']
    total_flux     = table['col7']
    total_flux_err = table['col8']
    Smaj_ax        = table['col15']*3600  # in arcseconds
    err_Smaj_ax    = table['col16']*3600
    Smin_ax        = table['col17']*3600
    err_Smin_ax    = table['col18']*3600
    PA             = table['col19']+90.0
    PA_err         = table['col20']
    RMS_col        = table['col43']
    #print 

    print('done inporting source list attributes ...')
    bins = int(0.25*(sqrt(len(peak_flux)))*2)

    # filter the sources based on thresholds
    Imap= []; Ifit=[]; Qmap=[]; Umap=[]; Pmap = []; good = []#weight = []
    goodsources=0
    sources = []
    infield = base

    for i in range(len(ra)): 
        try:
            print((' > ra[%d]: %.4f, dec[%d]: %.4f'%(i,ra[i],i,dec[i])))
            pixcrd = w.wcs_world2pix(ra[i],dec[i],0,0,0)
            y      = pixcrd[0]
            x      = pixcrd[1]
            print((' > xpix[%d]: %.4f, ypix[%d]: %.4f'%(i,x,i,y)))
            yp     = int(round(y))
            xp     = int(round(x))
            im_pix = 10800
            rms    = RMS_col[i]
            if (xp < im_pix) and (yp < im_pix): 
                try:
                    Iamp = Idata[xp,yp]
                except IndexError as ie:
                    print(ie)
                    try:
                        Iamp = Idata[0,0,xp,yp]
                    except IndexError as ie02:
                        print((' - (!) ', ie02))
                        Iamp = rms
                        pass
                wdata = 1. #pbcorr(xp,yp,cell_deg,bmaj)
                fluxcorr = safe_div(1.,wdata)
                snr = peak_flux[i]/rms
                print(('***SNR: %4.2f, rms: %4.2e, Iamp: %4.2e' %(snr, rms, Iamp)))
                if(snr > snrthreshold) and (Iamp > fluxthreshold):# and (fluxcorr < fluxcorrlimit):       
                    good.append(True) 
                    goodsources = goodsources + 1
                    print(("source: %5d %9.4f %7.4f %7.1f %7.1f %7.3f %9.3f %8.3f %7.2f %7.2f" \
                    % (goodsources,ra[i],dec[i],y,x,1000*rms,1000*total_flux[i],1000*peak_flux[i],\
                        Smaj_ax[i],Smin_ax[i])))
                else: 
                    good.append(False)
                Imap.append(Iamp)
                Ifit.append(peak_flux[i])
                ra_err = 3600.0*raerr[i]
                dec_err = 3600.0*decerr[i] 
                f = 1000*total_flux[i]*fluxcorr
                p = 1000*peak_flux[i]*fluxcorr
                dp = 1000*peak_flux_err[i]*fluxcorr
                df = 1000*total_flux_err[i]*fluxcorr
                if(total_flux_err[i] < 0.0) or (isnan(dp)): # replace errors with calculated values if AEGEAN fit error fails
                    (df,dp,ra_err,dec_err)=source_errors(Smaj_ax[i],Smin_ax[i],PA[i],total_flux[i],\
                                                        peak_flux[i],rms,beam)
                    ra_err = 3600.0*ra_err
                    dec_err = 3600.0*dec_err
                    dp = 1000*dp*fluxcorr
                    df = 1000*df*fluxcorr
                if(ra_err == -3600.0):
                    ra_err = -1.0
                if(dec_err == -3600.0):
                    dec_err = -1.0  
                fluxcorr_mJy = 1000*rms*fluxcorr
                sources.append(source(ra[i],ra_err,dec[i],dec_err,f,df,p,dp,\
                            Smaj_ax[i],err_Smaj_ax[i],Smin_ax[i],err_Smin_ax[i],PA[i],PA_err[i],\
                            fluxcorr_mJy,y,x))
                src_line='{:3.2f}&{:3.2f}&{:3.2f}&{:3.2f}&{:1.2e}&{:1.2e}&{:1.2e}&{:1.2e}&{:1.2e}&{:1.2e}&{:1.2e}&{:1.2e}&{:3.2f}&{:3.2f}&{:1.2e}&{:d}&{:d}\\\ \n\
                '.format(ra[i],ra_err,dec[i],dec_err,f,df,p,dp,Smaj_ax[i],err_Smaj_ax[i],Smin_ax[i],err_Smin_ax[i],PA[i],PA_err[i],fluxcorr_mJy,int(y),int(x))
        except Exception as ex002:
            print((' (!) > ', ex002))
            pass
        
    print(("\n\nFound %d good sources" % goodsources))
    print(("Rejected %d sources with map peak < %7.2e and snr < %4.1f, and pblimit < %4.2f " % \
        ((len(sources) - goodsources),fluxthreshold,snrthreshold,pblimit)))

    # open output files and write headers
    goodregfilename = outdir+'sources_pybdsm.reg'        # ds9 reg file of good sources
    badregfilename =  outdir+'badsources_pybdsm.reg'     # ds9 reg file of bad sources
    sourcefilename =  outdir+'sources_pybdsm.txt'        # output file for good sources
    sourcefilename_neg =  outdir+'neg_sources_pybdsm.txt'        # output file for good sources

    regfile = open(goodregfilename,'w')
    badfile = open(badregfilename,'w')
    sourcefile = open(sourcefilename,'w')
    neg_sourcefile_pybdsm = open(sourcefilename_neg,'w')

    regfile.write('# Region file format:  DS9 version 4.1 \n')
    regfile.write('global color=green width=1 \n')
    regfile.write('fk5\n')
    badfile.write('# Region file format:  DS9 version 4.1 \n')
    badfile.write('global color=red width=1 \n')
    badfile.write('fk5\n')
    sourcefile.write(  '#   id     ra    ra_err    dec   dec_err  i_flux   i_err   p_flux  p_err     rms      a   a_err      b   b_err     pa   pa_err    x       y')
    sourcefile.write('\n#        (deg)   (")     (deg)    (")         (mJy)           (mJy)        (mJy)       (")           (")           (deg)           (pix)')

    neg_sourcefile_pybdsm.write(  '#   id     ra    ra_err    dec   dec_err  i_flux   i_err   p_flux  p_err     rms      a   a_err      b   b_err     pa   pa_err    x       y')
    neg_sourcefile_pybdsm.write('\n#        (deg)   (")     (deg)    (")         (mJy)           (mJy)        (mJy)       (")           (")           (deg)           (pix)')

    print('---------------- writing to ds9 region files and output source file ---------------')

    if len(sources) > 0:
        quick(sources,0,len(sources)-1)     # sort the source list by right ascension
        j=0
        for i in range(len(sources)):
            if(good[i]):
                j=j+1
                regline = "ellipse(%9.4f,%9.4f,%7.2f\",%7.2f\",%6.1f) \n" %\
                (sources[i].ra,sources[i].dec,sources[i].Smaj,sources[i].Smin,sources[i].Spa)
                regfile.write(regline)
                sourceline = "\n%5d  %8.5f %5.2f  %8.5f %5.2f %8.4f %7.4f %8.4f %8.4f %7.4f %7.2f %6.2f %7.2f %6.2f %7.1f %5.1f %8.2f %8.2f" % \
                (j,sources[i].ra,sources[i].dra,sources[i].dec,sources[i].ddec,sources[i].flux,sources[i].df,sources[i].peak,sources[i].dp, \
                sources[i].rms,sources[i].Smaj,sources[i].dmaj,sources[i].Smin,sources[i].dmin,sources[i].Spa,sources[i].dpa,\
                sources[i].x,sources[i].y)
                sourcefile.write(sourceline)
            else:
                regline = "ellipse(%9.4f,%9.4f,%7.2f\",%7.2f\",%6.1f) \n" %\
                (sources[i].ra,sources[i].dec,sources[i].Smaj,sources[i].Smin,sources[i].Spa)
                badfile.write(regline)
        print(("wrote %d sources to file %s" % (j,sourcefilename)))
        regfile.close()
        badfile.close()
        sourcefile.close()

        t_pybdsm = Table.read(sourcefilename,format='ascii')

        try:
            iflux_pybdsm = t_pybdsm['i_flux']*1000.       # converts to mJy
            pflux_pybdsm = t_pybdsm['p_flux']*1000.
            ibins  = int(sqrt( len(iflux_pybdsm) ))*2
            pbins  = int(sqrt( len(pflux_pybdsm) ))*2
            
            mpl.hist(iflux_pybdsm,ibins,histtype='step', ec='gray',\
                    label = r'$\rm{S_{int}\, ,all\,sources\,PYBDSM}$',linewidth = 2.5 )
            mpl.hist(pflux_pybdsm,pbins,histtype='step', ec='blue',\
                    label = r'$\rm{P_{int}\, ,all\,sources\,PYBDSM}$',linewidth = 2.5 )
            mpl.xscale("log")
            mpl.xlim(10**1.,10**5)
            #mpl.ylim(0.,600.)
            mpl.xlabel('Flux / [mJy]')
            mpl.legend(prop={'size':10},loc='best',frameon = False,scatterpoints = 1)
            mpl.savefig(outdir+'histo_pybdsmaegean.pdf')
            #mpl.show()
        except:
            print(('(!) .. I flux may be infinite again ... check it in  {}'.format(sourcefilename)))
            pass
        
        pk_src = nanargmax(peak_flux)
        figure(figsize=(20,6))
        plot(src_id, peak_flux), plot(src_id, peak_flux, 'wo')
        plot(src_id[pk_src], peak_flux[pk_src], 'r*',\
            label=' src_id: %d (%s)\n Ipeak: %1.4e Jy/bm\n RA     : %2.4f deg\n Dec   : %2.4f deg\
        '%(src_id[pk_src],gcal_in, peak_flux[pk_src], ra[pk_src], dec[pk_src]) )
        legend(loc=0, fancybox=True, framealpha=0.5)
        xlabel('PyBDSM Source ID')
        ylabel('Peak Flux / [Jy/beam]')
        grid(ls=':', alpha=0.4)
        savefig(outdir+'srcs_vs_peakFlux.png')
        #show()
        
        #figure(figsize=(20,6))
        #subplot(131)
        #hist(peak_flux, bins, alpha=0.35), xlabel('Peak Flux')
        #xscale("log")
        #subplot(132)
        #scatter(ra, peak_flux)
        #xlabel('RA/deg'), ylabel('Peak Flux')
        #subplot(133)
        #scatter(dec, peak_flux)
        #xlabel('Dec/deg'), ylabel('Peak Flux')
        #tight_layout()
        #show()

    else:
        print('(!) *** No sources found! *** ')
    
    if show_fitted:
        img.show_fit(source_seds=True, gresid_image=False, mean_image=False,\
                    rms_image=True, ch0_flagged=True, ch0_image=True, smodel_image=True)
        show()
    
    return ()

def pk_attributes(spec_res,wdir):
    # for analysing sourcefinding calibrator results 
    cal_names, RAs, DECs, obs_mode, flx_pk, eflx_pk, flx_int, eflx_int = [],[],[],[],[],[],[],[]
    eRA, eDEC = [],[]
    for catsdir in os.listdir(wdir):
        if catsdir.startswith('SNR-thresh_isl-10-'+spec_res+'-1GC-'):# and catsdir.endswith("_noisyIQU.tbl"):
            catsdir = wdir+catsdir+'/'
            for src_cat in os.listdir(catsdir):
                if src_cat.startswith(spec_res+'-1GC-') and src_cat.endswith('-source-cat.txt'):
                    cat            = Table.read(catsdir+src_cat, format='ascii')
                    src            = src_cat.lstrip(spec_res+'-1GC-')
                    src            = src.rstrip('-source-cat.txt')
                    gcal           = src 
                    print(('Gain cal: {} source catalog found ..'.format(src)))
                    
                    src_id         = cat['col1']
                    Isl_id         = cat['col2']
                    ra             = cat['col3']   # don't remember why I had to add 360 deg??
                    raerr          = cat['col4']
                    dec            = cat['col5']
                    decerr         = cat['col6']
                    peak_flux      = cat['col9']
                    pk_flx2        = cat['col9']
                    peak_flux_err  = cat['col10']
                    total_flux     = cat['col7']
                    total_flux_err = cat['col8']
                    Smaj_ax        = cat['col15']*3600  # in arcseconds
                    err_Smaj_ax    = cat['col16']*3600
                    Smin_ax        = cat['col17']*3600
                    err_Smin_ax    = cat['col18']*3600
                    PA             = cat['col19']+90.0
                    PA_err         = cat['col20']
                    RMS_col        = cat['col43']
                    
                    pk_src         = nanargmax(peak_flux)
                    s_i            = pk_src
                    
                    cal_names.append(gcal), DECs.append(dec[s_i]), obs_mode.append(spec_res), flx_pk.append(peak_flux[s_i])
                    eflx_pk.append(peak_flux_err[s_i]), flx_int.append(total_flux[s_i]), eflx_int.append(total_flux_err[s_i])
                    eRA.append(raerr[s_i]), eDEC.append(decerr[s_i])
                    
                    if ra[s_i] > 0.:
                        RAs.append(ra[s_i])
                    else:
                        RAs.append(ra[s_i]+360.0)
                    
                    figure(figsize=(10,8))
                    plot(src_id, peak_flux, 'k:', alpha=0.3, label='_nolegend_')
                    plot(src_id, peak_flux, 'r.')
                    plot(src_id[s_i], peak_flux[s_i], '*', markersize=15,\
                        label=' src_id: %d\n Ipeak: %1.4f Jy/bm\n RA     : %2.4f deg\n Dec   : %2.4f deg\
                    '%(src_id[s_i], peak_flux[s_i], ra[s_i], dec[s_i]) )
                    
                    pk_flx2[pk_src]= 0.                 # removes brightest source to find 2nd brightest source easier
                    pk_src2        = nanargmax(pk_flx2) # finds 2nd brightest source
                    s_i            = pk_src2
                    
                    plot(src_id[s_i], peak_flux[s_i], '*', markersize=15,\
                        label=' src_id: %d\n Ipeak: %1.4f Jy/bm\n RA     : %2.4f deg\n Dec   : %2.4f deg\
                    '%(src_id[s_i], peak_flux[s_i], ra[s_i], dec[s_i]) )
                    title('Gain cal: {} ({} data)'.format(gcal, spec_res))
                    legend(loc=0, fancybox=True, framealpha=0.5)
                    xlabel('PyBDSM Source ID')
                    ylabel('Peak Flux / [Jy/beam]')
                    grid(ls=':', alpha=0.5)
                    try:
                        figname = catsdir+'Gaincal.%s.peakFlux.png'%gcal
                        savefig(figname)
                        print(('  >> Figure file: {} saved!!'.format(figname)))
                    except ValueError as voe:
                        print(voe)
                        print(('  (!) Figure file: {} NOT saved!!'.format(figname)))
                        pass
                    #show()
                    close()
    
    tt = atpy.Table(); tt_name = wdir+spec_res+'-gcal.catalog.tbl'
    tt.add_column('Source', cal_names)
    tt.add_column('RA', RAs, unit='deg', dtype=float32)
    tt.add_column('DEC', DECs, unit='deg', dtype=float32)
    tt.add_column('errRA', eRA, unit='deg', dtype=float32)
    tt.add_column('errDEC', eDEC, unit='deg', dtype=float32)
    tt.add_column('Obs.mode', spec_res)
    tt.add_column('I_peak', flx_pk, unit='Jy/beam', dtype=float32)
    tt.add_column('errI_peak', eflx_pk, unit='Jy/beam', dtype=float32)
    tt.add_column('I_int', flx_int, unit='Jy/beam', dtype=float32)
    tt.add_column('errI_int', eflx_int, unit='Jy/beam', dtype=float32)
    write_myipac(tt_name, tt)
    
    return()

def Area_ellipse(a,b):
    return (pi*a*b)

def substrctr_fts(fits_catalogue, perc_cut):
    # fits_catalogue is a **PYBDSF** produced Gaussian sky model components
    if perc_cut == None:
        perc_cut = 0.01
    tfits    = fits_catalogue
    t0       = atpy.Table()
    t0.read(tfits)
    #t0.describe()
    
    srcID,Isl_ID   = t0['Source_id'], t0['Isl_id']
    maj_ax, min_ax = t0['Maj_img_plane'], t0['Min_img_plane']   # in degrees
    RAs, DECs      = t0['RA'], t0['DEC']                        # 1-sigma errors in degrees
    eRAs, eDECs    = t0['E_RA'], t0['E_DEC']                    # in degrees
    Iint, Ipk      = t0['Total_flux'], t0['Peak_flux']          # in Jy/bm
    eIint, eIpk    = t0['E_Total_flux'], t0['E_Peak_flux']      # 1-sigma errors in Jy/bm
    Ipeak          = nanmax(Ipk)                                # in Jy/bm
    a, b           = maj_ax[argmax(Ipk)], min_ax[argmax(Ipk)]
    
    print((' === Secondary sources in the field. Percent cutoff: {}%, Ipk: {:.3} Jy/bm ==='.format(perc_cut*100.,Ipeak)))
    print('Source ID , Sep_from_peak[\'], RA[deg], Dec[deg], Ipk_Fraction[%%]')
    sep_Ipk = []; frc_Ipk = []; sub_src = []
    for i in range(len(srcID)):
        if Ipk[i] > perc_cut*Ipeak:
            di_pk     = separation(RAs[ nanargmax(Ipk) ], DECs[ nanargmax(Ipk) ],RAs[i], DECs[i])
            di_pk_min = di_pk*60.; frc_Ipk_i = 100.*(Ipk[i]/Ipeak)
            sep_Ipk.append(di_pk_min); frc_Ipk.append(frc_Ipk_i)
            sub_src.append([ srcID[i], RAs[i], DECs[i], Area_ellipse(maj_ax[i], min_ax[i]) ])
            print((' {} , {:.2f}  , {:.6f} , {:.6f} , {:.2f}%%'.format(srcID[i], di_pk_min, RAs[i], DECs[i], frc_Ipk_i )))
    return (array(sep_Ipk), array(frc_Ipk), array(sub_src))

#################################################################################
if __name__ == '__main__':
    run_this_code = True

    if run_this_code:
        do_pybdsf = False
        do_1k     = False           # perform sourcefinding on 1K obs
        do_4k     = True           # perform sourcefinding on 4K obs
        showfits  = False           # show plotted results
        get_attrs = False           # turn to get source attributes
        do_skymod = True
        do_prim   = False
        do_prerred_bpcals = False
        do_test   = False
        substrctr_txt = False
        plot_skymods = False

    test_dir = '/home/samuel/SARAO/MeerCals/sky_models/all_skymodels/'
    tfits0   = test_dir+'smooth.bcal-1543346311-J0408-6545.pybdsm.gaul.FITS'
    cats_dir = '/home/samuel/SARAO/MeerCals/sky_models/all_skymodels/'

    if plot_skymods:
        for f0 in os.listdir(cats_dir):
            if f0.endswith('.gaul.FITS'):
                targ_name = f0[:-17]
                print(('\n>> TARGET: {}'.format(targ_name)))
                ds, fs, spos = substrctr_fts(cats_dir+f0, None)
                plt.loglog(ds, fs, '.'), plt.xlabel('Separation from peak source [arcmin]')
                plt.ylabel('I$_{peak}$ fraction [%]')
                plt.grid(True,ls="-.", alpha=0.35)
                plt.show()

    if do_test:
        pwd       = '/scratch2/slegodi/raw_obs/self_calib/1532725253-5-phase.selfcal/DEEP_2.1532725253/DEEP_2_scal/'
        #fits_file = pwd+'DEEP_2_mfs.sc4-image.fits'
        fits_file = 'DEEP_2_mfs.sc4-image.fits'
        base      = 'DEEP_2_mfs.sc4-'
        cwd       = os.getcwd()
        os.chdir( pwd )
        print(('  **** Changed from directory\n {} to directory\n {}'.format(cwd, pwd)))
        find_sources(None, fits_file, base, pwd, show_fitted = showfits)
        os.chdir( cwd )
        print(('  **** Changed BACK to directory: {}'.format(cwd)))

    if get_attrs:
        obs_res  = ['4K', '1K']
        cwd      = "/home/samuel/SARAO/MeerCals/final.results/"     # parent working directory on xps15
        pk_attributes(obs_res[0],cwd)

    if do_prerred_bpcals:
        # - run Pybdsf for a list of preffered bcals (initially J1939-6342, J0408-6545, and J0252-7104)
        # - meant for testing and making local sky models of the above bpcals for purposes of flux calibration
        # - (!) images in "imdir" are selfcalibrated with the script "cals_selfcal.py"
        import bdsf, katdal
        pwd0  = "/scratch2/slegodi/raw_obs/Mcals/" 
        imdir = pwd0+'sourcefind_images/'
        #imdir = '/scratch2/slegodi/raw_obs/Mcals/msdir/preffered_bcals/tests/sfcal_images/'
        os.chdir( imdir )
        cwd = os.getcwd()
        print('\n============')
        print((' PWD: ', cwd))
        print('============')
        for fitsfile in os.listdir(cwd):
            if fitsfile.startswith('smooth.bcal-1543346311-') and fitsfile.endswith('.fits'):
                #bpcal = fitsfile.rstrip("_mfs.sc2.image.tt0.fits")
                #fitsfile = cwd+'/'+fitsfile
                bpcal    = fitsfile.rstrip(".fits")
                bpcal    = fitsfile.lstrip("smooth.bcal-1543346311-")
                find_sources(bpcal, fitsfile, bpcal,cwd, show_fitted = showfits)

    if do_pybdsf:
        ##################################
        # - runs on my python environment (new_env1) on com4 where Pybdsf, katdal are installed
        # - currently works best in the location of the source finding images
        # - so cd into that first:
        import bdsf, katdal
        os.chdir( '/scratch2/slegodi/raw_obs/Mcals/sourcefind_images/' )
        cwd = os.getcwd()
        print('\n============')
        print((' PWD: ', cwd))
        print('============')
        pwd0   = "/scratch2/slegodi/raw_obs/Mcals/"               # parent working directory on com4
        files0 = ["v6.1K-1GC-J1120+1420.fits", "v6.1K-1GC-J0240-2309.fits", "v6.4K-1GC-J1911-2006.fits",\
            "v6.4K-1GC-J0252-7104.fits", "v6.4K-1GC-J1744-5144.fits", "v6.4K-1GC-J1830-3602.fits"]
        ##################################
        
        if do_4k:
            observatns = {'15388': ["1538843813", "4K", [3,15]], '15275': ["1527564238", "4K", [0,1,6,14]], '15389': ["1538924761", "4K", [13,17,40]],\
            '15386': ["1538697656", "4K", [6,12,19,28,39,47]], '15392': ["1539238540", "1K", [0,2,16,26,39,50,54]],\
                '15433': ["1543346311", "4K", [1,2,4,16,24,25,28]], '15411': ["1541107285", "4K", [1,2,4,5,7,16]],\
                    '15448': ["1544842869", "4K", [0,1,3,8,17,30]], '15415': ["1541522024", "4K", [1,3,10]],\
                        '15410': ["1541095095", "4K", [1,2,5,24,33]], '15414': ["1541456887", "4K", [1,2,4,5,21]] }
        
            obs_in     = observatns['15433']
            g_id_start = 33
            g_id_stop  = 36
            
            slct_gcals = obs_in[2]
        
        #for g_id in np.arange(g_id_start,g_id_stop): # for iterating through the whole set of gaincals
        for g_id in slct_gcals:
            try:
                obs_ID,pwd,prefix0,specres,ms0,stim_code,in_dir,MSNAME,gc_0_split,gc_1_split,TARGETS,FDB,ref_ant,\
                    cal_radius,GCALS0,GCALS,TARGET,GCAL_ID,BPCAL,BCAL_ID,out_dir,primcal_ms,primcal_im = cal_inputs(g_id,obs_in)
                
                gcal    = GCALS[0][0]
                srcdir  = pwd+'sourcefind_images/'
                filestr = '/v6.all.fields.wsclean-1GC-'
                fileext = '-GCAL-image.fits '   # note the space after '.fits'
                
                if do_1k:
                    cp_comand = 'cp '+out_dir+filestr+gcal+fileext+srcdir+'v6.1K-1GC-'+gcal+'.fits'
                    src_find_file = "v6.1K-1GC-%s.fits"%gcal
                    error_log     = open(srcdir+'1K-src_find_errors.log','a',0)
                    
                if do_4k: 
                    if not do_skymod:
                        cp_comand = 'cp '+out_dir+filestr+gcal+fileext+srcdir+'v6.4K-1GC-'+obs_in[0]+'-'+gcal+'.fits'
                        src_find_file = "v6.4K-1GC-%s-%s.fits"%(obs_in[0],gcal)
                        error_log     = open(srcdir+'4K-src_find_errors.log','a',0)
                    
                    elif do_skymod:
                        filestr       = 'primcal-'
                        fileext       = '-image.fits '# note the space after '.fits'
                        
                        cp_comand     = 'cp '+out_dir+'/'+filestr+gcal+fileext+srcdir+'bcal-'+obs_in[0]+'-'+gcal+'.fits'
                        src_find_file = "bcal-%s-%s.fits"%(obs_in[0],gcal)
                        error_log     = open(srcdir+'bcal-src_find_errors.log','a',0)
                    
                print((' > Running: {}'.format(cp_comand)))
                os.system( cp_comand )
                
                find_sources(None,src_find_file,gcal, pwd, show_fitted = showfits)
            except Exception as ex1:
                print(('\n (!): ', ex1))
                error_log.write( '\n (!): {}'.format(ex1) )
                raise
                pass



    if substrctr_txt:
        test_dir  = '/home/samuel/MeerCals/sky_models/smooth.bcal-1543346311-J0408-6545/'
        txt_file  = test_dir+'smooth.bcal-1543346311-J0408-6545-source-cat.txt'
        txt       = open(txt_file, 'r')
        RAs, DECs = [], []; eRAs, eDECs = [], []
        Iint, Ipk = [], []; eIint, eIpk = [], []
        datalines = txt.readlines()
        for i in range (len(datalines)):
            l_i = datalines[i]
            if (l_i[0] !='#'):
                try:
                    RAs.append( float(l_i.split()[4]) ); eRAs.append( float(l_i.split()[5]) )
                    DECs.append( float(l_i.split()[6]) ); eDECs.append( float(l_i.split()[7]) )
                    Iint.append( float(l_i.split()[8]) ); eIint.append( float(l_i.split()[9]) )
                    Iint.append( float(l_i.split()[10])); eIint.append( float(l_i.split()[11]))
                except Exception as ex003:
                    print((' (!) On line {}: {}'.format(i, ex003)))

