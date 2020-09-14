"""
- Run script in casa session
- does image analysis (via CASA) and source finding (via PYBDSF)
- can be imported into other scripts
"""
########################################################################################
from astropy.io import fits 
import argparse, csv, os, time, sys, atpy
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
import matplotlib.pyplot as plt
from numpy import *

try:
    import pandas as pd  
except ImportError as ie1:
    print (" ! >" , ie1, "\n\n\n")
    pass

try:
    import py3gokatsdpimager as g 
    import image_analysis as a
except ImportError as ie2:
    print (" ! >" , ie2, "\n\n\n")
    print (" ! > run this part of code in python only, NOT casa!"), "\n\n\n"
    pass
########################################################################################
def get_bins(Arr):
    """
    Calculates number of histogram bins given 1d array
    """
    nbins = int(sqrt(len(Arr)))*2
    return (nbins)

def casa_imstats(IMS_list,logsdir=None): 
    """
    # Calculating image stats via CASA imstat()
    # results logged to fitsfile directory
    # input is a list of fitsfile paths or single fits
    """
    algs = ["biweight", "chauvenet","classic", "fit-half", "hinges-fences"]
    #algs = ['fit-half']
    stat_list = []
    if not type(IMS_list) == list:
        IMS_list = [IMS_list]
    
    for algo in algs:
        for im in IMS_list:
            if logsdir==None:
                logsdir = os.path.split(im)[0]+'-IMSTATS/'
                if not (os.path.isdir(logsdir)):
                    os.mkdir(logsdir)
                
            imbase = os.path.basename(im).rstrip('.fits')
            flog   = logsdir+imbase+'-imstat.%s.log'%(algo); i = 1
            while os.path.exists(flog):
                flog = flog[:-3]+'redo.%d.log'%(i)
                if not os.path.exists(flog):
                    break
                i += 1
            s_i    = imstat(imagename=im,algorithm=algo,logfile=flog,center='zero')
            stat_list.append( { str(imbase):['algorithm = %s'%algo,s_i] } )
        
    return(stat_list)

def getAvgStdev(Arr):
    avg, std = nanmean(Arr), nanstd(Arr)
    return (avg, std)

#def get_StokesQUV_atImax(inim, algo=None):
    #from astropy.io import fits
    #from astropy.wcs import WCS
    #import os
    #"""
    #Use CASA imstat() to get Stokes QUV at the position of I max flux.
    #- Input: fitsimage and imstat algorithm (optional)
    #- Output: frequency, max and rms flux.
    #"""
    #if algo == None:
        #algo = 'classic'
    
    #try:
        #if algo == "fit-half":
            #sts = imstat(inim, algorithm=algo,center='zero')
        #else:
            #sts = imstat(inim, algorithm=algo)
    #except Exception as ex001:
        #print(ex001)
        #pass
    
    #dd, hdr = fits.getdata(inim, header=True)
    
    

def get_maxflx(inim, algo=None):
    from astropy.io import fits
    import os
    """
    Use CASA imstat() to get max and rms flux.
    - Input: fitsimage and imstat algorithm (optional)
    - Output: frequency, max and rms flux.
    """
    
    if algo == None:
        algo = 'classic'
    
    try:
        if algo == "fit-half":
            sts = imstat(inim, algorithm=algo,center='zero')
        else:
            sts = imstat(inim, algorithm=algo)
    except Exception as ex001:
        print(ex001)
        pass
    
    dd, hdr = fits.getdata(inim, header=True)
    
    fmx  = sts['max'][0]
    frms = sts['rms'][0]
    freq = hdr['CRVAL3']  # freq in Hz
    src  = hrd['OBJECT']
    print(f'\n\n >> Input fits has source: {src}, reference freq: {freq}\n\n')
    return (freq, fmx, frms)

def chan_imstat(imdir, prefix, suffix):
    """
    Plots image statistics and prints them to file named 'prefix-suffix-imstats.csv'.
    Works best with a list of fits images inside directory'imdir' that have 
    'prefix' and 'suffix' specified.
    """
    
    from matplotlib import pyplot as plt
    import os
    
    freqs, fdmax, fdmin, fdavg, fdrms, fdsgm = [],[],[],[],[],[]
    outfile = open(imdir+'{}-{}-imstats.csv'.format(prefix,suffix),'w+',0)
    outfile.write('#frq, fmx, frms\n')
    for f in os.listdir(imdir):
        if f.startswith(str(prefix)) and f.endswith(str(suffix)):
            ff = imdir+'/'+f
            try:
                frq, fmx, frms = get_maxflx(ff)
                freqs.append(frq), fdmax.append(fmx), fdrms.append(frms)
                outfile.write('{},{},{}\n'.format(frq, fmx, frms))
            except Exception as ex:
                print(" (!) Exception: ", ex)
                pass
    
    plt.figure()
    plt.title(str(prefix))
    plt.scatter(freqs, fdmax)
    plt.grid(ls=':', alpha=0.5)
    plt.xlabel('Frequency / GHz')
    plt.ylabel('Peak Flux / Jy/bm')
    plt.tight_layout()
    plt.savefig(imdir+'{}-{}-fmax.png'.format(prefix,suffix[:-5]))
    
    plt.figure()
    plt.title(str(prefix))
    plt.scatter(freqs, fdrms)
    plt.grid(ls=':', alpha=0.5)
    plt.xlabel('Frequency / GHz')
    plt.ylabel('Flux RMS / Jy/bm')
    plt.tight_layout()
    plt.savefig(imdir+'{}-{}-frms.png'.format(prefix,suffix[:-5]))
    plt.show()
    outfile.close()
    return()

def plotSEDs(freq0, flux0, base0, freq1=None, flux1=None, base1=None,StokesQU=None):
    """
    Plot one or two SEDs on the same axes.
    input:
    - freq: frequency array/list in units of GHz
    - flux: flux list/array in units of Jy/beam
    - base: string name/reference for input data
    - feq1, flux1 and base1 are optional
    - StokesQU: bolean - True/False, True if input flux is Stokes Q and U
      in which case linear polarization and polarization angle will be plotted
      with the assumption that flux0=Q and flux1=U.
    
    Returns:
    - PNG plots
    """
    import matplotlib.pyplot as plt
    c = 2.998e+08
    
    freq0 = array(freq0)
    flux0 = array(flux0)
    
    plt.figure()
    plt.title(str(base0)+', '+str(base1))
    
    plt.scatter(freq0, flux0,marker='+', label="%s"%str(base0))
    
    try:
        freq1 = array(freq1)
        flux1 = array(flux1)
        plt.scatter(freq1, flux1, label="%s"%str(base1))
    except Exception as ex00:
        print(' ! > ', ex)
        pass
    
    plt.grid(ls=':', alpha=0.5)
    plt.xlabel('Frequency / GHz')
    plt.ylabel('Flux Density / Jy/bm')
    plt.tight_layout()
    plt.legend(loc='best')
    plt.savefig('SED-%s-vs-%s.png'%(base0,base1))
    #plt.show()
    #plt.close()
    
    plt.figure()
    plt.title(str(base0)+' vs. '+str(base1))
    
    plt.loglog(freq0, flux0,'+', label="%s"%str(base0))
    
    try:
        freq1 = array(freq1)
        flux1 = array(flux1)
        plt.loglog(freq1, flux1, '.', label="%s"%str(base1))
    except Exception as ex00:
        print(' ! > ', ex)
        pass
    
    plt.grid(ls=':', alpha=0.5)
    plt.xlabel('Frequency / GHz')
    plt.ylabel('Flux Density / Jy/bm')
    plt.tight_layout()
    plt.legend(loc='best')
    plt.savefig('logSED-%s-vs-%s.png'%(base0,base1))
    #plt.show()
    #plt.close()
    
    if StokesQU:
        Q  = flux0
        U  = flux1
        try:
            #P  = sqrt(Q*Q + U*U)
            #PA = 0.5*arctan2(U,Q)*(180.0/pi)
            P, PA = a.get_LinPol_PA2(Q, U)
            
            plt.figure()
            plt.title('Linear pol: '+str(base0)+', '+str(base1))
            plt.scatter(freq0,P, marker='.', label='lin.pol')
            plt.grid(ls=':', alpha=0.5)
            plt.xlabel('Frequency / GHz')
            plt.ylabel('Flux Density / Jy/bm')
            plt.tight_layout()
            plt.legend(loc='best')
            plt.savefig('linPol-SED-%s-vs-%s.png'%(base0,base1))
            
            plt.figure()
            plt.title('Pol. ang: '+str(base0)+', '+str(base1))
            plt.scatter(freq0,PA, label='PA', marker='+')
            plt.grid(ls=':', alpha=0.5)
            plt.xlabel('Frequency / GHz')
            plt.ylabel('Pol. ang / deg')
            plt.tight_layout()
            plt.legend(loc='best')
            plt.savefig('PA-SED-%s-vs-%s.png'%(base0,base1))
        except Exception as ex002:
            print(ex002)
            pass
    plt.show()
    return()

def SED_cat1cat2(cat1, cat2,base1=None,base2=None,StokesQU=None):
    """
    Compare SEDs of two catalogues.
    input:
    - cat1: str() - name of 1st catalogue (.csv file)
    - cat2: str() - name of 2nd catalogue (.csv file)
    - catalogue must have columns named "#frq" and " fmx"
    - "#frq" column from cat2 is converted to GHz
    - base1/base2: string name/reference for input data from cat1/cat2
    - uses my function:  plotSEDs(feq0, flux0, base0, feq1=None, flux1=None, base1=None)
    
    Returns:
    - PNG figure from function plotSEDs(feq0, flux0, base0, feq1=None, flux1=None, base1=None)
    """
    import pandas as pd
    import os, sys
    
    df1  = pd.read_csv(cat1); df2  = pd.read_csv(cat2)
    try:
        freq_col = '#frq'
        fmax_col = ' fmx'
        frms_col = ' frms'
        
        frq1 = list(df1[freq_col]); frq2 = array(list(df2[freq_col])) #/1.0e9
        fmx1 = list(df1[fmax_col]); fmx2 = list(df2[fmax_col])
        fe1 = list(df1[frms_col]); fe2 = list(df2[frms_col])
        
    except KeyError as ke00:
        print(' >> ', ke00)
        freq_col = str(input(' Enter frequency columnn name: '))
        fmax_col = str(input(' Enter flux max columnn name : '))
        frms_col = str(input(' Enter flux RMS columnn name : '))
    
    frq1 = list(df1[freq_col]); frq2 = array(list(df2[freq_col])) #/1.0e9
    fmx1 = list(df1[fmax_col]); fmx2 = list(df2[fmax_col])
    fe1 = list(df1[frms_col]); fe2 = list(df2[frms_col])
    
    if base1 == None:
        base1 = os.path.basename(cat1)[:-4]
    if base2 == None:
        base2 = os.path.basename(cat2)[:-4]
    
    plotSEDs(frq1, fmx1, base1, frq2, fmx2, base2,StokesQU)
    plotSEDs(frq1, fe1, base1+'-RMS', frq2, fe2, base2+'-RMS',StokesQU)
    
    return()

def MKCosBeam(rho, nu):
    """
    Calculate cosine beam shape (Condon & Ransom, Essential Radio Astronomy eq 3.95)
   
    Return power gain of circularly symmetric beam
    * rho   = offset from center (degrees)
    * nu    = Frequency (Hz)
    """
    ################################################################
    from math import radians, cos, pi
    
    #theta_b = radians(57.5/60) * (1.5e9/nu)
    theta_b = 0.0167261 * (1.5e9/nu)
    rhor = 1.18896*radians(rho)/theta_b
    gain = (cos(pi*rhor)/(1.-4.*(rhor**2)))**2
    return gain
# end MKCosBeam

def spindx(x,a,b,c):
    """
    Curved power law function for Stokes I.
    
    Inputs:
    - x : SED frequenc array.
    - a : Amplitude coefficient of the power law 
    - b : 'intercept' coefficient in the log-linear exponent of the power law
    - c : 'slope' coefficient in the log-linear exponent of the power law
    
    Output:
    - SED power law fit as an array.
    """
    
    import numpy as np
    xref = np.nanmin(x) #Reference Frequency
    arg  = (x/xref)
    return(a * arg**(b+ c * np.log(arg)))

def logSpec(x, a , b):
    """
    Produce log-linear function of array x.
    Input:
    - x : x-axis array (not log)
    - a : log-linear slope
    - b : log-linear intercept
    """
    logx = log(x)
    logS = a*logx + b
    return (S)

def linefunc(x,a,b):
    return(a*x+b)


def ensamble_catalogue(abs_filedir,outname,prefix=None,suffix=None, xcol=None, ycol=None,zcol=None):
    """
    Produce catalog of an ensamble of source SEDs.
    
    input:
    - abs_filedir : path to SED files (atpy readable files - viz fits/ipac/vot tables)
    - outname     : human readable reference name.
    - prefix      : prefix file pattern
    - suffix      : suffix file pattern
    - x, y, z     : column names for x, y, z columns
    
    Outputs:
    - csv catalog
    - RA and DEC columns may be wrong as they are based on input file name. CHECK them.
    """
    
    import glob, atpy
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.optimize import curve_fit
    
    if prefix == None:
        flist = glob.glob(str(abs_filedir)+f"/*{suffix}")
    elif suffix == None:
        flist = glob.glob(str(abs_filedir)+f"/{prefix}*")
    else:
        flist = glob.glob(str(abs_filedir)+f"/{prefix}*{suffix}")
    
    if xcol == None:
        xcol = 'Freq'
    if ycol == None: 
        ycol = 'Ipeak'
    if zcol == None:
        zcol = 'Vpeak'
    
    cat_out = str(outname)+"_sources_summary.csv"
    catalog = open(str(abs_filedir)+f"/{cat_out}",'w')
    catalog.write('#RA_deg,DECdeg,Imed,Istd,Alpha,Vmed,Vstd,a,b,c\n')
    for f0 in flist:
        t    = atpy.Table(f0)
        base = os.path.basename(f0).rstrip('_IQU_perchan.tbl')
        try:
            RA  = float(base.split('-',1)[0])
            DEC = -1.*float(base.split('-',1)[1])
        except ValueError as voe000:
            print(f' >> ValueError {voe000}')
            RA  = '--'
            DEC = '--'
        
        SI0  = t[ycol]
        freq = np.array(t[xcol])[np.where(SI0*2 > SI0)]
        SI   = np.array(t[ycol])[np.where(SI0*2 > SI0)]
        SV   = np.array(t[zcol])[np.where(SI0*2 > SI0)]
        
        try:
            logf = np.log(freq)
            logS = np.log(SI)
        except TypeError as ex002:
            t.describe()
            print(f' >> TypeError: {ex002}')
            print(f' >> freqs: {freq}')
            print(f' >> SI   : {SI}\n\n\n')
        
            logf = [], logS = []
            for i in range(len(freq)):
                logf.append(np.log(freq[i]))
                logS.append(np.log(SI[i]))
            logf = np.array(logf)
            logS = np.array(logS)
            pass
        
        popt,pcov   = curve_fit(linefunc,logf,logS)
        try:
            popt2,pcov2 = curve_fit(spindx,freq,SI)
        except RuntimeError as re00:
            print(f' >> RuntimeError: {re00}')
            popt2 = ['--','--','--']
        outstr = f'{RA},{DEC},{np.nanmedian(SI)},{np.nanstd(SI)},{popt[0]},{np.nanmedian(SV)},{np.nanstd(SV)},{popt2[0]},{popt2[1]},{popt2[2]}'
        catalog.write(f'{outstr}\n')
        print(outstr)
        
    
    catalog.close()
    return()


def box_imstat(fitsimage, RAc, Decc, boxsize=None):
    from astropy.wcs import WCS
    """
    Runs inpure python
    Returns an imstat() 'box' parameter centered at 'RAc, Decc' and 
    with a square box region size = 'boxsize'.
    RAc, Decc, boxsize are in units of degress.
    Returned object is a set of intergers in the form: x0,y0,x1,y1 where:
     - x0,y0 is the bottom-left pixel coord
     - x1,y1 is the top-right pixel coord
    """
    dd, hdr = fits.getdata(fitsimage, header=True)
    fwcs    = WCS(hdr)
    
    if boxsize == None:
        boxsize = 30./3600. # 30 arcsec converted to degrees. MeerKAt beam ~ 10 arcsec
    
    print ("boxsize: %f deg "%boxsize)
    
    r0 = boxsize/2.
    
    ra0  = RAc - r0
    dec0 = Decc - r0
    ra1  = RAc + r0
    dec1 = Decc + r0
    
    pix0 = fwcs.wcs_world2pix(ra0,dec0,0,0,0)
    pix1 = fwcs.wcs_world2pix(ra1,dec1,0,0,0)
    
    x0 = int(round(pix0[0]))
    y0 = int(round(pix0[1]))
    x1 = int(round(pix1[0]))
    y1 = int(round(pix1[1]))
    
    return (x0,y0,x1,y1)
    

########################################################################################
if __name__ == '__main__':
    compare_NED = False
    do_QU       = True
    
    if do_QU:
        pdir   = "/home/samuel/UHF/caltables_final/plots/" 
        gcalQ  = pdir+'J0252-7104_1569087057.Q.plane--image.fits-imstats.csv'; gcalU = pdir+'J0252-7104_1569087057.U.plane--image.fits-imstats.csv'
        bcalQ  =  pdir+'J1939-6342_1569087057.Q.plane--image.fits-imstats.csv'; bcalU = pdir+'J1939-6342_1569087057.U.plane--image.fits-imstats.csv'
        B0408Q = pdir+'J0408-6545_1569087057.Q.plane--image.fits-imstats.csv'; B0408U = pdir+'J0408-6545_1569087057.U.plane--image.fits-imstats.csv'
        pcalQ  = pdir+'J0521_1638_1569087057.Q.plane--image.fits-imstats.csv'; pcalU = pdir+'J0521_1638_1569087057.U.plane--image.fits-imstats.csv' 
        
        SED_cat1cat2(B0408Q, B0408U,base1='J0408-6545-Q',base2='J0408-6545-U',StokesQU=do_QU)
        SED_cat1cat2(gcalQ, gcalU,base1='J0252-7104-Q',base2='J0252-7104-U',StokesQU=do_QU)
        SED_cat1cat2(bcalQ, bcalU,base1='J1939-6342-Q',base2='J1939-6342-U',StokesQU=do_QU)
        SED_cat1cat2(pcalQ, pcalU,base1='3C138-Q',base2='3C138-U',StokesQU=do_QU)
    
    if compare_NED:
        SEDdir1 = '/home/samuel/SARAO/UHF/final_images/SEDs/'
        NEDSEDs = '/home/samuel/SARAO/UHF/final_images/NED-data/'
        for c1 in os.listdir(SEDdir1): 
            for c2 in os.listdir(NEDSEDs): 
                if c1.startswith(str(c2)[0:10]) and c1.endswith(str(c2)[-4:]): 
                    cat1 = SEDdir1+c1 
                    cat2 = NEDSEDs+c2 
                    print (f' cat1: {cat1}\n cat2: {cat2}\n...\n') 
                    try: 
                        SED_cat1cat2(cat1, cat2,base1=None,base2=None) 
                    except Exception as ex: 
                        print('>>>',ex) 
                        pass 
    
    """
    corrdir = "/home/samuel/SARAO/imagin_varification/2020-correlator-test/April/"
    catdir  = corrdir+"images/catalogues/"
    L4k     = corrdir+"images/1586074843_continuum_image_J2147-8132_IClean.fits" 
    L1k     = corrdir+"images/1585989326_continuum_image_J2147-8132_IClean.fits"
    U4k     = corrdir+"images/1586070058_continuum_image_J2147-8132_IClean.fits"
    U1k     = corrdir+"images/1586162571_continuum_image_J2147-8132_IClean.fits"
    L4k_v2  = corrdir+"images/1586305859_continuum_image_J2147-8132_IClean.fits"
    matches = catdir+"L1KxU4KxL4KxU1K.tbl"; summsX1 = catdir+"L1K-x-U4K-x-L4Kv2-x-U1K-x-SUMSS.tbl"
    summsX2 = catdir+"L1K-x-U4K-x-L4Kv2-x-U1K-x-SUMSS.tbl"
    tx      = atpy.Table(matches, type='ipac'); smx = atpy.Table(summsX2, type='ipac')

    wL4k    = L4k[:-5]+".wideband.fits"
    wL1k    = L1k[:-5]+".wideband.fits"
    wU4k    = U4k[:-5]+".wideband.fits"
    wU1k    = U1k[:-5]+".wideband.fits"
    wL4k_v2 = L4k_v2[:-5]+".wideband.fits"
    wims    = {L4k:wL4k, L4k_v2:wL4k_v2, L1k:wL1k, U4k:wU4k, U1k:wU1k}
    D_imsts = {}

    #############################
    do_srcfind = False
    do_xmatch  = False
    do_imstats = False
    ############################# 

    if do_srcfind or do_xmatch:
        do_imstats  = False 
        try:
            import py3gokatsdpimager as g 
            import image_analysis as a
        except ImportError as ie2:
            print (" ! >" , ie2)
            print (" ! > run this part of code in python only, NOT casa!")
            raise

    for im0 in wims:
        base = os.path.basename(im0)[0:10]
        fdat, fhdr = fits.getdata(im0, header=True)
        print( "{} imsize: {} pix, CTYPE1: {} CDELT1: {:1.4f}, CRPIX1: {:1.4f}, CRVAL1: {:1.4f}, CRVAL2: {:1.4f}, CLEANBMJ: {:1.4f} deg, CLEANBMN: {:1.4f} deg, CLEANBPA: {:1.4f} deg, CLEANNIT: {:d} iterationsprint ...".format(base,fhdr['NAXIS1'], fhdr['CTYPE1'],fhdr['CDELT1'],fhdr['CRPIX1'], fhdr['CRVAL1'],fhdr['CRVAL2'],fhdr['CLEANBMJ'], fhdr['CLEANBMN'],fhdr['CLEANBPA'],fhdr['CLEANNIT'],) )

    if do_srcfind or do_imstats:
        for im0 in wims:
            if not os.path.isfile(wims[im0]):
                g.generic_get_subbandImage(im0)
            if do_srcfind:
                g.find_sources(wims[im0], None)
            if do_imstats:
                algo = "fit-half"
                s_i  = imstat(imagename=wims[im0],algorithm=algo,center='zero')
                D_imsts['%s rms,std [Jy/beam]'%base] = [s_i['rms'][0], s_i['sigma'][0]]
                print ("%s rms: %f Jy/beam"%(base, s_i['rms']))
        if do_imstats:
            print (" Image RMSs: \n {}".format(D_imsts) )

    ra_L1k, dec_L1k = tx['RA_1'], tx['DEC_1']
    ra_U4k, dec_U4k = tx['RA_2'], tx['DEC_2']
    ra_L4k, dec_L4k = tx['RA_3'], tx['DEC_3']
    ra_U1k, dec_U1k = tx['RA_4'], tx['DEC_4']

    e_sra, e_sdec     = smx['e_RAJ2000'], smx['e_DEJ2000']
    e_L1kra, e_L1kdec = smx['E_RA_1']*3600., smx['E_DEC_1']*3600.
    e_L4kra, e_L4kdec = smx['E_RA_3']*3600., smx['E_DEC_3']*3600.
    e_U1kra, e_U1kdec = smx['E_RA_4']*3600., smx['E_DEC_4']*3600.
    e_U4kra, e_U4kdec = smx['E_RA_2']*3600., smx['E_DEC_2']*3600.

    avg_esra, std_esra   = getAvgStdev(e_sra)
    avg_esdec, std_esdec = getAvgStdev(e_sdec)
    avg_eL1ra, std_eL1ra   = getAvgStdev(e_L1kra)
    avg_eL1dec, std_eL1dec = getAvgStdev(e_L1kdec)
    avg_eL4ra, std_eL4ra   = getAvgStdev(e_L4kra)
    avg_eL4dec, std_eL4dec = getAvgStdev(e_L4kdec)
    avg_eU4ra, std_eU4ra   = getAvgStdev(e_U4kra)
    avg_eU4dec, std_eU4dec = getAvgStdev(e_U4kdec)
    avg_eU1ra, std_eU1ra   = getAvgStdev(e_U1kra)
    avg_eU1dec, std_eU1dec = getAvgStdev(e_U1kdec)

    print ( "\n\nSUMSS eRA avg : avg: {:1.4}, std: {:1.4} asec".format(avg_esra, std_esra) ) 
    print ( "SUMSS eDec avg: avg: {:1.4}, std: {:1.4} asec".format(avg_esdec, std_esdec) ) 
    print ( "L1k eRA avg   : avg: {:1.4}, std: {:1.4} asec".format(avg_eL1ra, std_eL1ra) ) 
    print ( "L1k eDec avg  : avg: {:1.4}, std: {:1.4} asec".format(avg_eL1dec, std_eL1dec) ) 
    print ( "L4k eRA avg   : avg: {:1.4}, std: {:1.4} asec".format(avg_eL4ra, std_eL4ra) ) 
    print ( "L4k eDec avg  : avg: {:1.4}, std: {:1.4} asec".format(avg_eL4dec, std_eL4dec) ) 
    print ( "U4k eRA avg   : avg: {:1.4}, std: {:1.4} asec".format(avg_eU4ra, std_eU4ra) ) 
    print ( "U4k eDec avg  : avg: {:1.4}, std: {:1.4} asec".format(avg_eU4dec, std_eU4dec) ) 
    print ( "U4k eRA avg   : avg: {:1.4}, std: {:1.4} asec".format(avg_eU1ra, std_eU1ra) ) 
    print ( "U4k eDec avg  : avg: {:1.4}, std: {:1.4} asec".format(avg_eU1dec, std_eU1dec) ) 

    if do_xmatch:
        L4x1k_idx, L4x1k_sep2d, L4x1k_sep3d, L4x1k_rax, L4x1k_decx, draL, ddecL = a.CrossMatch(ra_L4k, dec_L4k, ra_L1k, dec_L1k)
        U4x1k_idx, U4x1k_sep2d, U4x1k_sep3d, U4x1k_rax, U4x1k_decx, draU, ddecU = a.CrossMatch(ra_U4k, dec_U4k, ra_U1k, dec_U1k)
        LxU1k_idx, LxU1k_sep2d, LxU1k_sep3d, LxU1k_rax, LxU1k_decx, draLU1, ddecLU1 = a.CrossMatch(ra_L1k, dec_L1k, ra_U1k, dec_U1k)
        LxU4k_idx, LxU4k_sep2d, LxU4k_sep3d, LxU4k_rax, LxU4k_decx, draLU4, ddecLU4 = a.CrossMatch(ra_L4k, dec_L4k, ra_U4k, dec_U4k)
        
        avg_draL, std_draL       = nanmean(draL)*3600., nanstd(draL)*3600.
        avg_ddecL, std_ddecL     = nanmean(ddecL)*3600., nanstd(ddecL)*3600.
        avg_draU, std_draU       = nanmean(draU)*3600., nanstd(draU)*3600.
        avg_ddecU, std_ddecU       = nanmean(ddecU)*3600., nanstd(ddecU)*3600.
        avg_draLU1, std_draLU1   = nanmean(draLU1)*3600., nanstd(draLU1)*3600.
        avg_ddecLU1, std_ddecLU1 = nanmean(ddecLU1)*3600., nanstd(ddecLU1)*3600.
        avg_draLU4, std_draLU4   = nanmean(draLU4)*3600., nanstd(draLU4)*3600.
        avg_ddecLU4, std_ddecLU4 = nanmean(ddecLU4)*3600., nanstd(ddecLU4)*3600.
        
        print ( "\n\nL 4kx1k avg_draL, std_draL     : {:1.4f}, {:1.4f} asec".format(avg_draL, std_draL) )
        print ( "L 4kx1k avg_ddecL, std_ddecL   : {:1.4f}, {:1.4f} asec".format(avg_ddecL, std_ddecL) )
        print ( "UHF 4kx1k avg_draL, std_draL   : {:1.4f}, {:1.4f} asec".format(avg_draU, std_draU) )
        print ( "UHF 4kx1k avg_ddecU, std_ddecU : {:1.4f}, {:1.4f} asec".format(avg_ddecU, std_ddecU) )
        print ( "LxUHF 1k avg_dra, std_dra      : {:1.4f}, {:1.4f} asec".format(avg_draLU1, std_draLU1) )
        print ( "LxUHF 1k avg_ddec, std_ddec    : {:1.4f}, {:1.4f} asec".format(avg_ddecLU1, std_ddecLU1) )
        print ( "LxUHF 4k avg_dra, std_dra      : {:1.4f}, {:1.4f} asec".format(avg_draLU4, std_draLU4) )
        print ( "LxUHF 4k avg_ddec, std_ddec    : {:1.4f}, {:1.4f} asec\n\n".format(avg_ddecLU4, std_ddecLU4) )

        plt.figure(figsize=(6,15))
        plt.subplot(311)
        plt.hist(e_sra, bins=get_bins(e_sra), label='SUMSS eRA'); plt.grid(ls=':', alpha=0.5)
        plt.hist(e_L4kra, bins=get_bins(e_L4kra),linestyle=('--'), histtype='step',label='L4k eRA'); plt.grid(ls=':', alpha=0.5)
        plt.hist(e_L1kra, bins=get_bins(e_L1kra), histtype='step',label='L1k eRA'); plt.grid(ls=':', alpha=0.5)
        plt.hist(e_U1kra, bins=get_bins(e_U1kra),linestyle=('-.'), histtype='step',label='U1k eRA'); plt.grid(ls=':', alpha=0.5)
        plt.hist(e_U4kra, bins=get_bins(e_U4kra), histtype='step',label='U4k eRA'); plt.grid(ls=':', alpha=0.5)
        plt.legend(loc=0)
        
        plt.subplot(312)
        plt.hist(e_sdec, bins=get_bins(e_sdec),label='SUMSS eDEC'); plt.grid(ls=':', alpha=0.5)
        plt.hist(e_L1kdec, bins=get_bins(e_L1kdec),linestyle=('-.'), histtype='step',label='L1k eDec'); plt.grid(ls=':', alpha=0.5)
        plt.hist(e_L4kdec, bins=get_bins(e_L4kdec), histtype='step',label='L4k eDec'); plt.grid(ls=':', alpha=0.5)
        plt.legend(loc=0)
        
        plt.subplot(313)
        plt.hist(e_sdec, bins=get_bins(e_sdec),label='SUMSS eDEC'); plt.grid(ls=':', alpha=0.5)
        plt.hist(e_U1kdec, bins=get_bins(e_U1kdec),linestyle=('-.'), histtype='step',label='U1k eDec'); plt.grid(ls=':', alpha=0.5)
        plt.hist(e_U4kdec, bins=get_bins(e_U4kdec), histtype='step',label='U4k eDec'); plt.grid(ls=':', alpha=0.5)
        plt.xlabel('RA / Dec errors [arcsec]')
        plt.legend(loc=0)
        plt.savefig(catdir+"/e_posHists.png")

        plt.figure()
        plt.hist(L4x1k_sep2d.arcsec, bins=34, histtype='step', label='L: 4kx1k')
        plt.hist(U4x1k_sep2d.arcsec, bins=34, histtype='step', label='UHF: 4kx1k')
        plt.hist(LxU1k_sep2d.arcsec, bins=34, histtype='step', label='LxU: 1k')
        plt.hist(LxU4k_sep2d.arcsec, bins=34, histtype='step', label='LxU: 4k')
        plt.grid(ls=':', alpha=0.5)
        
        plt.xlabel('separation [arcsec]')
        plt.legend(loc=0)
        plt.savefig(catdir+"/Xmatch-sep2d.png")
        #plt.show(),plt.close()

        plt.figure()
        plt.title('Offsets in RA and Dec at \nL and UHF bands for cross \nmatched channelisations')
        plt.plot(draL*3600., ddecL*3600., 'r+', label='L: 4kx1k')
        plt.plot(draU*3600., ddecU*3600., 'b.', label='UHF: 4kx1k')
        plt.grid(ls=':', alpha=0.5)
        
        plt.xlabel('dRA [arcsec]')
        plt.ylabel('dDec [arcsec]')
        plt.legend(loc=0)
        plt.savefig(catdir+"/4kx1k-dRAvsDdec.png")
        #plt.show()

        plt.figure()
        plt.title('Offsets in RA and Dec for \ncross matched L and UHF bands \nat the same channelisations')
        plt.plot(draLU1*3600., ddecLU1*3600., 'm+', label='LxU: 1k')
        plt.plot(draLU4*3600., ddecLU4*3600., 'c.', label='LxU: 4k')
        plt.grid(ls=':', alpha=0.5)
        
        plt.xlabel('dRA [arcsec]')
        plt.ylabel('dDec [arcsec]')
        plt.legend(loc=0)
        plt.savefig(catdir+"/LxU-4k+1k-dRAvsDdec.png")
        plt.show()
        """
