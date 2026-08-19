"""
This code is import friendly.
"""

import katdal, os, sys, time, atpy, argparse, csv
from astropy.table import Table,Column
from astropy.io import fits
from astropy.wcs import WCS as w
from matplotlib import pyplot as plt
from matplotlib.pyplot import *
from astropy.wcs import WCS
import pandas as pd
from math import*
from numpy import*
from srcfind import*
import meerkatresults as myfns
import itertools
from astropy import units as u
from astropy.coordinates import SkyCoord
try:
    import py3gokatsdpimager as g 
except ImportError as ie2:
    print (" ! >" , ie2, "\n\n\n")
    print (" ! > run this part of code in python only, NOT casa!")
    print (" ! OR make sure that all imported modules exist in the working directory!")
    pass
####################################################
plt.close()
def safe_div(x,y):
    try:
        if y == 0.:
            return 0
    
    except ValueError:
        for yi in y:
            if yi == 0.:
                return 0
    return (x/y)

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

#def separation(ra1,dec1,ra2,dec2):
    #"""This function returns the angular separation between two points defined by
    #ra1,dec1, and ra2,dec2.   Input and output is in degrees """
    #pi = 3.1415928
    #deg2rad = pi/180.0
    #sep = pyasl.getAngDist(float(ra1), float(dec1), float(ra2), float(dec2))
    #return(sep)
    
def separation(RA1, Dec1, RA2, Dec2):
    """
    # Simple approximation of angular separation between two catalogues ( which may or not be cross-matched)
    # - assumes coordinates from the same reference frame
    # - pyasl.getAngDist() seems to be wrong in some instances!
    # - assumes inputs in deg and so converts Dec1 into rad in cos()
    """
    
    if hasattr(RA1, '__len__'):
        RA1  = array(RA1)
        Dec1 = array(Dec1)
        RA2  = array(RA2)
        Dec2 = array(Dec2)
    dRA  = (RA2-RA1)*cos(radians(Dec1))
    dDec = Dec2-Dec1
    
    return(dRA, dDec)


#def separation_asec(ra1,dec1,ra2,dec2):
    #"""This function returns the angular separation between two points defined by
    #ra1,dec1, and ra2,dec2.   Input is in degrees, and output arcsec """
    #pi = 3.1415928
    #deg2rad = pi/180.0
    #sep = pyasl.getAngDist(float(ra1), float(dec1), float(ra2), float(dec2))
    #return(sep*3600.)

def separation_asec(RA1, Dec1, RA2, Dec2):
    """
    # assumes coordinates from the same reference frame
    # simple approximation of angular separation
    # pyasl.getAngDist() seems to be wrong in some instances!
    # assumes inputs in deg and so converts Dec1 into rad in cos()
    """
    if hasattr(RA1, '__len__'):
        RA1  = array(RA1)
        Dec1 = array(Dec1)
        RA2  = array(RA2)
        Dec2 = array(Dec2)
    dRA  = (RA2-RA1)*cos(radians(Dec1))
    dDec = Dec2-Dec1
    
    return(dRA*3600., dDec*3600.)


def substrctr_fts(fits_catalogue, perc_cut=None):
    """
    Picks out secondary sources in a PYBDSF produced Gaussian sky model.
    - Inputs:
        fits_catalogue: FITS table of Gaussian source components from bdsf.process_image()
        perc_cut: cutoff of secondary source flux as a fraction of the peak flux in the field
    """
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

def pos_varCMC1xCMC2(ra1, dec1, x1, y1, ra2, dec2, x2, y2, base1, base2):
    """# ra, dec, x, y need to be arrays or lists
    # tested on the same image from cmc1 and cmc2 
    # basei defines the data in cmci"""
    
    n1 = list(range(len(ra1))); n2 = list(range(len(ra2)))
    x0, y0 = 5670.,5670.
    Dx, Dy, phi_d, Tht_asec, Tht_r = [],[],[],[],[]
    theta = []
    Dx2, Dy2, phi2_d, Tht_asec2, Tht_r2 = [],[],[],[],[]
    
    for j in n1:
        dx1, dy1 = x1[j]-x0, y1[j]-y0
        dx2, dy2 = x2[j]-x0, y2[j]-y0
        
        phi  = atan2(dy1,dx1)
        phi2 = atan2(dy2,dx2)
        phi_di = phi*raddeg
        phi2_di = phi2*raddeg 
        #print phi_di/cbmaj
        phi_asc = phi_di*3600.
        Dx.append(dx1), Dy.append(dy1), phi_d.append(phi_di), Tht_asec.append(phi_asc), Tht_r.append(phi)
        Dx2.append(dx2), Dy2.append(dy2), phi2_d.append(phi2_di)
    
    for i in n1:
        drai = (ra2[i] - ra1[i])**cos((dec1[i])*(pi/180.))
        ddeci= dec2[i] - dec1[i]
        thetai = atan2(ddeci,drai)*raddeg
        
        theta.append( thetai )
    
    Dra, Ddec = RA_Dec_offsets(ra1, dec1, ra2, dec2)
    Ddec = Ddec/cbmaj
    Dra  = Dra/cbmaj
    
    phi_d = array(phi_d)
    phi2_d = array(phi2_d)
    theta = array(theta)/cbmaj
    
    
    try:
        plt.figure()
        plt.subplot(211)
        plt.title('%s x %s' %(base1, base2))
        plt.plot(Dra,phi_d,'r+', label='$\phi_1$')
        plt.plot(Dra,phi2_d,'k.', label='$\phi_2$')
        plt.xlabel('$\Delta$RA [clean bmaj]')
        plt.ylabel('$\phi = \\arctan(\Delta y/\Delta x)$\n [deg]')
        plt.legend(loc=0, fancybox=True, framealpha=0.35)
        plt.grid(ls=':', alpha=0.5)
        plt.subplot(212)
        plt.plot(Ddec,phi_d,'r+', label='$\phi_1$')
        plt.plot(Ddec,phi2_d,'k.', label='$\phi_2$')
        plt.xlabel('$\Delta$Dec [clean bmaj]')
        plt.ylabel('$\phi$ [deg]')
        plt.legend(loc=0, fancybox=True, framealpha=0.35)
        plt.grid(ls=':', alpha=0.5)
        plt.tight_layout()
        plt.savefig('DRA-vs-phi-%sx%s.png'%(base1, base2))
        plt.show()
        plt.close()
    except Exception as ex001:
        print(ex001)
        plt.close()
        pass
    
    
    plt.title('%s x %s' %(base1, base2))
    plt.plot(Dra,Ddec,'.')
    plt.xlabel('$\Delta$RA / clean bmaj')
    plt.ylabel('$\Delta$Dec / clean bmaj')
    plt.grid(ls=':', alpha=0.5)
    plt.tight_layout()
    plt.savefig('DRA-vs-Ddec-%sx%s.png'%(base1, base2))
    plt.show()
    
    plt.title('%s x %s' %(base1, base2))
    plt.plot(theta, phi_d,'r+', label='$\phi_1$')
    plt.plot(theta, phi2_d,'k.', label='$\phi_2$')
    plt.xlabel('$\\theta = \\arctan(\Delta Dec/\Delta RA)$\n [clean bmaj]')
    plt.legend(loc=0, fancybox=True, framealpha=0.35)
    plt.ylabel('$\phi$ [deg]')
    plt.grid(ls=':', alpha=0.5)
    plt.tight_layout()
    plt.savefig('theta-vs-phi-%sx%s.png'%(base1, base2))
    plt.show()
    
    return ()

def get_SUMSSXcmcX(tablename):
    """
    Astrometry of a sources that have been X-matched between SUMSS data and another catalog.
    """
    # subscript SM = SUMSS data
    # subscript cmcX = cmcX data, X = 1 or 2, etc.
    
    t        = atpy.Table(tablename)
    ra_cmcX  = t['RA']; era_cmcX = t['E_RA']
    ra_SM    = t['_RAJ2000']; era_SM = t['e_RAJ2000']
    dec_cmcX = t['DEC']; edec_cmcX = t['E_DEC']
    dec_SM   = t['_DEJ2000']; edec_SM = t['e_DEJ2000']
    xpx_cmcX = t['Xposn']; ypx_cmcX = t['Yposn']
    xpx_SM   = t['Xpos']; ypx_SM = t['Ypos']
    
    Dra,Ddec = RA_Dec_offsets(ra_cmcX, dec_cmcX, ra_SM, dec_SM)
    
    print("> Position errors:")
    print((" - mean(Delta_RA) : {}, sd(Delta_RA) : {} [cbmaj]".format( nanmean(Dra)/cbmaj, nanstd(Dra)/cbmaj )))
    print((" - mean(Delta_Dec): {}, sd(Delta_Dec): {} [cbmaj]".format( nanmean(Ddec)/cbmaj, nanstd(Ddec)/cbmaj )))
    
    return (ra_cmcX, era_cmcX, ra_SM, era_SM, dec_cmcX, edec_cmcX, dec_SM, edec_SM, Dra, Ddec, xpx_cmcX, ypx_cmcX, xpx_SM, ypx_SM)    

def RA_Dec_offsets(RA1, Dec1, RA2, Dec2):
    """
    # assumes coordinates from the same reference frame
    # simple approximation of angular separation
    # pyasl.getAngDist() seems to be wrong in some instances!
    # assumes inputs in deg and so converts Dec1 into rad in cos()
    output:
        - delta(RA) and delta(Dec) in deg
    """
    
    if hasattr(RA1, '__len__'):
        RA1  = array(RA1)
        Dec1 = array(Dec1)
        RA2  = array(RA2)
        Dec2 = array(Dec2)
    dRA  = (RA2-RA1)*cos(radians(Dec1))
    dDec = Dec2-Dec1
    
    return(dRA, dDec)

def get_cmc2Xcmc1(tablename):
    t        = atpy.Table(tablename)
    ra_cmc1  = t['RA_1']; era_cmc1 = t['E_RA_1']
    ra_cmc2  = t['RA_2']; era_cmc2 = t['E_RA_2']
    dec_cmc1 = t['DEC_1']; edec_cmc1 = t['E_DEC_1']
    dec_cmc2 = t['DEC_2']; edec_cmc2 = t['E_DEC_2']
    xpx_cmc1 = t['Xposn_1']; ypx_cmc1 = t['Yposn_1']
    xpx_cmc2 = t['Xposn_2']; ypx_cmc2 = t['Yposn_2']
    theta    = []
    
    for i in range(len(ra_cmc1)):
        drai  = ra_cmc1[i] - ra_cmc2[i]
        ddeci = dec_cmc1[i] - dec_cmc2[i]
        theta.append( atan2(ddeci,drai)*raddeg )
        
    theta= array(theta)
    
    Dra, Ddec = RA_Dec_offsets(ra_cmc1, dec_cmc1, ra_cmc2, dec_cmc2)
    
    return (ra_cmc1, era_cmc1, ra_cmc2, era_cmc2, dec_cmc1, edec_cmc1, dec_cmc2, edec_cmc2, Dra, Ddec, xpx_cmc1, ypx_cmc1, xpx_cmc2, ypx_cmc2,theta)

def CrossMatch(RA1_arr, DEC1_arr, RA2_arr, DEC2_arr):
    """
    Cross matches two sets of sky coordinates (RA and Dec) in units of degrees.
    Finds the nearest on-sky matches of this coordinate in a set of catalog coordinates using 
    SkyCoord.match_to_catalog_sky()
    
        - Input : 
                2 sets of float or array-like RA and Dec in degrees that are converted to 
                astropy coordinates objects, crd1 and crd2, via SkyCoord().
                
        - Output: 
            idx    - (integer array) the indecies of crd2 that get the closest matches to elements in crd1 
            sep2d  - (astropy.coordinates.angle()) the on-sky distances between the matches.
            dist3d - the real 3D space distances between matches 
            (ignore dist3d since we don't specify LoS component of coordinates)
            RA2_arr[idx], DEC2_arr[idx] are the closest matched RA and Dec from crd2 to each RA and Dec of crd1.
            dra, ddec - delta(RA) and delta(Dec) in deg
    """
    #print(CrossMatch.__doc__)
    
    crd1 = SkyCoord(RA1_arr*u.deg, DEC1_arr*u.deg)
    crd2 = SkyCoord(RA2_arr*u.deg, DEC2_arr*u.deg)
    idx, sep2d, dist3d = crd1.match_to_catalog_sky(crd2)
    dra, ddec = RA_Dec_offsets(RA1_arr, DEC1_arr, RA2_arr, DEC2_arr)
    
    return (idx, sep2d, dist3d, RA2_arr[idx], DEC2_arr[idx], dra, ddec)

def get_LinPol_PA(Q, U, I=None, Qe=None, Ue=None):
    """
    Return the polarizzation angle (in deg) and fractional linear polarization
    given input Stokes Q, U, I, Q error, and U error 
    """
    import numpy as np
    
    PA = 0.5*np.arctan2(U,Q)*(180.0/pi)
    p0 = sqrt(Q*Q + U*U)/I
    
    if Qe != None:
        pe = sqrt(Qe*Qe + Ue*Ue)/I
    elif Qe == None:
        pe = nan
    return(p0,PA,pe)

def get_LinPol_PA2(Q, U):
    P, PA = [], []
    
    for qu in itertools.zip_longest(Q, U):
        q = qu[0]
        u = qu[1]
        
        if q != None and u != None:
            p = sqrt(q**2 + u**2)
            pa = 0.5*np.arctan2(u,q)*(180.0/pi)
        elif q == None or u == None:
            p = np.nan
            pa = np.nan
            pass
        print(f'P: {p}, PA: {pa} +++ Q: {q}, U:{u}')
        P.append(p)
        PA.append(pa)
    P, PA = np.array(P), np.array(PA)
    print (f'avg P: {nanmean(P)}, rms P: {nanstd(P)}, avgPA: {nanmean(PA)}')
    return(P, PA)
    
####################################################
if __name__ == '__main__':
    Xdir     = '/home/samuel/SARAO/imagin_varification/catalogues/SUMSS-xmatches/'
    tcmc1_1k = Xdir+'SUMSSxCMC1-1K-1565564738-J2147-8132.tbl'      # X-matched with SUMSS
    tcmc1_4k = Xdir+'SUMSSxCMC1-4K-1565651658-J2147-8132.tbl'      # X-matched with SUMSS
    tcmc2_1k = Xdir+'SUMSSxCMC2-1K-1565560531-J2147-8132.tbl'      # X-matched with SUMSS
    tcmc2_4k = Xdir+'SUMSSxCMC2-4K-1565646959-J2147-8132.tbl'      # X-matched with SUMSS
    tcmc2_1kXcmc1_1k = Xdir+'CMC2-1K-x-CMC1-1K-1565560531x1565564738.tbl' 
    tcmc2_1kXcmc1_4k = Xdir+'CMC2-1K-x-CMC1-4K-1565560531x1565651658.tbl' 
    tcmc2_4kXcmc1_1k = Xdir+'CMC2-4K-x-CMC1-1K-1565646959x1565564738.tbl' 
    tcmc2_4kXcmc1_4k = Xdir+'CMC2-4K-x-CMC1-4K-1565646959x1565651658.tbl' 

    cbmaj  = 2.483482E-03        # in deg
    cbmin  = 2.133589E-03        # in deg
    raddeg = 180./pi
    show_DraDdec = False
    show_rot = False
    do_test = True

    ##### POSITIONS #################
    # Xmatched with SUMSS
    print('\nCMC1 1k x SUMSS')
    ra_1kcmc1, era_1kcmc1, ra_1kcmc1xSM, era_1kcmc1xSM, dec_1kcmc1, edec_1kcmc1, dec_1kcmc1xSM, edec_1kcmc1xSM,\
        Dra_1kcmc1xSM, Ddec_1kcmc1xSM, xpx_1kcmc1, ypx_1kcmc1, xpx_1kcmc1xSM, ypx_1kcmc1xSM = get_SUMSSXcmcX(tcmc1_1k)

    print('\nCMC1 4k x SUMSS')
    ra_4kcmc1, era_4kcmc1, ra_4kcmc1xSM, era_4kcmc1xSM, dec_4kcmc1, edec_4kcmc1, dec_4kSM, edec_4kSM,\
        Dra_4kcmc1xSM, Ddec_4kcmc1xSM, xpx_4kcmc1, ypx_4kcmc1, xpx_4kcmc1xSM, ypx_4kcmc1xSM = get_SUMSSXcmcX(tcmc1_4k)

    print('\nCMC2 1k x SUMSS')
    ra_1kcmc2, era_1kcmc2, ra_1kcmc2xSM, era_1kcmc2xSM, dec_1kcmc2, edec_1kcmc2, dec_1kcmc2xSM, edec_1kcmc2xSM,\
        Dra_1kcmc2xSM, Ddec_1kcmc2xSM, xpx_1kcmc2, ypx_1kcmc2, xpx_1kcmc2xSM, ypx_1kcmc2xSM = get_SUMSSXcmcX(tcmc2_1k)

    print('\nCMC2 4k x SUMSS')
    ra_4kcmc2, era_4kcmc2, ra_4kcmc2xSM, era_4kcmc2xSM, dec_4kcmc2, edec_4kcmc2, dec_4kcmc2xSM, edec_4kcmc2xSM,\
        Dra_4kcmc2xSM, Ddec_4kcmc2xSM, xpx_4kcmc2, ypx_4kcmc2, xpx_4kcmc2xSM, ypx_4kcmc2xSM = get_SUMSSXcmcX(tcmc2_4k)

    ####### PLOTS ##################
    if show_rot:
        t1xt2 = tcmc2_4kXcmc1_1k
        base1 = 'CMC2-4K'; base2 = 'CMC1-4K' # 
        ra_2, era_2, ra_1, era_1, dec_2, edec_2, dec_1, edec_1, Dra2x1, Ddec2x1, xpx_2, ypx_2,xpx_1, ypx_1,\
            theta2x1 = get_cmc2Xcmc1(t1xt2)
        pos_varCMC1xCMC2(ra_1, dec_1, xpx_1, ypx_1, ra_2, dec_2, xpx_2, ypx_2, base1,base2)
        


    if show_DraDdec:
        plt.title('SUMSS x '+os.path.basename(tcmc1_1k)[-14:-4])
        plt.plot(Dra_1kcmc1xSM/cbmaj, Ddec_1kcmc1xSM/cbmaj, 'r+', label=os.path.basename(tcmc1_1k)[6:-15])
        plt.xlabel('$\Delta$RA/clean bmaj')
        plt.ylabel('$\Delta$Dec/clean bmaj')
        plt.legend(loc=0, fancybox=True, framealpha=0.35)
        plt.grid(ls=':', alpha=0.5)
        plt.tight_layout()
        plt.savefig('DRAvsDDEC_cmc11KXsumss.png')
        #plt.show()
        plt.close()


        #plt.title(os.path.basename(tcmc1_4k)[:-4])
        plt.plot(Dra_4kcmc1xSM/cbmaj, Ddec_4kcmc1xSM/cbmaj, 'k.', label=os.path.basename(tcmc1_4k)[6:-15])
        plt.xlabel('$\Delta$RA/clean bmaj')
        plt.ylabel('$\Delta$Dec/clean bmaj')
        plt.legend(loc=0, fancybox=True, framealpha=0.35)
        plt.grid(ls=':', alpha=0.5)
        plt.tight_layout()
        plt.savefig('DRAvsDDEC_cmc14KXsumss.png')
        #plt.show()
        plt.close()


        #plt.title(os.path.basename(tcmc2_1k)[:-4])
        plt.plot(Dra_1kcmc2xSM/cbmaj, Ddec_1kcmc2xSM/cbmaj, 'r+', label=os.path.basename(tcmc2_1k)[6:-15])
        plt.xlabel('$\Delta$RA/clean bmaj')
        plt.ylabel('$\Delta$Dec/clean bmaj')
        plt.legend(loc=0, fancybox=True, framealpha=0.35)
        plt.grid(ls=':', alpha=0.5)
        plt.tight_layout()
        plt.savefig('DRAvsDDEC_cmc21KXsumss.png')
        #plt.show()
        plt.close()


        #plt.title(os.path.basename(tcmc2_4k)[:-4])
        plt.plot(Dra_4kcmc2xSM/cbmaj, Ddec_4kcmc2xSM/cbmaj, 'k.', label=os.path.basename(tcmc2_4k)[6:-15])
        plt.xlabel('$\Delta$RA/clean bmaj')
        plt.ylabel('$\Delta$Dec/clean bmaj')
        plt.legend(loc=0, fancybox=True, framealpha=0.35)
        plt.grid(ls=':', alpha=0.5)
        plt.tight_layout()
        plt.savefig('DRAvsDDEC_cmc24KXsumss.png')

        #plt.savefig('DRAvsDDEC_cmcxXsumss.png')
        plt.show()
        
    if do_test:
        plt.figure()
        plt.subplot(211)
        plt.plot(Dra_4kcmc1xSM/cbmaj, Ddec_4kcmc1xSM/cbmaj, 'k.', label=os.path.basename(tcmc1_4k)[6:-15])
        plt.plot(Dra_4kcmc2xSM/cbmaj, Ddec_4kcmc2xSM/cbmaj, 'k+', label=os.path.basename(tcmc2_4k)[6:-15])
        plt.xlabel('$\Delta$RA/clean bmaj')
        plt.ylabel('$\Delta$Dec/clean bmaj')
        plt.legend(loc=0, fancybox=True, framealpha=0.35)
        plt.grid(ls=':', alpha=0.5)
        plt.tight_layout()
        plt.subplot(212)
        plt.plot(Dra_1kcmc2xSM/cbmaj, Ddec_1kcmc2xSM/cbmaj, 'r+', label=os.path.basename(tcmc2_1k)[6:-15])
        plt.plot(Dra_1kcmc1xSM/cbmaj, Ddec_1kcmc1xSM/cbmaj, 'r.', label=os.path.basename(tcmc1_1k)[6:-15])
        plt.xlabel('$\Delta$RA/clean bmaj')
        plt.ylabel('$\Delta$Dec/clean bmaj')
        plt.legend(loc=0, fancybox=True, framealpha=0.35)
        plt.grid(ls=':', alpha=0.5)
        plt.tight_layout()
        plt.show()
        plt.close()

