#!/usr/bin/env python3
"""
- Script (python3) to perform analysis on images produced via the katsdp imager (katsdpcontim)
- import friendly
"""
from __future__ import print_function
import katdal
import os, sys, time, atpy, argparse, csv, subprocess
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
from PyAstronomy import pyasl
import random as rand
####################################################
print (__doc__)
####################################################
# command line arguments:
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-grab_subbandIms', action='store_true', help='True/False: Extract sub band images from input cuboid.')
parser.add_argument('-do_SourceFinding', action='store_true', help='True/False: Do source finding via bdsf.process_image().')
parser.add_argument("-src_findIm", required=False, help='Fits file for sourcefinding. Required if find_sources = True!')
parser.add_argument('-do_Xmatch', action='store_true', help='True/False: Do source cross-matching via astropy.coordinates.match_to_catalog_sky().')
parser.add_argument('-check_rotation', action='store_true', help='True/False: Do image rotation analyses.')
parser.add_argument("-SUMSS_cat", required=False, help='SUMSS catalog file path')
parser.add_argument("-MKAT_cat", required=False, help='MeerKAT catalog file path')
parser.add_argument("-in_fitscube", required=False, help='Cuboid fits file from katsdpimager from which to get plane images via the "grab_subbandIms" parameter') 
parser.add_argument("-xfile", required=False, help='ATPY readable file containing all cross-matches between SUMMS and MeerKAT catalogs') 
parser.add_argument("-cwd", required=False, default=os.getcwd(), help="Path. Not required. Path to working directory.")
parser.add_argument("-test_descript", required=False, help="Required! Human readable short one word description of dataset/correlator test")

parser.parse_args()
args = parser.parse_args()

####################################################
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

def find_sources(fits_in,base=None, pwd=None, imagedir=None, show_fitted=None):
    from astropy.io import fits 
    import os
    from os.path import dirname, abspath
    
    """
    Function takes in a fits image (not a fits cube) and searches for sources
    using the PYBDSF task process_image(). Thresholds are specified within 
    this function but may later be added into the arguments.
    """
    
    if pwd == None:
        pwd = os.getcwd()+'/'
        
    if imagedir == None:
        imagedir = os.path.dirname(abspath(fits_in))
        
    if show_fitted == None:
        show_fitted = False
        
    try:
        import bdsf
        # PYBDSF source finding in python
        os.chdir(imagedir)
        print( f'>> Working in directory: {imagedir}' )
        fitsfile = os.path.basename(fits_in)      # assumes you're in the 'fits_in' parent directory
        
        if base == None:
            base = os.path.basename(fitsfile)[0:10]
        
        n_cores        = 16
        snrthreshold   = 10.0                # SNR limit for good sources
        isl_thresh     = 3.0                # SNR limit for source finding island region
        pix_thresh     = 1.5*isl_thresh      # SNR limit for source detection/identification
        fluxthreshold  = 1.0e-6
        pblimit        = 0.05  #this should be constant at about 5%, hence 0.05 
        
        out0    = imagedir+'/pybdsf.results/'
        outdir  = out0+base+'/'
        
        if not (os.path.isdir(out0)):
            os.mkdir(out0)
            os.mkdir(outdir)
        if not (os.path.isdir(outdir)):
            os.mkdir(outdir)
        
        print('> PYBDSF OUTPUT dir: {}\n'.format(outdir))

        try:
            img = bdsf.process_image(fitsfile, output_opts=True,output_all=True,shapelet_do=True,solnname = base,
                                     thresh='hard',thresh_isl=isl_thresh, thresh_pix=pix_thresh,adaptive_rms_box=True)
        except RuntimeError as re:
            print(re, '\n\n')
            fitsfile = str(input("Enter full path to fits image: "))
            
            img = bdsf.process_image(fitsfile, output_opts=True,output_all=True,shapelet_do=True,solnname = base,
                                     thresh='hard',thresh_isl=isl_thresh,thresh_pix=pix_thresh,adaptive_rms_box=True)

        meanimage   = outdir+base+'-meanimage.fits'
        chan0image  = outdir+base+'-chan0image.fits'
        rmsfile     = outdir+base+'-rms.fits'
        resfile     = outdir+base+'-residuals.fits'
        modimage    = outdir+base+'-modimage.fits'
        srcfilename = outdir+base+'-source-cat.txt'
        #polimage   = outdir+'base+-polimage.fits'
        
        img.export_image(outfile =modimage,img_format = 'fits',img_type = 'gaus_model',clobber =True)
        img.export_image(outfile =rmsfile,img_format = 'fits',img_type = 'rms',clobber =True)
        img.write_catalog(outfile=srcfilename,clobber =True,format='ascii',catalog_type ='gaul')
        
        #img.export_image(outfile =resfile,img_format = 'fits',img_type = 'gaus_resid',clobber =True)
        #img.export_image(outfile =meanimage,img_format = 'fits',img_type = 'mean',clobber =True)
        #img.export_image(outfile =chan0image,img_format = 'fits',img_type = 'ch0',clobber =True)
        #img.export_image(outfile =polimage,img_format = 'fits',img_type = 'pi',clobber =True)
        #img.export_image(outfile =psf_ratio_image,img_format = 'fits',img_type = 'psf_ratio',clobber =True)

        print('\n\n done exporting bdsf.process_image result files ... Making catalog.')

        hdulist_I      = fits.open(fitsfile)
        w              = WCS(hdulist_I[0].header)
        Idata          = hdulist_I[0].data
        Imax           = Idata.max()
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

        print('done inporting source list attributes ...')
        bins = int(0.25*(sqrt(len(peak_flux)))*2)

        # filter the sources based on thresholds
        Imap= []; Ifit=[]; Qmap=[]; Umap=[]; Pmap = []; good = []#weight = []
        goodsources=0
        sources = []
        infield = base

        for i in range(len(ra)): 
            try:
                print(' > ra[%d]: %.4f, dec[%d]: %.4f'%(i,ra[i],i,dec[i]))
                pixcrd = w.wcs_world2pix(ra[i],dec[i],0,0,0)
                y      = pixcrd[0]
                x      = pixcrd[1]
                print(' > xpix[%d]: %.4f, ypix[%d]: %.4f'%(i,x,i,y))
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
                            print(' - (!) ', ie02)
                            Iamp = rms
                            pass
                    wdata = 1. #pbcorr(xp,yp,cell_deg,bmaj)
                    fluxcorr = safe_div(1.,wdata)
                    snr = peak_flux[i]/rms
                    print('***SNR: %4.2f, rms: %4.2e, Iamp: %4.2e' %(snr, rms, Iamp))
                    if(snr > snrthreshold) and (Iamp > fluxthreshold):# and (fluxcorr < fluxcorrlimit):       
                        good.append(True) 
                        goodsources = goodsources + 1
                        print("source: %5d %9.4f %7.4f %7.1f %7.1f %7.3f %9.3f %8.3f %7.2f %7.2f" \
                        % (goodsources,ra[i],dec[i],y,x,1000*rms,1000*total_flux[i],1000*peak_flux[i],\
                            Smaj_ax[i],Smin_ax[i]))
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
                    if(total_flux_err[i] < 0.0) or (isnan(dp)): # replace errors with calculated values if fit error fails
                        (df,dp,ra_err,dec_err)=source_errors(Smaj_ax[i],Smin_ax[i],PA[i],total_flux[i], peak_flux[i],rms,beam)
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
                                Smaj_ax[i],err_Smaj_ax[i],Smin_ax[i],err_Smin_ax[i],PA[i],PA_err[i], fluxcorr_mJy,y,x))
                    src_line='{:3.2f}&{:3.2f}&{:3.2f}&{:3.2f}&{:1.2e}&{:1.2e}&{:1.2e}&{:1.2e}&{:1.2e}&{:1.2e}&{:1.2e}&{:1.2e}&{:3.2f}&{:3.2f}&{:1.2e}&{:d}&{:d}\\\ \n\
                    '.format(ra[i],ra_err,dec[i],dec_err,f,df,p,dp,Smaj_ax[i],err_Smaj_ax[i],Smin_ax[i],err_Smin_ax[i],PA[i],PA_err[i],fluxcorr_mJy,int(y),int(x))
            except Exception as ex002:
                print(' (!) > ', ex002)
                pass
            
        print("\n\nFound %d good sources" % goodsources)
        print("Rejected %d sources with map peak < %7.2e and snr < %4.1f, and pblimit < %4.2f " % \
            ((len(sources) - goodsources),fluxthreshold,snrthreshold,pblimit))

        # open output files and write headers
        goodregfilename = outdir+'sources_pybdsm.reg'        # ds9 reg file of good sources
        badregfilename =  outdir+'badsources_pybdsm.reg'     # ds9 reg file of bad sources
        sourcefilename =  outdir+'sources_pybdsm.csv'        # output file for good sources
        sourcefilename_neg =  outdir+'neg_sources_pybdsm.txt'        # output file for good sources

        regfile = open(goodregfilename,'w')
        badfile = open(badregfilename,'w')
        sourcefile = open(sourcefilename,'w')

        regfile.write('# Region file format:  DS9 version 4.1 \n')
        regfile.write('global color=green width=1 \n')
        regfile.write('fk5\n')
        badfile.write('# Region file format:  DS9 version 4.1 \n')
        badfile.write('global color=red width=1 \n')
        badfile.write('fk5\n')
        sourcefile.write(  '#id,ra,ra_err,dec,dec_err,i_flux,i_err,p_flux,p_err,rms,a,a_err,b,b_err,pa,pa_err,x,y')
        sourcefile.write('\n#[],(deg),("),(deg),("),(mJy),(mJy),(mJy),("),("),(deg),(pix)')

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
                    sourceline = "\n%5d,%8.5f,%5.2f,%8.5f,%5.2f,%8.4f,%7.4f,%8.4f,%8.4f,%7.4f,%7.2f,%6.2f,%7.2f,%6.2f,%7.1f,%5.1f,%8.2f,%8.2f" % \
                    (j,sources[i].ra,sources[i].dra,sources[i].dec,sources[i].ddec,sources[i].flux,sources[i].df,sources[i].peak,sources[i].dp, \
                    sources[i].rms,sources[i].Smaj,sources[i].dmaj,sources[i].Smin,sources[i].dmin,sources[i].Spa,sources[i].dpa,\
                    sources[i].x,sources[i].y)
                    sourcefile.write(sourceline)
                else:
                    regline = "ellipse(%9.4f,%9.4f,%7.2f\",%7.2f\",%6.1f) \n" %\
                    (sources[i].ra,sources[i].dec,sources[i].Smaj,sources[i].Smin,sources[i].Spa)
                    badfile.write(regline)
            print("wrote %d sources to file %s" % (j,sourcefilename))
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
                print('(!) .. I flux may be infinite again ... check it in  {}'.format(sourcefilename))
                pass
            
            figure(figsize=(10,6))
            subplot(121)
            scatter(raerr, decerr)
            xlabel('error: RA [deg]'), ylabel('error: Dec [deg]')
            grid(ls=':', alpha=0.4); tight_layout()
            subplot(122)
            scatter(PA, PA_err)
            xlabel('PA [deg]'), ylabel('error: PA [deg]')
            grid(ls=':', alpha=0.4); tight_layout()
            savefig(outdir+'posn_errs.png')
            #show()
            
            pk_src = nanargmax(peak_flux)
            figure(figsize=(20,6))
            plot(src_id, peak_flux), plot(src_id, peak_flux, 'wo')
            plot(src_id[pk_src], peak_flux[pk_src], 'r*',\
                label=' src_id: %d \n Ipeak: %1.4e Jy/bm\n RA     : %2.4f deg\n Dec   : %2.4f deg\
            '%(src_id[pk_src], peak_flux[pk_src], ra[pk_src], dec[pk_src]) )
            legend(loc=0, fancybox=True, framealpha=0.5)
            xlabel('PyBDSM Source ID')
            ylabel('Peak Flux / [Jy/beam]')
            grid(ls=':', alpha=0.4)
            savefig(outdir+'srcs_vs_peakFlux.png')
            #show()
            
            figure(figsize=(20,6))
            subplot(131)
            hist(peak_flux, bins, alpha=0.35), xlabel('Peak Flux')
            xscale("log")
            subplot(132)
            scatter(ra, peak_flux)
            xlabel('RA/deg'), ylabel('Peak Flux')
            subplot(133)
            scatter(dec, peak_flux)
            xlabel('Dec/deg'), ylabel('Peak Flux')
            tight_layout()
            savefig(outdir+'radec_vs_peakFlux.png')
            #show()

        else:
            print('(!) *** No sources found! *** ')
        
        if show_fitted:
            img.show_fit(source_seds=True, gresid_image=False, mean_image=False,\
                        rms_image=True, ch0_flagged=True, ch0_image=True, smodel_image=True)
            show()
        
        os.chdir(pwd)
    except Exception as ex008:
        print('ERROR: ', ex008)
        print('Exiting ...')
        raise
        pass
        
    return ()

def frac_freq(f1,f2,fc):
    df = abs(f1-f2)
    ff = df/fc
    dfM= df/1.0e6 # in MHz
    return (ff, dfM)

def pyasl_separation(ra1,dec1,ra2,dec2):
    """This function returns the angular separation between two points defined by
    ra1,dec1, and ra2,dec2.   Input and output is in degrees """
    from PyAstronomy import pyasl
    pi = 3.1415928
    deg2rad = pi/180.0
    sep = pyasl.getAngDist(float(ra1), float(dec1), float(ra2), float(dec2))
    return(sep)

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
    
    print(' === Secondary sources in the field. Percent cutoff: {}%, Ipk: {:.3} Jy/bm ==='.format(perc_cut*100.,Ipeak))
    print('Source ID , Sep_from_peak[\'], RA[deg], Dec[deg], Ipk_Fraction[%%]')
    sep_Ipk = []; frc_Ipk = []; sub_src = []
    for i in range(len(srcID)):
        if Ipk[i] > perc_cut*Ipeak:
            di_pk     = separation(RAs[ nanargmax(Ipk) ], DECs[ nanargmax(Ipk) ],RAs[i], DECs[i])
            di_pk_min = di_pk*60.; frc_Ipk_i = 100.*(Ipk[i]/Ipeak)
            sep_Ipk.append(di_pk_min); frc_Ipk.append(frc_Ipk_i)
            sub_src.append([ srcID[i], RAs[i], DECs[i], Area_ellipse(maj_ax[i], min_ax[i]) ])
            print(' {} , {:.2f}  , {:.6f} , {:.6f} , {:.2f}%%'.format(srcID[i], di_pk_min, RAs[i], DECs[i], frc_Ipk_i ))
    return (array(sep_Ipk), array(frc_Ipk), array(sub_src))


def get_subbandImage(fitscuboid, scanID, imagedir,time_offset):
    """
    # extracts fits planes from fits cuboid image
    # imagedir is the common image directory for sourcefinding images
    # single scan images will have the scan ID appended to the name, if scanID provided
    - returns 3 images at low, mid and high part of the full band of the input fits
    """
    if imagedir == None:
        imagedir = os.path.dirname(fitscuboid)
    
    if scanID != None:
        fitscuboid_old = os.path.basename(fitscuboid)
        fitscuboid_new = imagedir+fitscuboid_old[:-11]+'toff%.1f-scan%d.IClean.fits'%(time_offset,scanID)
        cp_command = 'cp %s %s'%(fitscuboid, fitscuboid_new)
        print('\n\n> Renaming image according to scan and time offset via:\n%s\n\n'%cp_command)
        #raw_input('continue??')
        os.system(cp_command)
        #fitscuboid_new = imagedir+fitscuboid_new
        fitscuboid     = fitscuboid_new
    
    #hdu      = fits.open(fitscuboid)
    #dd       = hdu[0].data
    #nx       = hdu[0].header['NAXIS1']
    #ny       = hdu[0].header['NAXIS2']
    freqs_dct = {}; freqs_ar = []
    dd, hdr = fits.getdata(fitscuboid, header=True)
    nx      = hdr['NAXIS1']
    ny      = hdr['NAXIS2']
    Nn      = hdr['NAXIS']
    outdir  = pwd+'/images/'
    if not (os.path.isdir(outdir)):
        os.mkdir(outdir)
        outdir = outdir+'/'
    
    low_bim  = dd[0,1]      # low band image at freqs ~ [8.989462890625E+08; 9.658212890625E+08]
    mid_bim  = dd[0,5]      # mid band image at freqs ~ [1.279088867188E+09; 1.367907226562E+09]
    hgh_bim  = dd[0,7]      # high band image at freqs ~ [1.465084960938E+09; 1.567278320312E+09]
    avg_wim  = dd[0,0]      # full band image
    slct_im  = dd[0,12]
    
    low_freq = 'FREQ0004'        # low band image at freqs ~ [8.989462890625E+08; 9.658212890625E+08]
    mid_freq = 'FREQ0014'        # mid band image at freqs ~ [1.279088867188E+09; 1.367907226562E+09]
    hgh_freq = 'FREH0021'        # high band image at freqs ~ [1.465084960938E+09; 1.567278320312E+09]
    avg_freq = 'CRVAL3'          # full band image
    slct_freq= 'FREQ0003'
    
    hdu_low  = fits.PrimaryHDU(data=low_bim)
    hdu_mid  = fits.PrimaryHDU(data=mid_bim)
    hdu_hgh  = fits.PrimaryHDU(data=hgh_bim)
    hdu_wide = fits.PrimaryHDU(data=avg_wim)
    hdu_slct = fits.PrimaryHDU(data=slct_im)
    HDUs = [[hdu_low,low_freq],[hdu_mid,mid_freq],[hdu_hgh,hgh_freq],\
        [hdu_wide,avg_freq],[hdu_slct,slct_freq]]
    
    try:
        for hdr_i in hdr:
            if hdr_i.startswith('FREQ'):
                freqs_dct[hdr_i] = hdr[hdr_i]
                freqs_ar.append(hdr[hdr_i])
    except Exception as ex006:
        print(ex006)
        pass
    
    for HDU_i in HDUs:
        try:
            HDU_i[0].header.set('NAXIS', hdr['NAXIS'])
            k_frq = 3; k_frq_key = HDU_i[1]
            for k in range(1, Nn+1):
                HDU_i[0].header.set('NAXIS%d'%k, nx)
                HDU_i[0].header.set('NAXIS%d'%k, ny)
                
                if hdr['CTYPE%d'%k] == 'SPECLNMF':
                    HDU_i[0].header.set('CTYPE%d'%k, 'FREQ', 'Axis type')
                else:
                    HDU_i[0].header.set('CTYPE%d'%k, hdr['CTYPE%d'%k], 'Axis type')
                
                HDU_i[0].header.set('CDELT%d'%k, hdr['CDELT%d'%k], 'Axis coordinate increment')
                HDU_i[0].header.set('CRPIX%d'%k, hdr['CRPIX%d'%k], 'Axis coordinate reference pixel')
                HDU_i[0].header.set('CROTA%d'%k, hdr['CROTA%d'%k], 'Axis coordinate rotation')
                HDU_i[0].header.set('CRVAL%d'%k, hdr['CRVAL%d'%k], 'Axis coordinate value at CRPIX')
            
            hdr_frq = hdr[k_frq_key]
            HDU_i[0].header.set('CRVAL%d'%k_frq, hdr_frq, 'Axis coordinate value at CRPIX')
            HDU_i[0].header.set('OBSRA' , hdr['OBSRA'] , 'Observed Right Ascension')
            HDU_i[0].header.set('OBSDEC', hdr['OBSDEC'], 'Observed declination')
            HDU_i[0].header.set('EPOCH', hdr['EPOCH'], 'Celestial coordiate epoch')
            HDU_i[0].header.set('EQUINOX', hdr['EQUINOX'], 'Celestial coordiate equinox')
            HDU_i[0].header.set('BUNIT', hdr['BUNIT'], 'Image pixel units')
            HDU_i[0].header.set('VELREF', hdr['VELREF'], '>256 radio, 1 LSR, 2 Hel, 3 Obs')
            HDU_i[0].header.set('CLEANBMJ', hdr['CLEANBMJ'], 'Convolving Gaussian major axis FWHM (deg)')
            HDU_i[0].header.set('CLEANBMN', hdr['CLEANBMN'], 'Convolving Gaussian minor axis FWHM (deg)')
            HDU_i[0].header.set('CLEANBPA', hdr['CLEANBPA'], 'Convolving Gaussian position angle (deg)')
            HDU_i[0].header.set('BMAJ', hdr['CLEANBMJ'], 'Clean/Convolving Gaussian major axis FWHM (deg)')
            HDU_i[0].header.set('BMIN', hdr['CLEANBMN'], 'Clean/Convolving Gaussian minor axis FWHM (deg)')
            HDU_i[0].header.set('BPA', hdr['CLEANBPA'], 'Clean/Convolving Gaussian position angle (deg)')
            HDU_i[0].header.set('CLEANNIT', hdr['CLEANNIT'], 'Number of Clean iterations')
            HDU_i[0].header.set('XPXOFF', hdr['XPXOFF'], 'x pixel offset')
            HDU_i[0].header.set('YPXOFF', hdr['YPXOFF'], 'y pixel offset')
            HDU_i[0].header.set('ALPHA', hdr['ALPHA'])
            #HDU_i[0].header.set('RFALPHA', hdr['RFALPHA'])
            HDU_i[0].header.set('ORIGIN', hdr['ORIGIN'], 'Software last writing file')
            HDU_i[0].header.set('TELESCOP', hdr['TELESCOP'], 'Telescope used')
            
            try:
                for hdr_i in hdr:
                    if hdr_i.startswith('FRE') and hdr_i.endswith(k_frq_key[-4:]):
                        HDU_i[0].header.set(hdr_i, hdr[hdr_i])
            except Exception as ex006:
                print(ex006)
                pass
            
            good_histr = ''
            try:
                for hdr_j in hdr:
                    if hdr_j.startswith('HISTORY'):
                        good_histr = good_histr+str(hdr[hdr_j]).encode('ascii')
                        HDU_i[0].header.add_history( good_histr.encode('ascii') )
            except Exception as ex006:
                print(ex006)
                pass
            
            HDU_i[0].header.set('HISTORY','')
            #HDU_i[0].header.set('', hdr[''], '')
        except KeyError as ex005:
            print(ex005)
            print(' > Writing out the BASIC header ONLY')
            HDU_i[0].header.set('NAXIS', hdr['NAXIS'])
            HDU_i[0].header.set('NAXIS1', nx)
            HDU_i[0].header.set('NAXIS2', ny)
            HDU_i[0].header.set('CTYPE1', hdr['CTYPE1'], 'Axis type')
            HDU_i[0].header.set('CDELT1', hdr['CDELT1'], 'Axis coordinate increment')
            HDU_i[0].header.set('CRPIX1', hdr['CRPIX1'], 'Axis coordinate reference pixel')
            HDU_i[0].header.set('CROTA1', hdr['CROTA1'], 'Axis coordinate rotation')
            HDU_i[0].header.set('CRVAL1', hdr['CRVAL1'], 'Axis coordinate value at CRPIX')
            HDU_i[0].header.set('CTYPE2', hdr['CTYPE2'], 'Axis type')
            HDU_i[0].header.set('CDELT2', hdr['CDELT2'], 'Axis coordinate increment')
            HDU_i[0].header.set('CRPIX2', hdr['CRPIX2'], 'Axis coordinate reference pixel')
            HDU_i[0].header.set('CROTA2', hdr['CROTA2'], 'Axis coordinate rotation')
            HDU_i[0].header.set('CRVAL2', hdr['CRVAL2'], 'Axis coordinate value at CRPIX')
            HDU_i[0].header.set('HISTORY', 'BASIC header ONLY - check original image (%s) header for full details'%fitscuboid)
            pass
    
    low_fileout = outdir+os.path.basename(fitscuboid)[:-5]+'.lowband.fits'
    hdu_low.writeto(low_fileout, overwrite=True)
    
    print('\nlow band file   : ', low_fileout)
    
    mid_fileout = outdir+os.path.basename(fitscuboid)[:-5]+'.midband.fits'
    hdu_mid.writeto(mid_fileout, overwrite=True)
    print('mid band file   : ', mid_fileout)
    
    hgh_fileout = outdir+os.path.basename(fitscuboid)[:-5]+'.hghband.fits'
    hdu_hgh.writeto(hgh_fileout, overwrite=True)
    print('high band file  : ', hgh_fileout)
    
    wide_fileout = outdir+os.path.basename(fitscuboid)[:-5]+'.wideband.fits'
    hdu_wide.writeto(wide_fileout, overwrite=True)
    print('wide band file  : ', wide_fileout, '\n\n')
    
    slct_fileout = outdir+os.path.basename(fitscuboid)[:-5]+'.freq%.3fGHz.fits' %(hdr_frq/(1.e9))
    hdu_slct.writeto(slct_fileout, overwrite=True)
    print(' > selected band centre: %.3f GHz'%(hdr_frq/(1.e9)))
    print(' > selcted band file   : ', slct_fileout, '\n\n')
    
    #print 'HEADER for: ', wide_fileout
    #print hdu_wide[0].header
    
    fitsdir = os.path.dirname(fitscuboid)
    #os.system('cp %s %s/lowband/'%(low_fileout, imagedir))
    #os.system('cp %s %s/midband/'%(mid_fileout, imagedir))
    #os.system('cp %s %s/highband/'%(hgh_fileout, imagedir))
    #os.system('cp %s %s/wideband/'%(wide_fileout, imagedir))
    
    return (hdr_frq)

def scaled_rms(im_rms, ants,Ttot):
    return(im_rms/sqrt(ants(ants-1.)*Ttot ))


def channelised_imstat(fitsfile, base=None, algo=None):
    from astropy.io import fits
    from matplotlib import pyplot as plt
    import os, time
    import numpy as np
    
    """
    Image statistics in a channelised fashion.
    - Runs in CASA session on a fits file
    - uses CASA's imstat()
    - input: 
        fitsfile: fits cube image, tested on one produced by katsdpimager/Obit MFImage
        base    : human readable reference 
        algo    : alpgorythm name to use in CASA's imstat() task (https://casa.nrao.edu/docs/TaskRef/imstat-task.html)
    output:
        stats0: statistics dictionary from imstat()
        freqs: frequencies as reported by imstat()
        freq0: frequencies as read in from fitsfile header
        fd_max: max flux density
        fd_avg: mean flux density
        fd_rms: flux density rms
        fd_sgm: flux density standard deviation
    """
    
    if base == None:
        base = os.path.basename(fitsfile)[:10]
    
    if algo == None:
        algo = 'classic'
        
    if base == '1590555655':#????
        mode = "c544M32kWB8s"
    elif base == '1590547878':#
        mode = "c544M4kWB8s"
    elif base == '1590381727':#
        mode = "c856M4k_n107M8s"
    elif base == '1590377455':#
        mode = "c856M32kWB8s"
    elif base == '1590129959':#
        mode = "c856M4kWB8s"
    
    elif "1585989326" in base:
        mode = "c8561kWB-April"
    elif "1586305859" in base:
        mode = "c8564kWB-April"
    elif "1586070058" in base:
        mode = "c5444kWB-April"
    elif "1586162571" in base:
        mode = "c5441kWB-April"
    
    else:
        mode = str(input("CBID is %s. Enter Observation mode (string) : "%base))
    
    show_figs      = True
    timestr        = time.strftime("%Y%m%d-%H%M%S")
    outfile        = open(fitsfile[:-4]+mode+'-'+algo+'-channelised_imstats.txt','w+')
    outfile2       = open(fitsfile[:-4]+mode+'-'+algo+'-channelised_imstats.csv','w')
    
    dd, hdr = fits.getdata(fitsfile, header=True)
    nplanes = hdr['NAXIS3']
    nx      = hdr['NAXIS1']
    ny      = hdr['NAXIS2']
    N0      = hdr['NAXIS']
    x0pix   = hdr['CRPIX1']
    y0pix   = hdr['CRPIX2']
    
    freqs, fd_max, fd_min, fd_avg, fd_rms, fd_sgm = [],[],[],[],[],[]
    stats0, freq0 = [], []
    
    outfile2.write('#mode,freqs,freq0,fd_max,fd_min,fd_avg,fd_rms,fd_sgm\n')
    for p in range(2,nplanes):
        pi     = p-1
        freqL  = hdr['FREL{:04d}'.format(pi)]
        freqC  = hdr['FREQ{:04d}'.format(pi)]
        freqH  = hdr['FREH{:04d}'.format(pi)]
        dfreq  = freqH - freqL
        freqCGHz = freqC/1.0e+9
        p_out  = fitsfile[:-4]+'FREQc_{:.3f}GHz.fits'.format(freqCGHz)
        p_data = fits.PrimaryHDU(data=dd[0,p])
        
        for k in range(1, N0+1):
                p_data.header.set('CTYPE%d'%k, hdr['CTYPE%d'%k], 'Axis type')
                p_data.header.set('CDELT%d'%k, hdr['CDELT%d'%k], 'Axis coordinate increment')
                p_data.header.set('CRPIX%d'%k, hdr['CRPIX%d'%k], 'Axis coordinate reference pixel')
                p_data.header.set('CROTA%d'%k, hdr['CROTA%d'%k], 'Axis coordinate rotation')
                p_data.header.set('CRVAL%d'%k, hdr['CRVAL%d'%k], 'Axis coordinate value at CRPIX')
        
        p_data.header.set('OBSRA' , hdr['OBSRA'] , 'Observed Right Ascension')
        p_data.header.set('OBSDEC', hdr['OBSDEC'], 'Observed declination')
        p_data.header.set('EPOCH', hdr['EPOCH'], 'Celestial coordiate epoch')
        p_data.header.set('EQUINOX', hdr['EQUINOX'], 'Celestial coordiate equinox')
        p_data.header.set('BUNIT', hdr['BUNIT'], 'Image pixel units')
        p_data.header.set('VELREF', hdr['VELREF'], '>256 radio, 1 LSR, 2 Hel, 3 Obs')
        p_data.header.set('CRVAL3', freqC, 'frequency of image plane (Hz)')
        p_data.header.set('CDELT3', dfreq, 'frequency span of image (Hz)')
        p_data.header.set('CLEANBMN', hdr['CLEANBMN'], 'Convolving Gaussian minor axis FWHM (deg)')
        p_data.header.set('CLEANBMJ', hdr['CLEANBMJ'], 'Convolving Gaussian major axis FWHM (deg)')
        p_data.header.set('CLEANBPA', hdr['CLEANBPA'], 'Convolving Gaussian position angle (deg)')
        p_data.header.set('OBJECT', hdr['OBJECT'], 'Name of object/target')
        p_data.writeto(p_out, overwrite=True)
        
        try:
            si      = imstat(imagename=p_out, algorithm=algo)
            si2     = imstat(imagename=p_out, algorithm=algo, )
            frq_si  = float(si['blcf'][-13:-3])
            fdmx_si = si['max'][0]
            fdmn_si = si['min'][0]
            fdavg_si= si['mean'][0]
            fdrms_si= si['rms'][0]
            fdsgm_si= si['sigma'][0]
            
            stats0.append( {str(p_out): si} )
            freqs.append( frq_si ), fd_max.append( fdmx_si ),fd_min.append( fdmn_si ) , fd_avg.append( fdavg_si )
            fd_rms.append( fdrms_si ), fd_sgm.append( fdsgm_si ), freq0.append(freqC)
            
            outfile.write("\nFull stats:\nFile: {}\nStats: {}\n".format(p_out, si))
            outfile2.write('{},{},{},{},{},{},{},{}\n'.format(mode, frq_si, freqC, fdmx_si, fdmn_si, fdavg_si, fdrms_si, fdsgm_si) )
        except Exception as ex:
            print (' ! Exception ! : ', ex)
            pass
    
    freqs, freq0, fd_max   = np.array(freqs),np.array(freq0),np.array(fd_max)
    fd_avg, fd_rms, fd_sgm = np.array(fd_avg),np.array(fd_rms),np.array(fd_sgm)    
    
    plt.figure()
    plt.title('Mode: {}, CBID: {}'.format(mode,base))
    plt.scatter(freq0, fd_max)
    plt.ylabel('max [Jy/bm]')
    plt.xlabel('Frequency [Hz]')
    plt.grid(True, ls = ':')
    plt.savefig(fitsfile[:-5]+'-imstat-Imax.png')
    
    if not show_figs:
        plt.close()
    
    plt.figure()
    plt.title('Mode: {}, CBID: {}'.format(mode,base))
    plt.scatter(freq0, fd_avg)
    plt.ylabel('mean [Jy/bm]')
    plt.xlabel('Frequency [Hz]')
    plt.grid(True, ls = ':')
    plt.savefig(fitsfile[:-5]+'-imstat-Iavg.png')
    
    if not show_figs:
        plt.close()
    
    plt.figure()
    plt.title('Mode: {}, CBID: {}'.format(mode,base))
    plt.scatter(freq0, fd_rms)
    plt.ylabel('rms [Jy/bm]')
    plt.xlabel('Frequency [Hz]')
    plt.grid(True, ls = ':')
    plt.savefig(fitsfile[:-5]+'-imstat-Irms.png')
    
    if not show_figs:
        plt.close()
    
    if show_figs:
        plt.show()
    
    outfile.close()
    outfile2.close()
    
    return (stats0, freqs, freq0, fd_min, fd_max, fd_avg, fd_rms, fd_sgm)


def generic_get_subbandImage(fitscuboid, slct_im_plane=None):
    """
    - extracts fits planes from fits cuboid image produced by the katsdpimager
    
    Inputs:
    - fitscuboid    : string, fits image from katsdpimager: "katsdpcontim"
    - slct_im_plane : int, user specified frequency plane in the fitscuboid
    Returns:
    - images at the low, mid and high parts of the full band of the input fits 
      as well as the MFS and image and an image at a specified frequency.
    """
    
    freqs_dct = {}; freqs_ar = []
    dd, hdr = fits.getdata(fitscuboid, header=True)
    nx      = hdr['NAXIS1']
    ny      = hdr['NAXIS2']
    nz      = hdr['NAXIS3']
    Nn      = hdr['NAXIS']
    outdir  = pwd+'/images/'
    if not (os.path.isdir(outdir)):
        os.mkdir(outdir)
        outdir = outdir+'/'
    
    low_bim  = dd[0,1]      # low band image at freqs ~ [8.989462890625E+08; 9.658212890625E+08]
    mid_bim  = dd[0,5]      # mid band image at freqs ~ [1.279088867188E+09; 1.367907226562E+09]
    hgh_bim  = dd[0,7]      # high band image at freqs ~ [1.465084960938E+09; 1.567278320312E+09]
    avg_wim  = dd[0,0]      # full band image
    if slct_im_plane == None:
        slct_im  = dd[0,nz-1]
    else:
        slct_im = slct_im_plane
        
    low_freq = 'FREQ0004'        # low band image at freqs ~ [8.989462890625E+08; 9.658212890625E+08]
    mid_freq = 'FREQ0014'        # mid band image at freqs ~ [1.279088867188E+09; 1.367907226562E+09]
    hgh_freq = 'FREH0021'        # high band image at freqs ~ [1.465084960938E+09; 1.567278320312E+09]
    avg_freq = 'CRVAL3'          # full band image
    slct_freq= 'FREQ0003'
    
    hdu_low  = fits.PrimaryHDU(data=low_bim)
    hdu_mid  = fits.PrimaryHDU(data=mid_bim)
    hdu_hgh  = fits.PrimaryHDU(data=hgh_bim)
    hdu_wide = fits.PrimaryHDU(data=avg_wim)
    hdu_slct = fits.PrimaryHDU(data=slct_im)
    HDUs = [[hdu_low,low_freq],[hdu_mid,mid_freq],[hdu_hgh,hgh_freq],\
        [hdu_wide,avg_freq],[hdu_slct,slct_freq]]
    
    try:
        for hdr_i in hdr:
            if hdr_i.startswith('FREQ'):
                freqs_dct[hdr_i] = hdr[hdr_i]
                freqs_ar.append(hdr[hdr_i])
    except Exception as ex006:
        print(ex006)
        pass
    
    for HDU_i in HDUs:
        try:
            HDU_i[0].header.set('NAXIS', hdr['NAXIS'])
            k_frq = 3; k_frq_key = HDU_i[1]
            for k in range(1, Nn+1):
                HDU_i[0].header.set('NAXIS%d'%k, nx)
                HDU_i[0].header.set('NAXIS%d'%k, ny)
                
                if hdr['CTYPE%d'%k] == 'SPECLNMF':
                    HDU_i[0].header.set('CTYPE%d'%k, 'FREQ', 'Axis type')
                else:
                    HDU_i[0].header.set('CTYPE%d'%k, hdr['CTYPE%d'%k], 'Axis type')
                
                HDU_i[0].header.set('CDELT%d'%k, hdr['CDELT%d'%k], 'Axis coordinate increment')
                HDU_i[0].header.set('CRPIX%d'%k, hdr['CRPIX%d'%k], 'Axis coordinate reference pixel')
                HDU_i[0].header.set('CROTA%d'%k, hdr['CROTA%d'%k], 'Axis coordinate rotation')
                HDU_i[0].header.set('CRVAL%d'%k, hdr['CRVAL%d'%k], 'Axis coordinate value at CRPIX')
            
            hdr_frq = hdr[k_frq_key]
            HDU_i[0].header.set('CRVAL%d'%k_frq, hdr_frq, 'Axis coordinate value at CRPIX')
            HDU_i[0].header.set('OBSRA' , hdr['OBSRA'] , 'Observed Right Ascension')
            HDU_i[0].header.set('OBSDEC', hdr['OBSDEC'], 'Observed declination')
            HDU_i[0].header.set('EPOCH', hdr['EPOCH'], 'Celestial coordiate epoch')
            HDU_i[0].header.set('EQUINOX', hdr['EQUINOX'], 'Celestial coordiate equinox')
            HDU_i[0].header.set('BUNIT', hdr['BUNIT'], 'Image pixel units')
            HDU_i[0].header.set('VELREF', hdr['VELREF'], '>256 radio, 1 LSR, 2 Hel, 3 Obs')
            HDU_i[0].header.set('CLEANBMJ', hdr['CLEANBMJ'], 'Convolving Gaussian major axis FWHM (deg)')
            HDU_i[0].header.set('CLEANBMN', hdr['CLEANBMN'], 'Convolving Gaussian minor axis FWHM (deg)')
            HDU_i[0].header.set('CLEANBPA', hdr['CLEANBPA'], 'Convolving Gaussian position angle (deg)')
            HDU_i[0].header.set('BMAJ', hdr['CLEANBMJ'], 'Clean/Convolving Gaussian major axis FWHM (deg)')
            HDU_i[0].header.set('BMIN', hdr['CLEANBMN'], 'Clean/Convolving Gaussian minor axis FWHM (deg)')
            HDU_i[0].header.set('BPA', hdr['CLEANBPA'], 'Clean/Convolving Gaussian position angle (deg)')
            try:
                HDU_i[0].header.set('CLEANNIT', hdr['CLEANNIT'], 'Number of Clean iterations')
            except Exception as ex0007:
                print('\n\n',ex0007)
                pass
            HDU_i[0].header.set('XPXOFF', hdr['XPXOFF'], 'x pixel offset')
            HDU_i[0].header.set('YPXOFF', hdr['YPXOFF'], 'y pixel offset')
            HDU_i[0].header.set('ALPHA', hdr['ALPHA'])
            #HDU_i[0].header.set('RFALPHA', hdr['RFALPHA'])
            HDU_i[0].header.set('ORIGIN', hdr['ORIGIN'], 'Software last writing file')
            HDU_i[0].header.set('TELESCOP', hdr['TELESCOP'], 'Telescope used')
            
            try:
                for hdr_i in hdr:
                    if hdr_i.startswith('FRE') and hdr_i.endswith(k_frq_key[-4:]):
                        HDU_i[0].header.set(hdr_i, hdr[hdr_i])
            except Exception as ex006:
                print(ex006)
                pass
            
            good_histr = ''
            try:
                for hdr_j in hdr:
                    if hdr_j.startswith('HISTORY'):
                        good_histr = good_histr+str(hdr[hdr_j]).encode('ascii')
                        HDU_i[0].header.add_history( good_histr.encode('ascii') )
            except Exception as ex006:
                print(ex006)
                pass
            
            HDU_i[0].header.set('HISTORY','')
            #HDU_i[0].header.set('', hdr[''], '')
        except KeyError as ex005:
            print(ex005)
            print(' > Writing out the BASIC header ONLY')
            HDU_i[0].header.set('NAXIS', hdr['NAXIS'])
            HDU_i[0].header.set('NAXIS1', nx)
            HDU_i[0].header.set('NAXIS2', ny)
            HDU_i[0].header.set('CTYPE1', hdr['CTYPE1'], 'Axis type')
            HDU_i[0].header.set('CDELT1', hdr['CDELT1'], 'Axis coordinate increment')
            HDU_i[0].header.set('CRPIX1', hdr['CRPIX1'], 'Axis coordinate reference pixel')
            HDU_i[0].header.set('CROTA1', hdr['CROTA1'], 'Axis coordinate rotation')
            HDU_i[0].header.set('CRVAL1', hdr['CRVAL1'], 'Axis coordinate value at CRPIX')
            HDU_i[0].header.set('CTYPE2', hdr['CTYPE2'], 'Axis type')
            HDU_i[0].header.set('CDELT2', hdr['CDELT2'], 'Axis coordinate increment')
            HDU_i[0].header.set('CRPIX2', hdr['CRPIX2'], 'Axis coordinate reference pixel')
            HDU_i[0].header.set('CROTA2', hdr['CROTA2'], 'Axis coordinate rotation')
            HDU_i[0].header.set('CRVAL2', hdr['CRVAL2'], 'Axis coordinate value at CRPIX')
            HDU_i[0].header.set('HISTORY', 'BASIC header ONLY - check original image (%s) header for full details'%fitscuboid)
            pass
    
    low_fileout = outdir+os.path.basename(fitscuboid)[:-5]+'.lowband.fits'
    hdu_low.writeto(low_fileout, overwrite=True)
    
    print('\nlow band file   : ', low_fileout)
    
    mid_fileout = outdir+os.path.basename(fitscuboid)[:-5]+'.midband.fits'
    hdu_mid.writeto(mid_fileout, overwrite=True)
    print('mid band file   : ', mid_fileout)
    
    hgh_fileout = outdir+os.path.basename(fitscuboid)[:-5]+'.hghband.fits'
    hdu_hgh.writeto(hgh_fileout, overwrite=True)
    print('high band file  : ', hgh_fileout)
    
    wide_fileout = outdir+os.path.basename(fitscuboid)[:-5]+'.wideband.fits'
    hdu_wide.writeto(wide_fileout, overwrite=True)
    print('wide band file  : ', wide_fileout, '\n\n')
    
    slct_fileout = outdir+os.path.basename(fitscuboid)[:-5]+'.freq%.3fGHz.fits' %(hdr_frq/(1.e9))
    hdu_slct.writeto(slct_fileout, overwrite=True)
    print(' > selected band centre: %.3f GHz'%(hdr_frq/(1.e9)))
    print(' > selcted band file   : ', slct_fileout, '\n\n')
    
    #fitsdir = os.path.dirname(fitscuboid)
    print(' > Output fits in: ', outdir)
    return ()


def MeerKAT_rms(Ants, Df, T0, SEFD, S_peak):
    """
    # Theoretical thermal Stokes I rms noise in Jy on MeerKAT (64, 13m Antennas)
    # from: https://drive.google.com/file/d/16LXJQ8KXUqWx2lo7GFm-fSIUfVL41w1k/view
    
    input:
    - Ants   : number of antennas
    - Df     : bandwidth for noise calculation
    - T0     : total integration time
    - SEFD   : system equivalent flux density (Jy), defined as the flux density of a radio source that doubles the system temperature.
    - S_peak : Known Peak flux (Jy)
    """
    
    if Df == None:
        Df = 856.e6     #MeerKAT full bandwidth in Hz
    if SEFD == None:
        SEFD  = 443.5
    
    NN_1  = Ants*(Ants-1.)
    S_rms = SEFD/sqrt( 2.*NN_1*Df*T0 )
    SNR   = S_peak/S_rms
    Dt0   = ( ( SEFD/(S_peak/10.) )**2 )/( Df*NN_1 )   # the integration time to get SNR = 10
    print('Given S_peak: {} Jy, {} antennas, SEFD: {} Jy, bandwith: {1.3e} Hz and {}s of integration:'.format(S_peak,Ants,SEFD,Df,T0))
    print('    I_rms   = {:1.4e} Jy, SNR = {:.2f}'.format(S_rms, SNR))
    print('    For SNR = 10, integration time is {:1.4e}s'.format(Dt0))
    RMSs  = {'I_rms': S_rms, 'SNR': SNR}
    return (RMSs)

def RMS_simple(num):
    return sqrt(sum(n*n for n in num)/len(num))

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
        freq_col = str(input(' Enter frequency columnn name : '))
        fmax_col = str(input(' Enter flux peak columnn name : '))
        frms_col = str(input(' Enter flux RMS columnn name  : '))
    
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

def MKCosBeam(RA0, DEC0, RA, DEC, nu):
    """
    Calculate cosine beam shape (Condon & Ransom, Essential Radio Astronomy eq 3.95) at 
    the position of a field source given a frequency 
   
    Return power gain of circularly symmetric beam
    * RA0, DEC0 = float(), float(), Phase center RA & DEC (degrees)
    e.g. RA0, DEC0 = 63.36, -80.0 deg for DEEP2
    * RA, DEC   = float(), 
    float(), Input RA & DEC (degrees) of field source
    * nu        = float(), Frequency (Hz)
    """
    ################################################################
    from math import radians, cos, pi
    import numpy as np
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    
    c1   = SkyCoord(RA*u.deg, DEC*u.deg, frame='icrs')
    c0   = SkyCoord(RA0*u.deg, DEC0*u.deg, frame='fk5')
    sep  = c1.separation(c0)
    rho  = sep.deg 
    
    #theta_b = radians(57.5/60) * (1.5e9/nu)
    theta_b = 0.0167261 * (1.5e9/nu)
    rhor = 1.18896*radians(rho)/theta_b
    gain = (cos(pi*rhor)/(1.-4.*(rhor**2)))**2
    
    #print(f' > Sep  : {rho} deg')
    #print(f'  - gain: {gain} deg')
    
    return gain
    
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
    from matplotlib import gridspec
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
    
    pbdir = str(abs_filedir)+f"/pbcorr"
    if not os.path.isdir(pbdir):
        os.mkdir(pbdir)
        
    RA0, DEC0 = 63.36, -80.0 #phase center in deg for DEEP2
    
    cat_out = str(outname)+"_sources_summary.csv"
    catalog = open(str(abs_filedir)+f"/{cat_out}",'w')
    catalog.write('#RA_deg,DECdeg,Imed,Istd,Alpha,Vmed,Vstd,a,b,c\n')
    for f0 in flist:
        t    = atpy.Table(f0)
        base = os.path.basename(f0).rstrip('_IQU_perchan.tbl')
        if prefix:
            base = base.lstrip(prefix)
        try:
            RA  = float(base.split('-',1)[0])
            DEC = -1.*float(base.split('-',1)[1])
        except ValueError as voe000:
            print(f'\n****\n >> ValueError: {voe000}')
            print(f'  > Input file: {f0}\n****\n')
            RA  = '--'
            DEC = '--'
        
        SI0  = t[ycol]
        freq = np.array(t[xcol])[np.where(SI0*2 > SI0)]
        SI   = np.array(t[ycol])[np.where(SI0*2 > SI0)]
        SV   = np.array(t[zcol])[np.where(SI0*2 > SI0)]
        
        pb_g = []
        for freqi in freq:
            pb_g.append( MKCosBeam(RA0, DEC0, RA, DEC, freqi) )
        pb_g = np.array(pb_g)
        
        try:
            logf     = np.log10(freq)
            logS     = np.log10(SI)
            logScorr = np.log10(SI/pb_g)
        except TypeError as ex002:
            t.describe()
            print(f' >> TypeError: {ex002}')
            print(f' >> freqs: {freq}')
            print(f' >> SI   : {SI}\n\n\n')
        
            logf = [], logS = [], 
            for i in range(len(freq)):
                logf.append(np.log10(freq[i]))
                logS.append(np.log10(SI[i]))
                logScorr.append( np.log10(SI[i]/pb_g[i]) )
            logf     = np.array(logf)
            logS     = np.array(logS)
            logScorr = np.array(logScorr)
            pass
        
        print(f'Source at RA, DEC: {RA}, {DEC}\n -- logScorr: {logScorr}\n')
        
        fig = plt.figure()
        gs = gridspec.GridSpec(2, 1, height_ratios = [3, 1]) 
        ax0 = plt.subplot(gs[0])
        ax0.set_title(f'Source at RA, DEC: {RA}, {DEC} deg')
        ax0.scatter(logf, logScorr, label='log$_{10}$(S/cos pb)')
        ax0.scatter(logf, logS)
        ax0.set_ylabel('log S')
        ax0.grid(ls=":"); ax0.legend(loc='best')
        
        ax1 = plt.subplot(gs[1])
        ax1.plot(logf, pb_g, label='cos pb')
        ax1.grid(ls=":"); ax1.legend(loc='best')
        ax1.set_ylabel('pb gain'); ax1.set_xlabel('log freq [Hz]')
        plt.tight_layout()
        plt.savefig(pbdir+f"/{str(RA)}{str(DEC)}_cosinepb_corr.png")
        plt.close()
        
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
    Runs in pure python
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


def CrossMatch(RA1_arr, DEC1_arr, RA2_arr, DEC2_arr):
    from astropy.coordinates import SkyCoord
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


def pos_varCMC1xCMC2(ra1, dec1, x1, y1, ra2, dec2, x2, y2, base1, base2,x0=None, y0=None):
    """
    Compare cross-matched source positions (RA and Dec) from different catalogs.
    - tested on the same image from cmc1 and cmc2 
    - also analyses sky rotation between the two coordinates
    
    input:
    - ra, dec, x, y need to be arrays or lists
    - x0, y0 (floats): pixel reference points for rotation analyses.
    - basei is a human readable string that defines the data in cmci.
    
    Returns:
    - Positions offset plots
    """
    
    n1 = list(range(len(ra1))); n2 = list(range(len(ra2)))
    
    if x0 == None or y0 == None:
        x0, y0 = 0.,0.
    
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
    
    phi_d  = array(phi_d)
    phi2_d = array(phi2_d)
    theta  = array(theta)/cbmaj
    
    
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

def im_rot(xfile,base,c0,cf,cref):
    """
    - ORIGINAL Code by Niruj Ramanujam (nramanujam@ska.ac.za)
    - Orignal doc string:
    L1k - 1585989326 imsize: 5799 pix
    L4kn- 1586305859 imsize: 5962 pix
    U1k - 1586162571 imsize: 3684 pix
    U4k - 1586070058 imsize: 3696 pix
    L4k - 1586074843 imsize: 5804 pix
    
    # approx, is valid only at equator
    # theta is PA of src2 from src1
    # phi is PA of src from image centre
    # no rotation => random; rotation => linear
    
    New doc string:
    
    Input:
    - xfile : file containing X-matched catalogs
    - base  : human readable text reference for the whole dataset
    
    Parameters (defined in function) to change/update:
    - sizes (image sizes in pix)
    - cell (resolution)
    - c0, cf (tags for 1st and last catalog parameter subscript, assuming all 
    corresponding columns are only differentiated by a numeric subscript)
    - cref (tag for reference table)
    """

    import atpy, os, sys, glob
    import numpy as N
    import pylab as pl
    from astropy.io import fits
    #pl.ion()

    f     = xfile #"L1K-x-U4K-x-L4Kv2-x-U1K-x-SUMSS.tbl"
    outdir= os.path.dirname(xfile)
    
    sizes = []
    cell  = []
    
    #sizes = N.asarray([5803, 3707, 5962, 3685, 5799])  # repeat l1k since l1k is ref since sumss has fat beam
    #cell = N.asarray([0.0003448295, 0.0005410507, 0.0003354017, 0.0005428091, 0.0003448295])
    #cell *= 3600

    t  = atpy.Table(f)
    #t.describe()  # will show column names
    b_tags,bases = [], [] #file metadata from input ATPY table (cross match table)
    with open(xfile, 'r') as xf:
        for l in xf:
            if 'INIMAGE' in l or 'FREQ0' in l or 'COMMENT' in l: 
                b_tags.append(l)
                print(f' >> l: {l}')
    
    print('\n\n')
    for k in range(len(b_tags)):
        try:
            obsID   = b_tags[k].split('"',1)[1][:10]
            refFreq = b_tags[k+1].split('=',1)[1].rstrip('\n')
            naxis1  = b_tags[k+2].split('NAXIS1',1)[1][3:8]
            naxis2  = b_tags[k+2].split('NAXIS2',1)[1][3:8]
            cellsize= abs(float(b_tags[k+2].split('CDELT1',1)[1][3:].split('"',1)[0].rstrip('CR')))
            sizes.append(int(naxis1)); cell.append(cellsize*3600.)
            print(f' > obs: {obsID}, freq: {refFreq} Hz, pix: {naxis1} x {naxis2}')
            bases.append([obsID,refFreq]) 
            print(f' >>> file {k}: {b_tags[k]}')
        except Exception as ex002: 
            print(f' Exception >>> {ex002}')
            try:
                if float(obsID):
                    refFreq = float(input(f"Reference frequency (in Hz) of obsID {obsID}: "))
                    naxis1  = float(input(f"'NAXIS1' value from obsID {obsID} image header: "))
                    naxis2  = float(input(f"'NAXIS2' value from obsID {obsID} image header: "))
                    cellsize= float(input(f"'CDELT1' abs.value from obsID {obsID} image header: "))
                    
                    sizes.append(int(naxis1)); cell.append(cellsize*3600.)
                    print(f' > obs: {obsID}, freq: {refFreq} Hz, pix: {naxis1} x {naxis2}')
                    bases.append([obsID,refFreq]) 
                    print(f' >>> file {k}: {b_tags[k]}')
            except Exception as ex001:
                print(f'{ex001}')
                pass
            
            pass 
    print('\n\n')
    
    sizes = N.array(sizes)
    xcen, ycen = sizes/2, sizes/2
    ras, decs  = [], []
    xpos, ypos = [], []
    
    print(f' - sizes: {sizes}')
    
    try:
        for i in range(cf):
            ras.append(t['RA_'+str(i+1)])
            decs.append(t['DEC_'+str(i+1)])
    except ValueError as ve000:
        print(' (!) 0: ', ve000)
        pass
    
    try:
        for i in range(cf):
            xpos.append(t['Xposn_'+str(i+1)])
            ypos.append(t['Yposn_'+str(i+1)])
    except ValueError as ve001:
        print(' (!) 1: ', ve001)
        pass
    
    try:
        ras.append(t['_RAJ2000'])
        decs.append(t['_DEJ2000'])
        xpos.append(t['Xpos'])
        ypos.append(t['Ypos'])
    except ValueError as ve002:
        print(' (!) 1: ', ve002)
        pass
    
    print (f'len(ras): {len(ras)}')
    # approx, is valid only at equator
    # theta is PA of src2 from src1
    # phi is PA of src from image centre
    # no rotation => random; rotation => linear
    
    #ra1, dec1 = ras[-1], decs[-1]
    #ra1, dec1 = ras[0], decs[0]
    print("\n\n *** \nWE MAY NEED A PROPER MeerKAT REFERENCE image for the rotation calculations\n ***\n\n")
    l1 = list(arange(0,cf,1))
    pl.figure(figsize=(10,15))
    nfigrows = N.round(cf/2.)+1
    try:
        for i in range(cf):
            li = list(zip(l1,[i]))[0]
            i0 = cref #li[0] #make sure that this is the image with the smallest beam (highest ref freq image)
            ii = li[1]
            print (f' > i0 x ii: {i0} x {ii}')
            ra1, dec1 = ras[ i0 ], decs[ i0 ]
            ra2, dec2 = ras[ii], decs[ii]
            dx, dy    = xpos[ii]-xcen[ii], ypos[i]-ycen[ii]
            theta = N.arctan2( dec2-dec1, (ra1-ra2)*cos(N.mean(dec1)/180*pi) )*180/pi
            phi   = N.arctan2(dy, dx)*180/pi
            phi   = N.where(phi<-90,phi+360,phi)
            dist  = N.sqrt( dx**2 + dy**2 )
            
            pl.subplot(nfigrows,2,ii+1)
            pl.scatter(phi, theta, s=dist*1.0/300, alpha=0.4)
            pl.title('obs {} x {}: {} x {}'.format(i0,ii, bases[i0],bases[ii]),fontsize=8)
            if ii in [1,3]: pl.ylabel('Phi (deg)')
            if ii in [3,4]: pl.xlabel('Theta (deg)')
    except IndexError as ie002:
        print(' (!) 2: ', ie002)
        pass

    pl.suptitle('PA of src from centre (phi) vs angle of offset (theta)')
    pl.tight_layout()
    pl.subplots_adjust(left=0.076, right=0.982, top=0.937, bottom=0.041,wspace=0.186,hspace=0.223)
    pl.savefig(outdir+'/%s-corrtest_rotation_phi_theta.png'%base)
    pl.show()
    
    return ()

def quick_stats(csv_file, base=None, do_rms=1, do_mean=0, do_max=1, do_min=1, inmark=None):
    """
    Plot RMS SEDs per mode.
    inputs:
    csv_file - str(), path to csv file
    base - str(), human readable reference to the data 
    do_rms - bolean, enable plot of rms SEDs
    do_mean - bolean, enable plot of mean SEDs
    do_max - bolean, enable plot of max SEDs 
    do_min - bolean, enable plot of min SEDs
    inmark - str(), custom plot marker
    """
    import pandas as pd
    import random as r
    
    m   = ['o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X']
    if inmark == None:
        inmark = r.choice(m)
    if base == None:
        base = ''
        
    df  = pd.read_csv(csv_file)
    frq = array(list(df['freq0']))
    rms = array(list(df['fd_rms']))
    avg = array(list(df['fd_avg']))
    Fmx = array(list(df['fd_max']))
    Fmn = array(list(df['fd_min']))
    mode= list(df['#mode'])[0]
    
    rms_avg = nanmean(rms)
    fd_mean = nanmean(avg)
    fd_max  = nanmean(Fmx)
    fd_min  = nanmean(Fmn)
    
    std_rms = nanstd(rms)
    std_avg = nanstd(avg)
    std_max = nanstd(Fmx)
    std_min = nanstd(Fmn)
    
    if do_mean:
        plt.subplot(111)
        plt.title('Stokes I mean for each mode')
        plt.scatter(frq, avg, marker=inmark, label='{}-{}, avg: {:.2e}$\pm{:.2e}$ Jy/bm'.format(mode,base, fd_mean, std_avg))
        plt.legend(loc=0)
        plt.ylabel('Mean [Jy/bm]', fontsize=20)
        plt.xlabel('Frequency [Hz]', fontsize=20)
        plt.grid(True, ls = ':')
        #plt.tight_layout()
    
    
    if do_rms:
        plt.subplot(311)
        plt.title('Stokes I RMS for each mode')
        plt.scatter(frq, rms, marker=inmark, label='{}-{}, avg: {:.2e}$\pm{:.2e}$ Jy/bm'.format(mode,base, rms_avg, std_rms))
        plt.legend(loc=0)
        plt.ylabel('rms [Jy/bm]', fontsize=20)
        if not do_max:
            plt.xlabel('Frequency [Hz]', fontsize=20)
        plt.grid(True, ls = ':')
        #plt.tight_layout()
        print('{}, {:.2e}, {:.2e}, {:.2e}'.format(mode, fd_min, fd_max, rms_avg)) 
    
    if do_max:
        plt.subplot(313)
        plt.title('Stokes I max for each mode')
        plt.scatter(frq, Fmx, marker=inmark, label='{}-{}, avg: {:.2e}$\pm{:.2e}$ Jy/bm'.format(mode,base, fd_max, std_max))
        plt.legend(loc=0)
        plt.ylabel('Max [Jy/bm]', fontsize=20)
        plt.xlabel('Frequency [Hz]', fontsize=20)
        plt.grid(True, ls = ':')
        #plt.tight_layout()
    
    if do_min:
        plt.subplot(312)
        plt.title('Stokes I min for each mode')
        plt.scatter(frq, Fmn, marker=inmark, label='{}-{}, avg: {:.2e}$\pm{:.2e}$ Jy/bm'.format(mode,base, fd_min, std_min))
        plt.legend(loc=0)
        plt.ylabel('Min [Jy/bm]', fontsize=20)
        if not do_max:
            plt.xlabel('Frequency [Hz]', fontsize=20)
        plt.grid(True, ls = ':')
        #plt.tight_layout()
    
    return(fd_mean, fd_max, rms_avg)

def RA_Dec_offsets(RA1, Dec1, RA2, Dec2, plotit=None):
    """
    # assumes coordinates from the same reference frame
    # simple approximation of angular separation
    # pyasl.getAngDist() seems to be wrong in some instances!
    # assumes inputs in deg and so converts Dec1 into rad in cos()
    output:
        - delta(RA) and delta(Dec) in deg
        - delta(RA) vs delta(Dec) plot if plotit = True
    """
    
    if hasattr(RA1, '__len__'):
        RA1  = array(RA1)
        Dec1 = array(Dec1)
        RA2  = array(RA2)
        Dec2 = array(Dec2)
    dRA  = (RA2-RA1)*cos(radians(Dec1))
    dDec = Dec2-Dec1
    
    if plotit and hasattr(dRA, '__len__'):
        dRA_avg  = nanmean(dRA)
        dDec_avg = nanmean(dDec)
        dRA_std  = nanstd(dRA)
        dDec_std = nanstd(dDec)
        
        plt.scatter(dRA, dDec, marker='+', alpha=0.5,
                    label=base+': avg RA = {:.2e} $\pm$ {:.2e}\n   avg Dec = {:.2e} $\pm$ {:.2e}'.format(dRA_avg, dRA_std, dDec_avg,dDec_std) ) 
        plt.legend(loc=0)  
        plt.ylabel('dDEC [deg]', fontsize=15)  
        plt.xlabel('dRA [deg]', fontsize=15)  
        plt.grid(True, ls = ':') 
        plt.tight_layout()
        plt.savefig('offsets-RA-Dec.png')
        plt.show()
    
    return(dRA, dDec)

def Get_posnRMS(self, min_axs, maj_axs, SNR):
    beamsize = sqrt( min_axs**2 + maj_axs**2 )
    pos_RMS  = beamsize/(2.*SNR)
    
    return (pos_RMS)

def corr_test_offsets(tLxL0=None, tUxU0=None):
    """
    Plot MeerKAT L-band to L-band RA vs Dec offsets.
    Tested with Offsets between two meerKAT observations of the Jul and May 2020 
    correlator imaging test observations. 
    
    input:
    - tLxL0 : L to L band xmatch full file path
    - tUxU0 : U to U band xmatch full file path
    
    """
    
    if tLxL0 == None or tUxU0 == None:
        catdir = "/home/samuel/SARAO/imagin_varification/2020-correlator-test/July/"
        tLxL   = atpy.Table(catdir+'L4kJul-x-L4kMay.tbl')
        tUxU   = atpy.Table(catdir+'U4kJul-x-U4kMay.tbl')  
    
    else:
        catdir = os.path.dirname(tLxL0)
        tLxL   = atpy.Table(tLxL0)
        tUxU   = atpy.Table(tUxU0)  
    print('\n\ndRA_avg, dRA_std, dDec_avg, dDec_std, avg_sep, std_sep .. [units in asec]')
    plt.figure()
    plt.title('RA vs Dec offsets', fontsize=15) 
    for tb in ([tLxL, 'LxL'],[tUxU, 'UxU']):
        t = tb[0]; base = tb[1] 
        dRA_i,dDec_i = RA_Dec_offsets( t['RA_1'], t['DEC_1'], t['RA_2'], t['DEC_2'] ) 
        dRA_i,dDec_i = dRA_i*3600, dDec_i*3600
        
        dRA_avg  = nanmean(dRA_i)
        dDec_avg = nanmean(dDec_i)
        dRA_std  = nanstd(dRA_i)
        dDec_std = nanstd(dDec_i)
        sep_i    = t['Separation']
        #a, b     = t['']
        avg_sep  = nanmean(sep_i)
        std_sep  = nanstd(sep_i)
        
        print(f'{dRA_avg}, {dRA_std}, {dDec_avg}, {dDec_std}, {avg_sep}, {std_sep}')
        
        #plt.subplot(311)
        plt.scatter(dRA_i, dDec_i, marker='+', alpha=0.5,
                    label=base+': avg RA = {:.2e} $\pm$ {:.2e}\"\n   avg Dec = {:.2e} $\pm$ {:.2e}\"'.format(dRA_avg, dRA_std, dDec_avg,dDec_std) ) 
        plt.legend(loc=0)  
        plt.ylabel('dDEC [arcsec]', fontsize=15)  
        plt.xlabel('dRA [arcsec]', fontsize=15)  
        plt.grid(True, ls = ':')  
        
        ##plt.subplot(312)
        #plt.scatter(sep_i,dRA_i, marker='+', alpha=0.5,label=base+' avg sep: {:.2e} $\pm$ {:.2e}\"'.format(avg_sep,std_sep) ) 
        #plt.legend(loc=0)  
        ##plt.xlabel('xmatch speration [arcsec]', fontsize=15)  
        #plt.ylabel('dRA [arcsec]', fontsize=15)  
        #plt.grid(True, ls = ':') 
        
        ##plt.subplot(313)
        #plt.scatter(sep_i,dDec_i,marker='+', alpha=0.5) 
        #plt.legend(loc=0)  
        #plt.xlabel('xmatch speration [arcsec]', fontsize=15)  
        #plt.ylabel('dDec [arcsec]', fontsize=15)  
        #plt.grid(True, ls = ':') 
    plt.tight_layout()
    plt.savefig(catdir+'/offsets-RA-Dec.png')
    plt.show()
    print('\n\n')
    return ()


########################################################################
if __name__ == '__main__':
    process_start  = time.time()
    timestr        = time.strftime("%Y%m%d-%H%M%S")
    do_subbands    = args.grab_subbandIms
    in_fitscube    = args.in_fitscube
    do_SourceFind  = args.do_SourceFinding
    src_findIm     = args.src_findIm
    do_Xmatch      = args.do_Xmatch
    check_rotation = args.check_rotation
    pwd            = args.cwd #"/scratch2/slegodi/raw_obs/im_verify/" #os.getcwd()
    Xfile          = args.xfile
    test_descript  = args.test_descript
    
    mymarkers      = ['o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X']
    
    print('\n > START TIME        : ',timestr)
    print(' > PWD               : ',pwd)
    print(' > cube              : ',in_fitscube)
    print(' > Source find image : ',src_findIm)
    print(' > Do Source finding : ',do_SourceFind)
    print(' > Do subband image  : ',do_subbands)
    print(' > Do Xmatching      : ',do_Xmatch)
    print(' > Check im rotation : ',check_rotation)
    
    print('\n\n')
    if do_subbands:
        generic_get_subbandImage(in_fitscube)
        
    if do_SourceFind:
        if not src_findIm:
            src_findIm = str(input('Enter full file path to source finding image: '))
        find_sources(src_findIm)
        
    if check_rotation:
        print(f' > X-match files parent directory: {os.path.dirname(Xfile)}')
        iref = int(input('int(), Reference catalog subscript in x-match file "%s": '%os.path.basename(Xfile)))
        nXcat= int(input('int(), total number of catalogs x-matched in file "%s" : '%os.path.basename(Xfile)))
        im_rot(Xfile, test_descript, 1, nXcat, iref)
        corr_test_offsets(None, None)


"""
To do:
    - get total flux of central source from source finding catalogs
    - image rms via CASA
    - get theoretical rms (based on actual data size - time, bw etc)
    - ratio of the two
    - OR else scale image rms by sqrt of time on src and num baselines
    
    - compare May to April and Jul to April position offsets instead of comp to SUMMS
    - plot for l and u
    - plot may-apr and jul-apr offsets, apr L-band image is the reference 
    - tabulate mean and rms of offset in ra and dec

"""
