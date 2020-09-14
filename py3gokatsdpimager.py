"""
- Script (python3) to perform imaging via the katsdp imager: katsdpcontim and also do analysis on produced images
- import friendly
- requires commnad line arguments so use option -h to see options
- launches docker with name given my parameter: descript+timestamp
- to find docker logs do the following on the command line:
 $ docker logs -f <descript+timestamp>
- to list all docker sessions, do:
 $ docker ps
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
parser.add_argument('-screen_on', action='store_true', help='True/False: Run the entire script in a detacthed screen')
parser.add_argument("-in_rdb", required=False, help='Required! - RDB file tag, e.g. rdb1. See script for a list of RDBs for this') 
parser.add_argument("-pol", required=False, default='HH, VV', help="Not required - Polarisation to image ( H or V or HH, VV )")
parser.add_argument("-toff", required=False, default=0, help="Float. Not required - Offset to add to all timestamps, in seconds")
parser.add_argument("-maxInt_m", required=False,default=0.001, help="Float. Not required - Maximum integration (min)")
parser.add_argument("-solPInt", required=False, default=0.1333, help="Float. Not required - phase SC Solution interval (min)")
parser.add_argument('-im_FOV',  required=False, default=3.0,    help='Float, image field of view in degrees.')
parser.add_argument("-stokes", type=str,required=False, default='I', help="String. Not required - Stokes to process")

parser.add_argument('-do_image', action='store_true', help='True/False: Turn on the katsdpcontim imager')
parser.add_argument('-image_targ', action='store_true', help='True/False: Image/analyse only the science target')
parser.add_argument('-image_gcal', action='store_true', help='True/False: Image/analyse only the gain calibrators')
parser.add_argument('-image_bpcal', action='store_true', help='True/False: Image/analyse only the bandpass/primary calibrators')

parser.add_argument('-get_planes', action='store_true',\
    help='True/False: Split the lowband, midband, and highband planes from the cuboid image produced by katsdpcontim imager')
parser.add_argument('-do_srcfind', action='store_true', help='True/False: Turn source finding function')

parser.add_argument('-do_low', action='store_true', help='True/False: Do imaging and/or source finding on the lowband image')
parser.add_argument('-do_mid', action='store_true', help='True/False: Do imaging and/or source finding on the midband image')
parser.add_argument('-do_high', action='store_true', help='True/False: Do imaging and/or source finding on the highband image')
parser.add_argument('-do_wide', action='store_true', help='True/False: Do imaging and/or source finding on the wideband image')
parser.add_argument('-do_fullband', action='store_true', help='True/False: Do imaging and/or source finding on the full band image')
parser.add_argument('-selct_band', action='store_true', help='True/False: Do imaging and/or source finding on the selected band image. Select according to fitscuboid structure.')

parser.add_argument('-select_chans', action='store_true',\
    help='True/False: select data according to scan numbers (see script for details - default will choose according to frequencies)')
parser.add_argument('-select_freqs', action='store_true', help='True/False: select data according to frequency range (see script for details)')

parser.add_argument('-per_scan_images', action='store_true',help='True/False: image each scan separately')
parser.add_argument('-Allscan_image', action='store_true',help='True/False: image all scans into one image')

parser.add_argument('-use_smoothed_image', action='store_true',\
    help='True/False: select a smoothed fits image. Currently set to images with suffix smooth2sumss.fits')
parser.add_argument("-CWD", required=False, default=os.getcwd(), help="Path. Not required. Path to working directory")


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
        isl_thresh     = 5.0                # SNR limit for source finding island region
        pix_thresh     = 1.5*isl_thresh      # SNR limit for source detection/identification
        fluxthreshold  = 1.0e-06
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

def separation(ra1,dec1,ra2,dec2):
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

def smoothimage(image_in, bmaj, bmin, bpa, base):
    # this runs in a casa session 
    # bmaj, bmin, bpa in the form: Xunit
    # base = human readable reference/name string
    outname = image_in[:-5]+'.%s.image'%base
    outfits = outname[:-6]+'.fits'
    imsmooth(imagename=image_in,kernel="g",major=bmaj,minor=bmin,pa=bpa,targetres=True,kimage="",scale=-1.0,region="",box="",\
        chans="",stokes="",mask="",outfile=outname,stretch=False,overwrite=True,beam="")
    
    exportfits(imagename=outname,fitsimage=outfits,velocity=False,optical=False,bitpix=-32,minpix=0,maxpix=-1,\
        overwrite=True,dropstokes=False,stokeslast=True,history=True,dropdeg=False)
    
    print('> Input image        : %s'%image_in)
    print('> Smoothed casa image: %s'%outname)
    print('> Smoothed fits image: %s'%outfits)
    print('> Smoothed beam      : bmaj: {}, bmin: {}, bpa: {}'.format(bmaj, bmin, bpa))
    
    return ()

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
    
    low_fileout = fitscuboid[:-5]+'.lowband.fits'
    hdu_low.writeto(low_fileout, overwrite=True)
    
    print('\nlow band file   : ', low_fileout)
    
    mid_fileout = fitscuboid[:-5]+'.midband.fits'
    hdu_mid.writeto(mid_fileout, overwrite=True)
    print('mid band file   : ', mid_fileout)
    
    hgh_fileout = fitscuboid[:-5]+'.hghband.fits'
    hdu_hgh.writeto(hgh_fileout, overwrite=True)
    print('high band file  : ', hgh_fileout)
    
    wide_fileout = fitscuboid[:-5]+'.wideband.fits'
    hdu_wide.writeto(wide_fileout, overwrite=True)
    print('wide band file  : ', wide_fileout, '\n\n')
    
    slct_fileout = fitscuboid[:-5]+'.freq%.3fGHz.fits' %(hdr_frq/(1.e9))
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
        algo = 'fit-half'
        
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
        mode = str(raw_input("CBID is %s. Enter Observation mode (string) : "%base))
    
    show_figs      = True
    timestr        = time.strftime("%Y%m%d-%H%M%S")
    outfile        = open(fitsfile[:-4]+mode+'-'+algo+'-channelised_imstats.txt','w+')
    outfile2       = open(fitsfile[:-4]+mode+'-'+algo+'-channelised_imstats.csv','w')
    
    dd, hdr = fits.getdata(fitsfile, header=True)
    nplanes = hdr['NAXIS3']
    nx      = hdr['NAXIS1']
    ny      = hdr['NAXIS2']
    N0      = hdr['NAXIS']
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
            si      = imstat(imagename=p_out, algorithm=algo, center='zero')
            #si      = imstat(imagename=p_out)
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


def generic_get_subbandImage(fitscuboid):
    """
    - extracts fits planes from fits cuboid image
    - fitscuboid from katsdpimager: "katsdpcontim"
    - returns 3 images at low, mid and high part of the full band of the input fits
    """
    
    freqs_dct = {}; freqs_ar = []
    dd, hdr = fits.getdata(fitscuboid, header=True)
    nx      = hdr['NAXIS1']
    ny      = hdr['NAXIS2']
    nz      = hdr['NAXIS3']
    Nn      = hdr['NAXIS']
    
    low_bim  = dd[0,1]      # low band image at freqs ~ [8.989462890625E+08; 9.658212890625E+08]
    mid_bim  = dd[0,5]      # mid band image at freqs ~ [1.279088867188E+09; 1.367907226562E+09]
    hgh_bim  = dd[0,7]      # high band image at freqs ~ [1.465084960938E+09; 1.567278320312E+09]
    avg_wim  = dd[0,0]      # full band image
    slct_im  = dd[0,nz-1]
    
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
    
    low_fileout = fitscuboid[:-5]+'.lowband.fits'
    hdu_low.writeto(low_fileout, overwrite=True)
    
    print('\nlow band file   : ', low_fileout)
    
    mid_fileout = fitscuboid[:-5]+'.midband.fits'
    hdu_mid.writeto(mid_fileout, overwrite=True)
    print('mid band file   : ', mid_fileout)
    
    hgh_fileout = fitscuboid[:-5]+'.hghband.fits'
    hdu_hgh.writeto(hgh_fileout, overwrite=True)
    print('high band file  : ', hgh_fileout)
    
    wide_fileout = fitscuboid[:-5]+'.wideband.fits'
    hdu_wide.writeto(wide_fileout, overwrite=True)
    print('wide band file  : ', wide_fileout, '\n\n')
    
    slct_fileout = fitscuboid[:-5]+'.freq%.3fGHz.fits' %(hdr_frq/(1.e9))
    hdu_slct.writeto(slct_fileout, overwrite=True)
    print(' > selected band centre: %.3f GHz'%(hdr_frq/(1.e9)))
    print(' > selcted band file   : ', slct_fileout, '\n\n')
    
    fitsdir = os.path.dirname(fitscuboid)
    print(' > Output fits in: ', fitsdir)
    return ()

def casa_regrid(in_image, im_template):
    # (!) to be ran in casa session
    td_temp = imregrid(imagename=im_template, template="get")
    imregrid(imagename=in_image, output=in_image[:-5]+'.regrid.fits', template=td_temp, overwrite=True)
    return()

def RADEC_HHMMSS(ra_deg, dec_deg):
    ra_deg = float(ra_deg); dec_deg = float(dec_deg)
    posn   = sk(ra=ra_deg,dec=dec_deg,unit='deg',frame='icrs')
    ra_hms = posn.ra.hms; dec_dms = posn.dec.dms
    
    RA_HHMMSS  = '%02d:%02d:%07.4f'%(ra_hms[0],ra_hms[1],ra_hms[2])
    DEC_DDMMSS = '%02d:%02d:%07.4f'%(dec_dms[0],abs(dec_dms[1]),abs(dec_dms[2]))
    radec_str  = '{}, {}'.format(RA_HHMMSS, DEC_DDMMSS)
    return(radec_str, RA_HHMMSS, DEC_DDMMSS)

def MeerKAT_rms(Ants, Df, T0, SEFD, S_peak):
    # Stokes I rms noise in Jy on MeerKAT (64, 13m Antennas)
    # from: https://drive.google.com/file/d/16LXJQ8KXUqWx2lo7GFm-fSIUfVL41w1k/view
    if Df == None:
        Df = 856.e6     #MeerKAT full bandwidth in Hz
    if SEFD == None:
        SEFD  = 443.5
    
    NN_1 = Ants*(Ants-1.)
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

def renamer(dir0, pref, suff, new_suff):
    # dir0 is full puth to the parent directory
    import os, sys
    for pth_i in os.listdir(dir0):
        if pth_i.startswith(str(pref)) and pth_i.endswith(str(suff)):
            pth_i  = dir0+pth_i
            pth_ii = pth_i+str(new_suff)
            os.system('mv %s %s'%(pth_i, pth_ii))
    return ()
    
########################################################################################################################
if __name__ == '__main__':
    process_start  = time.time()
    timestr        = time.strftime("%Y%m%d-%H%M%S")
    pwd            = args.CWD#"/scratch2/slegodi/raw_obs/im_verify/" #os.getcwd()
    pwd0           = "/scratch2/slegodi/raw_obs/im_verify/"
    print('\n > PWD: ',pwd)
    ######################################################

    # the selected bandwidth needs to be divisible by 8
    full_band  = slice(0,4096); full_fs   = (0.856e+9,1.711791e+9)   # Hz
    wide_band  = slice(204,3892);  all_fs = (0.8880791015625e+9, 1.711477539062e+9)   # Hz
    low_band   = slice(570,770);   low_fs = (0.97512109375e+9, 1.01691796875e+9)      # Hz
    mid_band   = slice(2250,2450); mid_fs = (1.32621484375e+9, 1.36801171875e+9)      # Hz
    hgh_band   = slice(2900,3100); hgh_fs = (1.4620546875e+9, 1.5038515625e+9)        # Hz
    bandwdth0  = 0.0417968e+9 # in Hz: (208.984 kHz)*200 channels
    fov        = args.im_FOV          # Radius of field to image (deg)
    maxInt_time= args.maxInt_m
    solPInt    = args.solPInt
    do_all_ObsIDs = 0 # (!) this should almost always be False/0

    # RDB (full) links, numelrically named according to obs date
    target = 'J2147-8132'
    bpcal  = 'J1939-6342'
    gcal   = 'J2147-7536'

    ### L-band obs ###
    rdb0 = {'rdb': "http://archive-gw-1.kat.ac.za:7480/1563312811/1563312811_sdp_l0.full.rdb",\
        'tgt_scans':[6,12],'bcal_scans':[1,14],'gcal_scans':[3,9], 'dumps':'2s'}
    rdb1 = {'rdb': "http://archive-gw-1.kat.ac.za:7480/1564254957/1564254957_sdp_l0.full.rdb",\
        'tgt_scans':[5,9],'bcal_scans':[1,11],'gcal_scans':[3,7], 'dumps':'8s'}
    rdb2 = {'rdb': "http://archive-gw-1.kat.ac.za:7480/1564259663/1564259663_sdp_l0.full.rdb",\
        'tgt_scans':[5,9],'bcal_scans':[1,11],'gcal_scans':[3,7], 'dumps':'8s'}
    rdb3 = {'rdb': "http://archive-gw-1.kat.ac.za:7480/1564263655/1564263655_sdp_l0.full.rdb",\
        'tgt_scans':[6,10],'bcal_scans':[1,11],'gcal_scans':[4,8], 'dumps':'2s'}

    # New data
    rdb4 = {'rdb': "http://archive-gw-1.kat.ac.za:7480/1565560531/1565560531_sdp_l0.full.rdb",\
        'tgt_scans':[5,9], 'bcal_scans':[1], 'gcal_scans':[3,7], 'dumps':'8s'}
    rdb5 = {'rdb': "http://archive-gw-1.kat.ac.za:7480/1565564738/1565564738_sdp_l0.full.rdb",\
        'tgt_scans':[5,9], 'bcal_scans':[1,11], 'gcal_scans':[3,7], 'dumps':'8s'}
    rdb6 = {'rdb': "http://archive-gw-1.kat.ac.za:7480/1565646959/1565646959_sdp_l0.full.rdb",\
        'tgt_scans':[4,8], 'bcal_scans':[0,10], 'gcal_scans':[2,6], 'dumps':'8s'}
    rdb7 = {'rdb': "http://archive-gw-1.kat.ac.za:7480/1565651658/1565651658_sdp_l0.full.rdb",\
        'tgt_scans':[5,9], 'bcal_scans':[1,11], 'gcal_scans':[3,7], 'dumps':'8s'}
    rdb8 = {'rdb': "http://archive-gw-1.kat.ac.za:7480/1566764293/1566764293_sdp_l0.full.rdb",\
        'tgt_scans':[1], 'bcal_scans':arange(start=13, stop=42, step=1, dtype=int),'gcal_scans':[1], 'dumps':'8s'}

    rdb11= {'rdb': "http://archive-gw-1.kat.ac.za:7480/1547977303/1547977303_sdp_l0.full.rdb",\
        'tgt_scans':arange(5,33,4), 'bcal_scans':[1,47,93],'gcal_scans':arange(3,155,4), 'dumps':'8s'}

    ### UHF obs ###
    #rdb10 = {'rdb': "http://archive-gw-1.kat.ac.za:7480/1551961010/1551961010_sdp_l0.full.rdb",\
        #'tgt_scans':arange(5,148,4), 'bcal_scans':[1],'gcal_scans':arange(3,143,4), 'dumps':'8s'}

    rdb9 = {'rdb': "http://archive-gw-1.kat.ac.za:7480/1551961010/1551961010_sdp_l0.full.rdb",
            'tgt_scans':arange(5,21,4), 'bcal_scans':[1],'gcal_scans':arange(3,143,4), 'dumps':'8s'}

    rdb10= {'rdb': "http://archive-gw-1.kat.ac.za:7480/1563925636/1563925636_sdp_l0.full.rdb",
            'tgt_scans':[37,39,61,91], 'bcal_scans':[7,3,1],'gcal_scans':[9,11,13,35], 'dumps':'8s'}

    rdb12= {'rdb': "http://archive-gw-1.kat.ac.za:7480/1567882150/1567882150_sdp_l0.full.rdb",
            'tgt_scans':concatenate([arange(5,46,4),arange(59,100,4), arange(115,156,4), arange(169,210,4), arange(218,224,4)]),
            'bcal_scans':[1,49,51,53,55,103,109,159,163,213,214],'gcal_scans':arange(3,220,4), 'dumps':'8s'}

    rdb13= {'rdb': "http://archive-gw-1.kat.ac.za:7480/1569521762/1569521762_sdp_l0.full.rdb",
            'tgt_scans':list(arange(5,30,4))+list(arange(37,62,4))+list(arange(69,94,4))+list(arange(101,126,4)),
            'bcal_scans':[1],'gcal_scans':arange(3,128,4), 'dumps':'8s'}

    rdb14={'rdb': "http://archive-gw-1.kat.ac.za:7480/1571007678/1571007678_sdp_l0.full.rdb",
        'tgt_scans':list(arange(5,50,4)),'bcal_scans':[1,29,53],'gcal_scans':list(arange(3,52,4)),
        'dumps':'4s'}

    #rdb15= {'rdb': "http://archive-gw-1.kat.ac.za:7480/1569087057/1569087057_sdp_l0.full.rdb",
            #'tgt_scans':concatenate([arange(5,46,4),arange(59,100,4), arange(115,156,4), arange(169,210,4), arange(218,224,4)]),
            #'bcal_scans':[1,51,55,105,49,53,103,109,159,163,213,214],'gcal_scans':arange(3,220,4), 'dumps':'8s'}

    rdb15= {'rdb': "http://archive-gw-1.kat.ac.za:7480/1569087057/1569087057_sdp_l0.full.rdb",
            'tgt_scans':[5],
            'bcal_scans':[1,51,55,105,49,53,103,109,159,163,213,214],'gcal_scans':arange(3,220,4), 'dumps':'8s'}

    rdb16= {'rdb': "/scratch2/slegodi/raw_obs/im_verify/1571058370_sdp_l0.fixed.full.rdb",
        'tgt_scans':[],'bcal_scans':[],'gcal_scans':[]}

    rdb17= {'rdb':"http://archive-gw-1.kat.ac.za:7480/1568825523/1568825523_sdp_l0.full.rdb",
            'tgt_scans':[6,10,13,15,18,20],'bcal_scans':[1,37,73,107],'gcal_scans':[2,4,8,12,17,22]}

    rdb18= {'rdb': "https://archive-gw-1.kat.ac.za/1573490035/1573490035_sdp_l0.full.rdb?token=eyJ0eXAiOiJKV1QiLCJhbGciOiJFUzI1NiJ9.eyJzdWIiOiJkZXYiLCJwcmVmaXgiOlsiMTU3MzQ5MDAzNSJdLCJzY29wZXMiOlsicmVhZCJdLCJpc3MiOiJrYXQtYXJjaGl2ZS5rYXQuYWMuemEiLCJhdWQiOiJhcmNoaXZlLWd3LTEua2F0LmFjLnphIiwiaWF0IjoxNTc0NDE3MTAwLCJleHAiOjE1NzUwMjE5MDB9.CfNElL8DIx7ipHIQWADig8uSqzcKTwZSolMkbFyU_o-cwZ4r2m30aOAk_wPSPsDXHcAWcNHHrfKxW-R-y8uO-w",'tgt_scans':[1],'bcal_scans':[1],'gcal_scans':[1]}
    
    #rdb19= {'rdb':"/scratch2/slegodi/raw_obs/ScienceVerification/1562441876_sdp_l0.full.rdb",
            #'tgt_scans':[5,9,13,17,63,67,71,119,125,129],'bcal_scans':[1,43,85,131],'gcal_scans':[3,7,11,61,65,69,123,127,133]}
            
    rdb19= {'rdb':"https://archive-gw-1.kat.ac.za/1562441876/1562441876_sdp_l0.full.rdb?token=eyJ0eXAiOiJKV1QiLCJhbGciOiJFUzI1NiJ9.eyJpc3MiOiJrYXQtYXJjaGl2ZS5rYXQuYWMuemEiLCJhdWQiOiJhcmNoaXZlLWd3LTEua2F0LmFjLnphIiwiaWF0IjoxNTgxNDA4NDExLCJwcmVmaXgiOlsiMTU2MjQ0MTg3NiJdLCJleHAiOjE1ODIwMTMyMTEsInN1YiI6InNsZWdvZGlAc2thLmFjLnphIiwic2NvcGVzIjpbInJlYWQiXX0.L7wnQZ_8Ku3zwjU90FRuzb0zzpTRDz2vk6vfJ4RaXGhyI_Cuq0apJ2TyivjbqJrXEb950vj74OdF9eA9ZlPpdw",
            'tgt_scans':[5,9,13,17,63,67,71,119,125,129],'bcal_scans':[1,43,85,131],'gcal_scans':[3,7,11,61,65,69,123,127,133]}
    
    rdb20= {'rdb': 'https://archive-gw-1.kat.ac.za/1586074843/1586074843_sdp_l0.full.rdb?token=eyJ0eXAiOiJKV1QiLCJhbGciOiJFUzI1NiJ9.eyJpc3MiOiJrYXQtYXJjaGl2ZS5rYXQuYWMuemEiLCJhdWQiOiJhcmNoaXZlLWd3LTEua2F0LmFjLnphIiwiaWF0IjoxNTg2MjQ3NDAwLCJwcmVmaXgiOlsiMTU4NjA3NDg0MyJdLCJleHAiOjE1ODY4NTIyMDAsInN1YiI6InNsZWdvZGlAc2thLmFjLnphIiwic2NvcGVzIjpbInJlYWQiXX0.Un-UhbhvZh1Ra3avJYKnSlJLehYki1pFfqHFylWGQvd0TzIRty5ENci04-5CT5XcHPxcRwPYjouGxqGNF8tgeg',
            'tgt_scans' : [5,9],
            'bcal_scans': [1,11],
            'gcal_scans': [3,7]}
    
    rdb21= {'rdb': 'https://archive-gw-1.kat.ac.za/1586305859/1586305859_sdp_l0.full.rdb?token=eyJ0eXAiOiJKV1QiLCJhbGciOiJFUzI1NiJ9.eyJpc3MiOiJrYXQtYXJjaGl2ZS5rYXQuYWMuemEiLCJhdWQiOiJhcmNoaXZlLWd3LTEua2F0LmFjLnphIiwiaWF0IjoxNTg3NzE5OTU1LCJwcmVmaXgiOlsiMTU4NjMwNTg1OSJdLCJleHAiOjE1ODgzMjQ3NTUsInN1YiI6InNsZWdvZGlAc2thLmFjLnphIiwic2NvcGVzIjpbInJlYWQiXX0.lcpfdxwYDj_XUKnKfmcmyyziXf_tC_a8HKffKAY69y77695GGBP86Nq_VwPsAJ8w0TiydFyvrAa0bfBrD2xhLQ',
            'tgt_scans' : [5,9],
            'bcal_scans': [1,11],
            'gcal_scans': [3,7]}
    
    rdb22={'rdb': '/scratch2/slegodi/raw_obs/correlator_test/rdbs/1586305859_sdp_l0.full.rdb',
           'tgt_scans' : [5,9],
            'bcal_scans': [1,11],
            'gcal_scans': [3,7]}
    
    rdb23={'rdb': 'https://archive-gw-1.kat.ac.za/1596163435/1596163435_sdp_l0.full.rdb?token=eyJ0eXAiOiJKV1QiLCJhbGciOiJFUzI1NiJ9.eyJpc3MiOiJrYXQtYXJjaGl2ZS5rYXQuYWMuemEiLCJhdWQiOiJhcmNoaXZlLWd3LTEua2F0LmFjLnphIiwiaWF0IjoxNTk2NTQ2Mjk5LCJwcmVmaXgiOlsiMTU5NjE2MzQzNSJdLCJleHAiOjE1OTcxNTEwOTksInN1YiI6InNsZWdvZGlAc2thLmFjLnphIiwic2NvcGVzIjpbInJlYWQiXX0.Vedw9sq3x5angxvDJF1bRjU5WZhUjKsvIYgrucjlaUo1HRNBKod7J9C69jH79KDbdm65Ici4fnbKfSTszNq1-g',
           'tgt_scans' : [4,8],
            'bcal_scans': [1,11],
            'gcal_scans': [3,7]}

    RDBs = {'rdb0': [rdb0, '',''],
            'rdb1': [rdb1,'rdb1',''],
            'rdb2': [rdb2,'rdb2',''],
            'rdb3': [rdb3,'rdb3',''],
            'rdb4': [rdb4,'CMC2-1K',''],
            'rdb5': [rdb5,'CMC1-1K',''],
            'rdb6': [rdb6,'CMC2-4K',''],
            'rdb7': [rdb7,'CMC1-4K',''],
            'rdb8': [rdb8,'CMC2-4K',''],
            'rdb9': [rdb9, '4K-UHF-DEEP2','1551961010'],
            'rdb10':[rdb10, '4K-UHF-CalImaging','1563925636'],
            'rdb11':[rdb11, '4K-L-DEEP2','1547977303'],
            'rdb12':[rdb12, 'obs1-1567882150','1567882150'],
            'rdb13':[rdb13,'S190814bv-50pc-Initial-Sky-Map','1569521762'],
            'rdb14':[rdb14,'MV-01_Epoch_2','1571007678'],
            'rdb15':[rdb15,'obs2-1569087057','1569087057'],
            'rdb16':[rdb16,'narrowband_imaging_test_on_G330.89-0.36','1571058370'],
            'rdb17':[rdb17,'Niruj_UHF-ionosphere','1568825523'],
            'rdb18':[rdb18,'Niruj-J0408-6545.32k.Stability.Obs','1573490035'],
            'rdb19':[rdb19, 'SV-data-FRB180924', '1562441876'],
            'rdb20':[rdb20, 'corrtest-L4k', '1586074843'],
            'rdb21':[rdb21, 'corrtest-L4kv2', '1586305859'],
            'rdb22':[rdb22, 'corrtest-L4kv2', '1586305859'],
            'rdb23':[rdb23, 'corrtest-c544M48s', '1596163435']}
            

    print('\n\nObservations availbale:')
    for key, value in RDBs.items():
        print(key, ': ', RDBs[key][0]['rdb'])
    ## ++++++++++++++++++++++++++++++++++++++++++++++
    # Select an observation:
    rdb_in   = RDBs[str(args.in_rdb)][0]
    descript = RDBs[str(args.in_rdb)][1]+'.%s'%timestr
    obsID    = RDBs[str(args.in_rdb)][2]

    rdbfile0    = rdb_in['rdb']
    do_image    = args.do_image     # enable katsdppipeline imaging
    do_srcfind  = args.do_srcfind      # enable source finding
    get_planes  = args.get_planes     # enable extraction of low, mid, high, and wide band images from cubeoid image
    image_targ  = args.image_targ
    image_bpcal = args.image_bpcal
    image_gcal  = args.image_gcal
    select_freq = args.select_freqs
    select_chan = args.select_chans

    if do_all_ObsIDs:
        obsID = ['1565560531', '1565564738', '1565646959', '1565651658']
    elif os.path.isfile(rdbfile0):
        obsID = os.path.basename(rdbfile0)[0:10]

    r_int    = list(range(0, 1000, 1))
    rnd_int  = rand.choice(r_int)
    descript = descript+str(rnd_int)
    print('\n >>> OBS ID = {}, descript/docker NAME: {}\nrdb_in: {}\n'.format(obsID,descript, rdb_in))
    
    ## ++++++++++++++++++++++++++++++++++++++++++++++

    # (!) 1 = True; 0 = False

    screen_on   = args.screen_on
    do_low      = args.do_low
    do_mid      = args.do_mid
    do_high     = args.do_high
    do_wide     = args.do_wide
    do_fullband = args.do_fullband
    selct_band  = args.selct_band
    toff        = args.toff
    stokes      = args.stokes

    use_smoothed_image = args.use_smoothed_image
    do_allscans        = args.Allscan_image          # image all scans
    do_singlescans     = args.per_scan_images        # image each scan separately

    ## ++++++++++++++++++++++++++++++++++++++++++++++
    if do_singlescans:
        tgt_scans  = rdb_in['tgt_scans']
        bcal_scans = rdb_in['bcal_scans']; gcal_scans = rdb_in['gcal_scans']
        do_allscans= False

    elif do_allscans:
        tgt_scans  = [rdb_in['tgt_scans']]
        bcal_scans = [rdb_in['bcal_scans']]; gcal_scans = [rdb_in['gcal_scans']]
        do_singlescans = False


    print('\n > do_allscans     = ', do_allscans)
    print(' > do_singlescans  = ', do_singlescans)
    ## ++++++++++++++++++++++++++++++++++++++++++++++

    if do_low:
        slct_chan  = low_band
        basename   = 'low_band'
        freqs      = low_fs
        maxFbw     = (bandwdth0)/mean(freqs)
    elif do_mid:
        slct_chan  = mid_band
        basename   = 'mid_band'
        freqs      = mid_fs
        maxFbw     = (bandwdth0)/mean(freqs)
    elif do_high:
        slct_chan  = hgh_band
        basename   = 'high_band'
        freqs      = hgh_fs
        maxFbw     = (bandwdth0)/mean(freqs)
    elif do_wide or selct_band:
        slct_chan  = wide_band
        basename   = 'wide_band'
        freqs      = all_fs
        maxFbw     = (bandwdth0)/mean(freqs)

    elif do_fullband:
        slct_chan  = full_band
        basename   = 'full_band'
        freqs      = full_fs
        maxFbw     = (bandwdth0)/mean(freqs)

    if screen_on:
        scrn_sessn = 'screen -L -m ' 
        #scrn_sessn = 'screen -L -d -m ' 
    if not screen_on:
        scrn_sessn = ''

    if select_chan:
        select_name = 'channels'
        select_vals = slct_chan
    else:
        select_name = 'freqrange'
        select_vals = freqs
        select_freq = True

    configs_disk = pwd0+"/configs/"
    parent_dsk   = pwd+"/tri_band/%s/"%basename

    if not (os.path.isdir(parent_dsk)):
        parent_dsk = pwd+'/'+basename+'/'
        try:
            os.mkdir(parent_dsk)
        except OSError as ose001:
            print(ose001)
            print('... passing')
            pass
    else:
        parent_dsk   = pwd+"/tri_band/%s/"%basename

    pols    = args.pol #'HH,VV'
    if pols in ['HH,VV', 'HH, VV']:
        pols_name = 'HHVV'
    elif pols in ['V', 'H']:
        pols_name = pols
    elif not pols:
        pols_name = 'HHVV'
        pols      = 'HH,VV'
        
    if image_targ or do_srcfind:
        slct_scans  = tgt_scans
    elif image_bpcal:
        slct_scans = bcal_scans
    elif image_gcal:
        slct_scans = gcal_scans
    else:
        slct_scans  = tgt_scans
    ## ++++++++++++++++++++++++++++++++++++++++++++++
    
    if do_image: 
        print(' \n> Starting katsdp imager ....')
        mydocker    = "sdp-docker-registry.kat.ac.za:5000/katsdpcontim:latest"
        off_docker  = "katsdpcontim-timeoffset" #"katacomb-time_offset:latest"
        corrsX      = 'cross'                # correlation products
        minflux0    = 0.003                  # Minimum Clean component (Jy) (not yet calibrated since there's no flux calibration)
        minFList0   = [1.0e-7]               # Minimum flux density to CLEAN per Self-cal cycle after first
        minFlux_psc = 1.0e-5                 # Min. peak phase self cal 
        minFlux_asc = 0.01*minFlux_psc       # Min. peak amp self cal 
        Niter0      = int(5.0e+5)            # Maximum # of I CLEAN comp.
        autoCen0    = 1.0e+25                # Auto center min flux density
        maxRealtime0= 86400.0                # Wall time after which MFImage will begin to shut down
        maxPSCLoop0 = 4                      # Max. number of phase selfcal loops
        targs0      = 'J0408-6545'           # target to image
        solAInt0    = 8./60.                 # A&P SC Solution interval (min)
        
        if float(toff) != 0.:
            run_offdocker = True
        elif float(toff) == 0.:
            run_offdocker = False
        
        for scanx in slct_scans:
            if do_singlescans:
                if float(toff) != 0.:
                    toff       = float(toff)
                    mount_disk = parent_dsk+"/FOV%s-deg.SCAN-%s-off-%.1fs-pol-%s"%(fov,scanx,toff,pols_name)
                else:
                    mount_disk = parent_dsk+"/FOV%s-deg.SCAN-%s-pol-%s"%(fov,scanx,pols_name)
                scanx      = [scanx]
            elif do_allscans:
                mount_disk = parent_dsk+"FOV%s-deg.multiscan-pol-%s"%(fov,pols_name)
                if type(scanx) == ndarray:
                    scanx = scanx.tolist()
            
            print('\n > imaging scan(s) (scanx=) ', scanx)
            
            if not (os.path.isdir(mount_disk)):
                os.mkdir(mount_disk)
                
            params_txt = open(mount_disk+'/docker-%s_run.txt'%(descript),'w')
            if os.path.isfile(rdbfile0):
                os.system( 'cp %s %s'%(rdbfile0,mount_disk) )
                rdbfile1 = mount_disk+'/'+os.path.basename(rdbfile0)
                rdbfile  = rdbfile1 #os.path.basename(rdbfile0)
            elif not os.path.isfile(rdbfile0):
                rdbfile = rdbfile0
            print(' > RDB in: ',rdbfile)
            os.system('cp %s*.yaml %s/'%(configs_disk, mount_disk))        # copy moded config files to the mounted disk
            print('> The following CONFIG files are in mounted disk: {}'.format(mount_disk))
            os.system('ls -lrt %s/*.yaml'%mount_disk)
            os.system('chmod a+w %s/.'%mount_disk)                       # allow docker user login to write to the mounted drive
            
            if run_offdocker:
                docker_in  = off_docker
                toff       = float(toff)
                RUN_comnd  = scrn_sessn+"docker run --user $USER --name {} --runtime=nvidia -v {}:/scratch {} continuum_pipeline.py {} -w /scratch -o /scratch --select \"scans={}; corrprods=\'{}\'; pol=\'{}\'; {}={}\" --uvblavg \"avgFreq=1; chAvg=4; FOV={}; maxInt={}\" --mfimage \"solPInt={};solAInt={};prtLv=2; FOV={}; Stokes=\'{}\'; maxFBW={}; minFlux={}; minFList={}; minFluxPSC={}; minFluxASC={}; Niter={}; autoCen={:.1f}; maxRealtime={}; maxPSCLoop={}\"".format(descript,mount_disk,docker_in,rdbfile,scanx,corrsX,pols,select_name,select_vals,fov,maxInt_time,solPInt,solAInt0,fov,stokes,maxFbw,minflux0,minFList0,minFlux_psc,minFlux_asc, Niter0, autoCen0,maxRealtime0,maxPSCLoop0) 
                
            else:
                docker_in  = mydocker
                RUN_comnd  = scrn_sessn+"docker run --user $USER --name {} --runtime=nvidia -v {}:/scratch {} continuum_pipeline.py {} -w /scratch -o /scratch --select \"scans={}; corrprods=\'{}\'; pol=\'{}\'; {}={}\" --uvblavg \"avgFreq=1; chAvg=4; FOV={}; maxInt={}\" --mfimage \"solPInt={};solAInt={};prtLv=2; FOV={}; Stokes=\'{}\'; maxFBW={}; minFlux={}; minFList={}; minFluxPSC={}; minFluxASC={}; Niter={}; autoCen={:.1f}; maxRealtime={}; maxPSCLoop={}\"".format(descript,mount_disk,docker_in,rdbfile,scanx,corrsX,pols,select_name,select_vals,fov,maxInt_time,solPInt,solAInt0,fov,stokes,maxFbw,minflux0,minFList0,minFlux_psc,minFlux_asc, Niter0, autoCen0,maxRealtime0,maxPSCLoop0) 
            
            if args.in_rdb == rdb12:
                RUN_comnd  = scrn_sessn+"docker run --user $USER --name {} --runtime=nvidia -v {}:/scratch {} continuum_pipeline.py {} -w /scratch -o /scratch --select \"targets=\'DEEP_2\'; scans=\'track\'; corrprods=\'{}\'; pol=\'{}\'; {}={}\" --uvblavg \"avgFreq=1; chAvg=4; FOV={}; maxInt={}\" --mfimage \"solPInt={};solAInt={};prtLv=2; FOV={}; Stokes=\'{}\'; maxFBW={}; minFlux={}; minFList={}; minFluxPSC={}; minFluxASC={}; Niter={}; autoCen={:.1f}; maxRealtime={}; maxPSCLoop={}\"".format(descript,mount_disk,docker_in,rdbfile,corrsX,pols,select_name,select_vals,fov,maxInt_time,solPInt,solAInt0,fov,stokes,maxFbw,minflux0,minFList0,minFlux_psc,minFlux_asc, Niter0, autoCen0,maxRealtime0,maxPSCLoop0) 
            
            if obsID in ['1573490035','1568825523', '1569521762', '1571007678', '1569087057', '1571058370']:
                if os.path.isfile(rdbfile0):
                    print((' >>>>>>\n Using local RDB: %s\n  >>>>>>'%rdbfile))
                    rdbfile = '/scratch/'+rdbfile
                    RUN_comnd  = scrn_sessn+"docker run --user $USER --name {} --runtime=nvidia -v {}:/scratch {} continuum_pipeline.py {} -w /scratch -o /scratch --select \"scans=\'track\'; corrprods=\'{}\'; pol=\'{}\'; {}={}\" --uvblavg \"avgFreq=1; chAvg=4; FOV={}; maxInt={}\" --mfimage \"solPInt={};solAInt={};prtLv=2; FOV={}; Stokes=\'{}\'; maxFBW={}; minFlux={}; minFList={}; minFluxPSC={}; minFluxASC={}; Niter={}; autoCen={:.1f}; maxRealtime={}; maxPSCLoop={}\"".format(descript,mount_disk,docker_in,rdbfile,corrsX,pols,select_name,select_vals,fov,maxInt_time,solPInt,solAInt0,fov,stokes,maxFbw,minflux0,minFList0,minFlux_psc,minFlux_asc, Niter0, autoCen0,maxRealtime0,maxPSCLoop0) 
                elif not os.path.isfile(rdbfile0):
                    RUN_comnd  = scrn_sessn+"docker run --user $USER --name {} --runtime=nvidia -v {}:/scratch {} continuum_pipeline.py {} -w /scratch -o /scratch --select \"targets=\'{}\';scans=\'track\'; corrprods=\'{}\'; pol=\'{}\'; {}={}\" --uvblavg \"avgFreq=1; chAvg=4; FOV={}; maxInt={}\" --mfimage \"solPInt={};solAInt={};prtLv=2; FOV={}; Stokes=\'{}\'; maxFBW={}; minFlux={}; minFList={}; minFluxPSC={}; minFluxASC={}; Niter={}; autoCen={:.1f}; maxRealtime={}; maxPSCLoop={}\"".format(descript,mount_disk,docker_in,rdbfile,targs0,corrsX,pols,select_name,select_vals,fov,maxInt_time,solPInt,solAInt0,fov,stokes,maxFbw,minflux0,minFList0,minFlux_psc,minFlux_asc, Niter0, autoCen0,maxRealtime0,maxPSCLoop0) 
            
            print(('\nDoing {} imaging. Selected channels: {}'.format(basename, slct_chan)))
            print(('RUN command :\n'+RUN_comnd+'\n'))
            print(f'****\n To track docker logs via terminal, do:\n\n docker logs -f {descript}\n\n****')
            
            params_txt.write('****\n To track docker logs via terminal, do:\n\n docker logs -f {}\n\n****'.format(descript))
            params_txt.write('Docker RUN command:\n'+RUN_comnd+'\n')
            os.system(RUN_comnd)
            
            if screen_on:
                print('Live screen sessions:')
                os.system('screen -ls')
                
            if do_allscans:
                break


    for obsIDx in obsID:
        if get_planes:
            print(' \n> extracting low, mid, and high band images ....')
            full_imdir = "/scratch2/slegodi/raw_obs/im_verify/tri_band/wide_band/"      # parent directory of full band image 
            slct_scans = tgt_scans
            srcimdir   = pwd+'/tri_band/images/wideband/'                               # source finding image directory
            
            for scanx in slct_scans:
                if float(toff) != 0.:
                    print('> time-offset: {} s'.format(toff))
                    if pols_name:
                        if pols_name == 'HHVV':
                            try:
                                imdir = pwd+'/tri_band/wide_band/FOV2.0-deg.SCAN-%d-off-%.1fs-pol-%s/'%(scanx,float(toff),pols_name)
                            except IOError as ioe002:
                                print('(!)', ioe002)
                                imdir = pwd+'/tri_band/wide_band/FOV2.0-deg.SCAN-%d-off-%.1fs/'%(scanx,float(toff))  # directory of katsdppipeline original images 
                                pass
                        else:
                            imdir = pwd+'/tri_band/wide_band/FOV2.0-deg.SCAN-%d-off-%.1fs-pol-%s/'%(scanx,float(toff),pols_name)  # directory of katsdppipeline original images
                    else:
                        imdir = pwd+'/tri_band/wide_band/FOV2.0-deg.SCAN-%d-off-%.1fs/'%(scanx,float(toff))  # directory of katsdppipeline original images 
                else:
                    imdir  = full_imdir+'FOV{}.0-deg/'.format(fov)  # directory of katsdppipeline original images 
                    
                filesx = []; subdirsx = []; pathsx = []
                
                for root, dirs, files in os.walk(imdir):
                    pathsx.append([root, dirs, files])
                
                print("Image directory for the funtion get_subbandImage():\n", imdir)
                for px in pathsx:
                    for subx in px:
                        try:
                            for xfile in os.listdir(subx):
                                xbase  = os.path.basename(xfile)
                                if xfile.startswith(obsIDx) and xfile.endswith('IClean.fits'):
                                    fitsname = subx+'/'+xfile
                                    xfilei   = xfile
                                    fitsfile = fitsname
                                elif xfile.endswith('.fits'):
                                    print('\nobsID: ', obsIDx)
                                    print('Other fits files found: ',subx+'/'+xfile)
                        except Exception as ex003:
                            #print ex003
                            pass
                
                try:
                    print('\nORIGINAL FILE: %s'%fitsfile)
                    hdr_frq = get_subbandImage(fitsfile, scanx, srcimdir,float(toff))
                except Exception as ioe001:
                    print('\n',ioe001)
                    
                    print('Possible file locations for file: ', xfilei)
                    os.system('ls -td $(locate %s)'%xfilei)
                    if do_all_ObsIDs:
                        pass
                    else:
                        fitsfile = str(input('Enter full path of cubeoid fitsfile: '))
                        print('\nORIGINAL FILE: %s'%fitsfile)
                        hdr_frq = get_subbandImage(fitsfile, scanx, srcimdir,float(toff))
                        
    for obsIDx in obsID:
        if do_srcfind:
            print(' \n> Starting source finding ....')
            test        = False
            not_test    = True
            show_model  = False  # show PYBDSF model fits
            
            slct_scans  = tgt_scans
            imdir0      = pwd+'/tri_band/images/'
            imdir00     = imdir0+'wideband/'
            fits00      = str(input(' >> full path to source finding fits image: ')) #'1564254957_continuum_image_J2147-8132_IClean.fits' # for testing
            fdata, hdr  = fits.getdata(fits00, header=True)
            targetField = hdr['OBJECT']
            
            # (!) fits files for twsting
            fitswide    = fits00[:-5]+'.wideband.fits'
            fitslow     = fits00[:-5]+'.lowband.fits'         # Assumes that file names match 
            fitsmid     = fits00[:-5]+'.midband.fits'
            fitshgh     = fits00[:-5]+'.hghband.fits'
            
            if not_test:
                if do_low:
                    band = 'lowband'
                elif do_mid:
                    band = 'midband'
                elif do_high:
                    band = 'hghband'
                elif do_wide:
                    band = 'wideband'
                elif selct_band:
                    band = 'freq%.3fGHz' %(hdr_frq/(1.e9))
                
                imdir    = imdir0+'/%s/'%band
                if not os.path.isdir(imdir):
                    print(' > making directory: ', imdir)
                    os.mkdir(imdir)
                
            elif test:
                imdir    = '/scratch2/slegodi/raw_obs/WSCLEAN.DEEP2.cubing/'     # for testing
                fitsname = 'DEEP_2.Iave.fits'                                    # test image
                fitsfile = imdir+fitsname
                print('Running on test image: ',fitsfile)
            
            for scanx in slct_scans:
                try:
                    fits00   = imdir00+obsIDx+'_continuum_image_%s_toff%.1f-scan%d.IClean.%s.fits'%(targetField,float(toff),scanx,band)
                    os.system('mv %s %s'%(fits00, imdir))
                except Exception as ex004:
                    print(ex004)
                    pass
                
                for srcfits in os.listdir(imdir):
                    if use_smoothed_image:
                        im_suffix = str(scanx)+'.IClean.'+band+'.smooth2sumss.fits'     # (smoothed images) Imagename suffix for source finding input image
                    if not use_smoothed_image:
                        im_suffix = str(scanx)+'.IClean.'+band+'.fits'                  # Imagename suffix for source finding input image
                        
                    
                    if srcfits.startswith(obsIDx) and srcfits.endswith(im_suffix) and '_toff%.1f'%float(toff) in srcfits:
                        print('> Loading fits "{}" into sourcefinding function'.format(imdir+srcfits))
                        fitsname = srcfits
                        try:
                            os.chdir(imdir)
                            fitsfile = fitsname
                            print('> fitsfile for sourcefinding: ', imdir+fitsfile)
                            find_sources(fitsfile,None, pwd, imdir, show_model)
                            os.chdir(pwd)
                        except Exception as re:
                            print(re, '\n\n')
                            pass
            os.chdir(pwd)

    #params_txt.close()
    process_end = time.time()
    duration_h  = (process_end - process_start)/3600.0
    duration_m  = (process_end - process_start)/60.0
    print("\n\n\n Total run time: %2.12f hrs (%2.12f min)" % (duration_h,duration_m)) 

