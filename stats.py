from matplotlib import pyplot as plt
import pandas as pd
from math import*
import os, sys, time, atpy, argparse, csv, subprocess
from astropy.io import fits
from numpy import*
import random as r
#####################################################################################################
mymarkers   = ['o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X']
def quick_stats(csv_file, base=None, do_rms=1, do_mean=1, do_max=1, do_min=0, inmark=None):
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
        plt.subplot(312)
        plt.title('Stokes I max for each mode')
        plt.scatter(frq, Fmx, marker=inmark, label='{}-{}, avg: {:.2e}$\pm{:.2e}$ Jy/bm'.format(mode,base, fd_max, std_max))
        plt.legend(loc=0)
        plt.ylabel('Max [Jy/bm]', fontsize=20)
        #plt.xlabel('Frequency [Hz]', fontsize=20)
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
        
    
    if do_mean:
        plt.subplot(313)
        plt.title('Stokes I mean for each mode')
        plt.scatter(frq, avg, marker=inmark, label='{}-{}, avg: {:.2e}$\pm{:.2e}$ Jy/bm'.format(mode,base, fd_mean, std_avg))
        plt.legend(loc=0)
        plt.ylabel('Mean [Jy/bm]', fontsize=20)
        plt.xlabel('Frequency [Hz]', fontsize=20)
        plt.grid(True, ls = ':')
        #plt.tight_layout()
    
    
    return(fd_mean, fd_max, rms_avg)

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

def Get_posnRMS(self, min_axs, maj_axs, SNR):
    beamsize = sqrt( min_axs**2 + maj_axs**2 )
    pos_RMS  = beamsize/(2.*SNR)
    
    return (pos_RMS)

def JulxMay_corr_test_offsets():
    """
    Plot MeerKAT L-band to L-band RA vs Dec offsets.
    Offsets are between two meerKAT observations of the Jul and May 2020 
    correlator imaging test observations. 
    
    """
    catdir = "/home/samuel/SARAO/imagin_varification/2020-correlator-test/July/"
    tLxL   = atpy.Table(catdir+'L4kJul-x-L4kMay.tbl')
    tUxU   = atpy.Table(catdir+'U4kJul-x-U4kMay.tbl')  
    
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
    plt.show()
    
    return ()


#####################################################################################################
indir  = '/home/samuel/SARAO/imagin_varification/2020-correlator-test/May/analysis/'
algo   = 'fit-half' #'classic'
xmatch = indir+"SUMMSxL4kWB-cat1xL32kWB-cat2xL32kNB-cat3xU4kWB-cat4xU32kWB-cat5.fits"
xmatch2= indir+"SUMMSxL4kWBcmc2xL4kWBcmc1.tbl" #SUMMS x 1590129959 x 1586305859
xmatch3= indir+"L4kWBcmc2xL4kWBcmc1.tbl" #1590129959 x 1586305859
xmatch4= indir+"L4kWBcmc2xL4kWBcmc1-10asecradius.tbl"

cats1  = {'1':'L4kWB','2':'L32kWB','3':'L32kNB','4':'U4kWB','5':'U32kWB'}
cats2  = {'2':'L4kWB-May (1590129959)', '3':'L4kWB-April (1586305859)'}
cats3  = {'1':'L4kWB-May (1590129959)', '2':'L4kWB-April (1586305859)'}
Nxcats = 6

Tx  = atpy.Table(xmatch) 
Tx2 = atpy.Table(xmatch2)
Tx3 = atpy.Table(xmatch3)
Tx4 = atpy.Table(xmatch4)

maj1 = Tx4['Maj_1']; min1 = Tx4['Min_1']                                                                                                                                                                                                                              
maj2 = Tx4['Maj_2']; min2 = Tx4['Min_2']                                                                                                                                                                                                                              
snr1 = Tx4['Peak_flux_1']/Tx4['Isl_rms_1']                                                                                                                                                                                                                            
snr2 = Tx4['Peak_flux_2']/Tx4['Isl_rms_2']  
posrms1 = Get_posnRMS(min1, max1, snr1)
posrms2 = Get_posnRMS(min2, maj2, snr2)

plt.figure()
plt.scatter(posrms1*3600., posrms2*3600.)
plt.plot(posrms1*3600., posrms1*3600., 'r', label='y=x')
plt.legend(loc=0)
plt.ylabel('posn RMS1 [arcsec]', fontsize=15)
plt.xlabel('posn RMS2 [arcsec]', fontsize=15)
plt.grid(True, ls = ':')
plt.show() 

ra0, dec0   = Tx['_RAJ2000'], Tx['_DEJ2000']
ra02, dec02 = Tx2['_RAJ2000'], Tx2['_DEJ2000']

incat = cats1
inTab = Tx

plt.figure()
for key, val in incat.items():
    i           = int(key)
    ra00, dec00 = inTab['_RAJ2000'], inTab['_DEJ2000']
    try:
        ra_i, dec_i  = inTab['RA_%d'%i], inTab['DEC_%d'%i]
        dRA_i,dDec_i = RA_Dec_offsets(ra00, dec00, ra_i, dec_i)
        inmark       = r.choice(mymarkers)
        mode         = val
        dRA_avg      = nanmean(dRA_i)*3600
        dDec_avg     = nanmean(dDec_i)*3600
        
        plt.title('Position offsets between each mode and SUMMS', fontsize=20)
        plt.scatter(dRA_i*3600, dDec_i*3600,
                    marker=inmark,
                    label='{}: {:.2e}\", {:.2e}\"'.format(mode, dRA_avg, dDec_avg))
        plt.legend(loc=0)
        plt.ylabel('dDEC [arcsec]', fontsize=20)
        plt.xlabel('dRA [arcsec]', fontsize=20)
        plt.grid(True, ls = ':')
    except Exception as ex:
        print(' (!) ',ex)
        pass

for key, val in cats3.items(): 
    i           = int(key) 
    try: 
        dRA_i,dDec_i = RA_Dec_offsets(Tx3['RA_%d'%i], Tx3['DEC_%d'%i], Tx3['RA_%d'%(i+1)], Tx3['DEC_%d'%(i+1)]) 
        dRA_avg      = nanmean(dRA_i)*3600 
        dDec_avg     = nanmean(dDec_i)*3600 
        #plt.title('May L4kWB (1590129959) cross-matched with the same mode from April (1586305859)') 
        plt.scatter(dRA_i*3600, dDec_i*3600, marker='+', alpha=0.5, 
                    label='L4kWB-May x April: {:.2e}\", {:.2e}\"'.format(dRA_avg, dDec_avg)) 
        plt.legend(loc=0) 
        plt.ylabel('dDEC [arcsec]', fontsize=20) 
        plt.xlabel('dRA [arcsec]', fontsize=20) 
        plt.grid(True, ls = ':') 
         
    except Exception as ex: 
        print ('>>> (!)', ex) 
        pass                   

#plt.show()

plt.figure()
for fi in os.listdir(indir):
    if fi.startswith('15') and fi.endswith('channelised_imstats.csv'):
        fi = indir+fi
        try:
            quick_stats(fi) 
        except Exception as ex:
            print (f' ! : {ex}')
            print (f' ! input file: {fi}')
            pass
plt.show()
