# We will search for pulsations over a range of frequencies around the known pulsation period.
#!/bin/bash
# -*- coding: utf-8 -*-
# written by Tong
import numpy as np
import matplotlib.pyplot as plt
from stingray.pulse.search import epoch_folding_search, z_n_search
from stingray.pulse.search import search_best_peaks
from stingray.stats import fold_detection_level, z2_n_detection_level
from matplotlib import ticker
from astropy.io import fits
from scipy.optimize import curve_fit
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.timeseries import LombScargle
import stingray as sr
from stingray.events import EventList
from stingray.lightcurve import Lightcurve
from stingray import Lightcurve, Crossspectrum, sampledata,Powerspectrum,AveragedPowerspectrum
import hawkeye as hawk
import NGC104.plot_pXS as plot_pXS
import GC.localCVs as localCVs
import NGC104.CV_model as CV_model
import NGC104.timing_comb as tdata
def read_data(dataname,path):
    ecf=90
    bin_len = 12.8
    net_p=0.94
    useidlw=[9502, 9500, 9501,9854,9503,9892,9893, 9504]
    (src_evt_use,epoch_info_use)=tdata.load_data(dataname=dataname,ecf=90,ifpath=path,ifobsID=[])
    if src_evt_use.ndim<2:time=[]
    else:time=src_evt_use[:,0]

    return time

def get_Z2(dataname,ep=1):
    warnlist=[]
    nharm = 1

    events = EventList()
    path_M31='/Users/baotong/Desktop/M31XRB/M31HRC_txt/txt_all_obs_p90/';useid=[]
    path_GC = '/Users/baotong/Desktop/period_M28/txt_all_obs_p90/';useid=[]
    path_LW='/Users/baotong/Desktop/period_LW/txt_all_obs/'
    path_CDFS=f'/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep{ep}/'

    time=read_data(dataname=dataname,path=path_CDFS)
    if len(time)<2:return 0
    events.time=time
    df_min=1e-7
    frequencies=np.arange(2e-5,1e-3,df_min)
    nbin=20
    ntrial = (frequencies[-1] - frequencies[0]) / df_min
    freq, zstat = z_n_search(events.time, frequencies, nbin=nbin, nharm=nharm)
    z_detlev = z2_n_detection_level(n=1, epsilon=0.001, ntrial=len(freq)*1000)
    z_detlev2 = z2_n_detection_level(n=1, epsilon=0.01, ntrial=len(freq)*1000)
    z_detlev3 = z2_n_detection_level(n=1, epsilon=0.1, ntrial=len(freq)*1000)
    # ---- PLOTTING --------
    plt.figure()
    plt.plot(freq, (zstat - nharm), label='$Z_2$ statistics')
    plt.axhline(z_detlev - nharm, label='$Z^2_1$ det. lev. 99.9%',color='r',linestyle='--')
    plt.axhline(z_detlev2 - nharm, label='$Z^2_1$ det. lev. 99%',color='g',linestyle='--')
    plt.axhline(z_detlev3 - nharm, label='$Z^2_1$ det. lev. 90%',color='g',linestyle='--')
    # plt.plot(freq, efstat - nbin + 1, color='gray', label='EF statistics', alpha=0.5)

    # plt.axvline(1/period, linestyle='--',color='r', lw=1, alpha=0.5, label='Correct frequency')
    plt.xlim([frequencies[0], frequencies[-1]])
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Statistics - d.o.f.')
    plt.legend()
    plt.semilogx()
    # plt.figure(figsize=(15, 5))
    # plt.plot(freq, (zstat - nharm), label='$Z_2$ statistics')
    # plt.plot(freq, efstat - nbin + 1, color='gray', label='EF statistics', alpha=0.5)

    # plt.axvline(1/period, color='r', lw=3, alpha=0.5, label='Correct frequency')
    # plt.xlabel('Frequency (Hz)')
    # plt.ylabel('Statistics - d.o.f. (Zoom)')
    # plt.ylim([-15, 15])
    # _ = plt.xlim([frequencies[0], frequencies[-1]])
    figurepath=f'/Users/baotong/Desktop/CDFS/figure_Z2/ep{ep}/'
    plt.savefig(figurepath+f'Z2_{dataname}.pdf',bbox_inches='tight', pad_inches=0.05)
    if np.max(zstat)>z_detlev3:
        return int(dataname)
    else:
        return 0
    # plt.show()
if __name__=='__main__':
    # read_data()
    for ep in [1,2,3,4]:
        warnlist = []
        for i in range(1,1000):
            a=get_Z2(str(i),ep)
            if a>0:
                warnlist.append(a)

        np.savetxt(f'/Users/baotong/Desktop/CDFS/figure_Z2/ep{ep}/warnlist90.txt',warnlist)