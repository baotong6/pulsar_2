'''
Author: baotong && baotong@smail.nju.edu.cn
Date: 2023-06-20 14:58:30
LastEditors: baotong && baotong@smail.nju.edu.cn
LastEditTime: 2024-01-17 15:24:26
FilePath: /pulsar/hawkeye/Z2.py
Description: 

Copyright (c) 2023 by baotong, All Rights Reserved. 
'''
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
    useidlw=[7285]
    (src_evt_use,epoch_info_use)=tdata.load_data(dataname=dataname,ecf=90,ifpath=path,ifobsID=useidlw)
    if src_evt_use.ndim<2:time=[]
    else:time=src_evt_use[:,0]
    return time

def get_Z2(dataname,pathin=None,ep=1):
    warnlist=[]
    nharm = 1
    events = EventList()
    time=read_data(dataname=dataname,path=pathin)
    # time=np.loadtxt(pathin+dataname)[:,0]
    print('cts=',len(time))
    if len(time)<2:return 0
    events.time=time
    df_min=1e-8
    frequencies=np.arange(1/10000,1/500,df_min)
    nbin=20
    ntrial = (frequencies[-1] - frequencies[0]) / df_min
    freq, zstat = z_n_search(events.time, frequencies, nbin=nbin, nharm=nharm)
    z_detlev = z2_n_detection_level(n=1, epsilon=0.001, ntrial=len(freq)*1)
    z_detlev2 = z2_n_detection_level(n=1, epsilon=0.01, ntrial=len(freq)*1)
    z_detlev3 = z2_n_detection_level(n=1, epsilon=0.1, ntrial=len(freq)*1)
    # ---- PLOTTING --------
    plt.figure()
    plt.plot(freq, (zstat - nharm), label='$Z_2$ statistics')

    # plt.scatter([0.0001727,0.0001727*2,0.0001727*3,0.0001727*4,0.0001727*5,0.0001727*6,0.0001727*7,0.0001727*8,0.0001727*9,0.0001727*10],
    #             [400,300,200,100,50,50,50,40,40,40],
    #             s=100,marker='^',color='black')
    plt.axhline(z_detlev - nharm, label='$Z^2_1$ det. lev. 99.9%',color='r',linestyle='--')
    plt.axhline(z_detlev2 - nharm, label='$Z^2_1$ det. lev. 99%',color='g',linestyle='--')
    plt.axhline(z_detlev3 - nharm, label='$Z^2_1$ det. lev. 90%',color='g',linestyle='--')
    # plt.plot(freq, efstat - nbin + 1, color='gray', label='EF statistics', alpha=0.5)
    print('Period=',1/freq[zstat.argmax()])
    # plt.axvline(1/period, linestyle='--',color='r', lw=1, alpha=0.5, label='Correct frequency')
    plt.xlim([frequencies[0], frequencies[-1]])
    plt.xlabel('Frequency (Hz)')
    plt.ylabel(r'$Z^2$ power - d.o.f.')
    plt.legend()
    plt.semilogx()
    plt.show()
    # plt.figure(figsize=(15, 5))
    # plt.plot(freq, (zstat - nharm), label='$Z_2$ statistics')
    # plt.plot(freq, efstat - nbin + 1, color='gray', label='EF statistics', alpha=0.5)

    # plt.axvline(1/period, color='r', lw=3, alpha=0.5, label='Correct frequency')
    # plt.xlabel('Frequency (Hz)')
    # plt.ylabel('Statistics - d.o.f. (Zoom)')
    # plt.ylim([-15, 15])
    # _ = plt.xlim([frequencies[0], frequencies[-1]])
    # figurepath=f'/Users/baotong/Desktop/CDFS/figure_Z2/ep{ep}/'
    # plt.savefig(figurepath+f'Z2_{dataname}.pdf',bbox_inches='tight', pad_inches=0.05)
    if np.max(zstat)>z_detlev3:
        return dataname
    else:
        return 0
    # plt.show()
if __name__=='__main__':
    get_Z2('74',pathin='/Users/baotong/Desktop/period_M31XRB/txt_all_obs_p90_HRC/',ep=1)
    # read_data()
    # for ep in [1,2,3,4]:
    #     warnlist = []
    #     for i in range(1,1000):
    #         a=get_Z2(str(i),ep)
    #         if a>0:
    #             warnlist.append(a)
    #
    #     np.savetxt(f'/Users/baotong/Desktop/CDFS/figure_Z2/ep{ep}/warnlist90.txt',warnlist)
    # get_Z2(dataname='110',pathin='/Users/baotong/Desktop/period_M31XRB/M31HRC_txt/txt_all_obs_p90/')
    # get_Z2(dataname='evt_allobs_GTI100.txt',pathin='/Users/baotong/swift/ztf_src1/output/')

    # id=['00035071001','00035071002','00035071003','00035071004','00035071005',
    #     '00035071006','00035071007','00035071008','00035071009','00035071010',
    #     '00035071011']
    # for i in range(len(id)):
    #     get_Z2(dataname='srcevt.txt',pathin='/Users/baotong/swift/ztf_src1/output/{}/'.format(id[i]))



