#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import string
from datetime import datetime
from scipy.interpolate import lagrange
from scipy import interpolate
from scipy.optimize import curve_fit
import pandas as pd
from astropy.stats import poisson_conf_interval
import scipy
import hawkeye as hawk

def load_data(dataname,ecf=90):
    path_Tuc='/Users/baotong/Desktop/period_Tuc/txt_all_obs_p{0}/'.format(ecf)
    # path_Tuc='/Users/baotong/Desktop/period_NGC3201/txt_all/txt_all_obs_p{0}/'.format(ecf)
    # path_Tuc = f'/Users/baotong/eSASS/data/raw_data/47_Tuc/txt/txt_merge_psf{ecf}_0.2_5/'
    path = path_Tuc
    dataname = '{0}.txt'.format(dataname)
    epoch_file = path + 'epoch_src_' + dataname
    src_evt=np.loadtxt(path+dataname)
    epoch_info=np.loadtxt(epoch_file)
    if epoch_info.ndim==1:epoch_info=np.array([epoch_info])
    CR=hawk.plot_longT_V(src_evt=src_evt, bkg_file=None,epoch_info=epoch_info)
    print(epoch_info[:,2])
    CR/=ecf/100.

    (useid, epoch_info_use)=hawk.choose_obs(epoch_info,flux_info=CR,
                                            flux_filter=30,expT_filter=10000,
                                            if_flux_high=0, if_expT_high=1,obsID=[16528])

    # [953,955,956, 2736, 3385,2738,16527,15747, 16529,15748]
    # [78, 953,   954,   955,   956,  2735,  3384,  2736,  3385,  2737,3386,  2738,  3387,
    # [16527,15747, 16529, 17420, 15748, 16528]
    src_evt_use =hawk.filter_obs(src_evt, useid)
    print(useid)
    return (src_evt_use,epoch_info_use)

def get_lc_frombkgimg(srcID,src_evt_use,epoch_info_use,ecf,bin_len):
    obsIDlist=epoch_info_use[:,2].astype('int')
    blank = np.zeros(len(obsIDlist)) + 1
    record=0
    for i in range(len(obsIDlist)):
        if src_evt_use.ndim < 2 or len(src_evt_use) < 10:
            blank[i] = 0
            continue
        obsid=obsIDlist[i]
        path='/Users/baotong/eSASS/data/raw_data/47_Tuc/txt/txt_psf{0}_{1}/'.format(ecf,obsid)
        src_info=np.loadtxt(path+'src_info.txt')
        bkg_cts_est_list = src_info[:, 2]
        bkg_cts_est = bkg_cts_est_list[np.where(src_info[:, 0] == srcID)][0]
        print(bkg_cts_est * bin_len / (epoch_info_use[:,1][i] - epoch_info_use[:,0][i]))
        time = src_evt_use[:, 0][np.where(src_evt_use[:,-1]==obsid)]
        lc = hawk.get_hist(time, len_bin=bin_len,tstart=epoch_info_use[:,0][i],tstop=epoch_info_use[:,1][i])
        lc.counts = lc.counts - bkg_cts_est * bin_len / (epoch_info_use[:,1][i] - epoch_info_use[:,0][i])
        lc.counts[np.where(lc.counts<0)]=0
        if i==0:
            lc_all=lc
            record=1
        elif record==0:
            lc_all=lc
        else:
            lc_all=lc_all.join(lc)

    return lc_all

def main_process():
    dataname='182'
    bin_len = 300
    (src_evt_use,epoch_info_use)=load_data(dataname=dataname,ecf=90)
    # lc=get_lc_frombkgimg(int(dataname),src_evt_use,epoch_info_use,ecf=90,bin_len=bin_len)
    figurepath = '/Users/baotong/Desktop/aas/pXS_Tuc/figure/'
    period =31180.4009/2
    net_p = 0.90

    time = src_evt_use[:, 0]

    time=hawk.filter_energy(src_evt_use[:,0],src_evt_use[:,1],[500,8000])
    hawk.plot_longT_V(src_evt=src_evt_use, bkg_file=None,epoch_info=epoch_info_use,iffold=True,p_test=period,shift=0.0)
    # plt.close()
    hawk.phase_fold(time=time,epoch_info=epoch_info_use,net_percent=net_p,p_test=period,outpath=figurepath,bin=100,shift=0.83,
                    label=dataname,text='Seq.162',save=0,show=1)

    # plt.hist(time,bins=300,histtype='step')
    # plt.show()
    lc=hawk.get_hist(time,len_bin=bin_len,tstart=epoch_info_use[:,0][0],tstop=epoch_info_use[:,1][-1])
    T_tot=epoch_info_use[:,1][-1]-epoch_info_use[:,0][0]
    freq = np.arange(1 / T_tot, 0.5 / bin_len, 1 / (10* T_tot))
    freq = freq[np.where(freq > 1 / 50000.)]

    figurepath='/Users/baotong/Desktop/aas/pXS_Tuc/figure/'
    (FP, out_period, max_NormLSP)=hawk.get_LS(lc.time,lc.counts,freq=freq,outpath=figurepath, outname=str(dataname),save=0,show=1)
    print('Period=',format(out_period))
    hawk.plot_singleobs_lc(lc,period=period,ifsin=0,figurepath='/Users/baotong/Desktop/aas/pXS_Tuc/figure/',dataname=dataname,save=0,show=1)
if __name__=='__main__':
    main_process()