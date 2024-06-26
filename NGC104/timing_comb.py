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
import mcmc.gptLS as gptLS

def load_data(dataname,ecf=90,ifpath=None,ifobsID=[]):
    # path_Tuc='/Users/baotong/Desktop/period_Tuc/txt_startover/txt_all_obs_p{0}/'.format(ecf)
    # path_Tuc = '/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep4/'
    path_Tuc = '/Users/baotong/Desktop/period_Tuc/txt_startover/txt_all_obs_p{0}/'.format(ecf)
    path_M31='/Users/baotong/Desktop/M31XRB/M31HRC_txt/txt_all_obs_p90/'
    # path_Tuc='/Users/baotong/Desktop/period_omg/txt_all_obs_p{0}/'.format(ecf)
    # path_Tuc='/Users/baotong/Downloads/'
    # path_Tuc='/Users/baotong/Desktop/period_NGC3201/txt_all/txt_all_obs_p{0}/'.format(ecf)
    # path_Tuc = f'/Users/baotong/eSASS/data/raw_data/47_Tuc/txt/txt_merge_psf{ecf}_0.2_5/'
    if not ifpath:
        path = path_Tuc
    if ifpath:path=ifpath
    dataname = '{0}.txt'.format(dataname)
    epoch_file = path + 'epoch_src_' + dataname
    # dataname = 'transient_evt.txt'
    # epoch_file = path + 'SgrA_S_epoch.txt'
    src_evt=np.loadtxt(path+dataname)
    epoch_info=np.loadtxt(epoch_file)
    if epoch_info.ndim==1:epoch_info=np.array([epoch_info])
    if src_evt.ndim==1:src_evt=np.array([src_evt])

    if len(src_evt)<2:
        print('empty')
        return (np.array([]),np.array([]))
    else:
        CR=hawk.plot_longT_V(src_evt=src_evt, bkg_file=None,epoch_info=epoch_info)
        # CR/=ecf/100.
        # useobsID=epoch_info[:,2][0:50].astype('int')
        (useid, epoch_info_use)=hawk.choose_obs(epoch_info,flux_info=CR,
                                                flux_filter=280.,expT_filter=20000,
                                                if_flux_high=0, if_expT_high=1,obsID=ifobsID)
        print(epoch_info_use[:,3])
        if len(useid)==0 or len(epoch_info_use)==0: return ([],[])
        src_evt_use =hawk.filter_obs(src_evt, useid)
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

def main_process(path,dataname,period=None,ifobsID=None):
    ecf=90
    bin_len = 10
    net_p=0.95
    (src_evt_use,epoch_info_use)=load_data(dataname=dataname,ecf=90,ifpath=path,ifobsID=ifobsID)
    if src_evt_use==[] or epoch_info_use==[]:return None
    # lc=get_lc_frombkgimg(int(dataname),src_evt_use,epoch_info_use,ecf=90,bin_len=bin_len)
    figurepath = '/Users/baotong/Desktop/aas/GCall/figure/fold/'
    time=hawk.filter_energy(src_evt_use[:,0],src_evt_use[:,1],[20,8000])
    print('counts=',len(time))
    ## if filter time ##
    # time=hawk.filter_time_t1t2(time,t1=50000,t2=130000)
    # hawk.plot_Z2(time,freq=np.linspace(1/1e-4,1/500,10000))
    # hawk.plot_longT_V(src_evt=src_evt_use, bkg_file=None,epoch_info=epoch_info_use,
    #                   iffold=True,p_test=period,shift=0.,show=1)
    # hawk.phase_fold(time=time,epoch_info=epoch_info_use,net_percent=net_p,p_test=period,outpath=figurepath,bin=15,shift=0.,
    #                 label=dataname,text='Seq.{}'.format(dataname),save=0,show=1)
    hawk.phase_fold(time=time,epoch_info=epoch_info_use,net_percent=net_p,p_test=period,
                    outpath=figurepath,bin=30,shift=0.3,
                    label=dataname,textdef='#{0} (M28), P={1:.2f}s'.format(dataname,period,len(time)),
                    save=0,show=1)
    # plt.hist(time,bins=300,histtype='step')
    # plt.show()
    lc=hawk.get_hist(time,len_bin=bin_len,tstart=epoch_info_use[:,0][0],tstop=epoch_info_use[:,1][-1])
    # print(lc.time,lc.counts)
    T_tot=epoch_info_use[:,1][-1]-epoch_info_use[:,0][0]
    freq = np.arange(1 / T_tot, 0.5/ bin_len, 1 / (10* T_tot))
    freq = freq[np.where(freq > 1 / 2000.)]
    # outname=f'{dataname}_{ifobsID[0]}'
    # (FP, out_period, max_NormLSP)=hawk.get_LS(lc.time,lc.counts,freq,
    #                                           outpath='/Users/baotong/Desktop/period_M31XRB/fig_LS_HRC/',
    #                                           outname=f'{dataname}',save=0,show=1)

    # print(out_period)
    # (freq_grid, pgram)=gptLS.lomb_scargle(lc.time,lc.counts,freq=freq)
    # hawk.plot_singleobs_lc(lc,period=period,ifsin=0,
    #                        figurepath='/Users/baotong/Desktop/aas/pXS_Tuc/figure/',
    #                        shift=0.7,dataname=dataname,save=0,show=1)
if __name__=='__main__':
    # path_M31='/Users/baotong/Desktop/period_M31XRB/txt_all_obs_p90_HRC/'
    useid=[]
    path_GC = '/Users/baotong/Desktop/period_M28/txt_all_obs_p90/';useid=[]
    path_LW='/Users/baotong/Desktop/period_LW/txt_all_obs/'
    path_CDFS='/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep4/'
    path_NSC='/Users/baotong/Desktop/period/txt_all_obs_IG/'
    path_M31='/Users/baotong/Desktop/period_M31XRB/M31ACIS_txt/txt_all_obs_p90/'
    # path_M31='/Users/baotong/nustar/simulation/'
    dataname=21
    main_process(path=path_M31,dataname=dataname,period=22834.083,ifobsID=[])
    # for dataname in np.arange(214,215,1):
    #     for id in [1912,5925, 6177, 5926, 6202, 5927, 5928, 7283, 7284, 7285, 7286, 8526, 8527, 8528, 8529, 8530, 9825, 9826, 9827, 9828, 9829,10838,10683,10684,10882,10883,10884,10885,10886,11808,11809,12110,12111,12112,12113,12114,13178,13179,13180,13227,13228,13229,13230,13231,13278,13279,13280,13281]:
    #         main_process(path=path_M31,dataname=dataname,period=4521.41,ifobsID=[id])
