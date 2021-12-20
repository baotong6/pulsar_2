#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import string
from scipy.interpolate import lagrange
from scipy import interpolate
from scipy.optimize import curve_fit
import pandas as pd
from astropy.stats import poisson_conf_interval
import scipy
import hawkeye as hawk
import NGC104.timing_comb as time_all

def get_HR(softband,hardband,src_evt):
    time_soft=hawk.filter_energy(src_evt[:,0],src_evt[:,1],softband)
    time_hard=hawk.filter_energy(src_evt[:,0],src_evt[:,1],hardband)
    cts_soft=len(time_soft)
    cts_hard=len(time_hard)
    return (cts_soft,cts_hard)

def load_data(dataname,ecf=90):
    path_Tuc='/Users/baotong/Desktop/period_Tuc/txt_all_obs_p{0}/'.format(ecf)
    path = path_Tuc
    dataname = '{0}.txt'.format(dataname)
    epoch_file = path + 'epoch_src_' + dataname
    src_evt=np.loadtxt(path+dataname)
    epoch_info=np.loadtxt(epoch_file)
    CR=hawk.plot_longT_V(src_evt=src_evt, bkg_file=None,epoch_info=epoch_info)
    # print(epoch_info[:,2])
    CR/=ecf/100.
    # (useid, epoch_info_use)=hawk.choose_obs(epoch_info,flux_info=CR,
    #                                         flux_filter=2e-3,expT_filter=1000,
    #                                         if_flux_high=True,if_expT_high=True,obsID=16527)
    (useid, epoch_info_use)=hawk.choose_obs(epoch_info,flux_info=CR,
                                            flux_filter=2e-3,expT_filter=1000,
                                            if_flux_high=True, if_expT_high=True,obsID=[15747])
    src_evt_use =hawk.filter_obs(src_evt, useid)
    print(useid)
    return (src_evt_use,epoch_info_use)

def plot_HR_lc():
    (src_evt_use,epoch_info_use)=load_data(dataname='211',ecf=90)
    len_bin=1000
    expT=epoch_info_use[:,-1];tstart=epoch_info_use[:,0][0];tstop=epoch_info_use[:,1][0]
    time = src_evt_use[:, 0]
    time_s = hawk.filter_energy(src_evt_use[:, 0], src_evt_use[:, 1], [500, 2000])
    time_h = hawk.filter_energy(src_evt_use[:, 0], src_evt_use[:, 1], [2000, 8000])
    lc_s=hawk.get_hist(time_s, len_bin, tstart=tstart, tstop=tstop)
    lc_h=hawk.get_hist(time_h, len_bin, tstart=tstart, tstop=tstop)
    plt.plot(lc_s.time,lc_s.counts)
    plt.show()
    plt.plot(lc_h.time,lc_h.counts)
    plt.show()
    HR=(lc_h.counts-lc_s.counts)/(lc_s.counts+lc_h.counts)
    print(HR)
    plt.scatter(lc_s.time,HR)
    plt.show()

if __name__=="__main__":
    plot_HR_lc()

