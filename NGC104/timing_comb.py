#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import string
#import correct as correct
from datetime import datetime
from scipy.interpolate import lagrange
from scipy import interpolate
from scipy.optimize import curve_fit
import pandas as pd
from astropy.stats import poisson_conf_interval
import scipy
import hawkeye as hawk

def load_data(dataname):
    path_Tuc='/Users/baotong/Desktop/period_Tuc/txt_all_obs_p90/'
    path = path_Tuc
    dataname = '{0}.txt'.format(dataname)
    epoch_file = path + 'epoch_src_' + dataname
    src_evt=np.loadtxt(path+dataname)
    epoch_info=np.loadtxt(epoch_file)
    useindex=np.where(epoch_info[:,3]>1000)
    epoch_info_use=epoch_info[useindex]
    useid=epoch_info_use[:,2]
    print(useid)
    src_evt_use =hawk.filter_obs(src_evt, useid)
    return (src_evt_use,epoch_info_use)

def main_process():
    (src_evt,epoch_info_use)=load_data(dataname='182')
    period =31055.9006;net_p = 0.8
    bin_len=500.
    time = src_evt[:, 0]
    hawk.plot_longT_V(src_evt=src_evt, bkg_file=None,epoch_info=epoch_info_use,iffold=True,p_test=period,shift=0.0 )
    plt.close()
    hawk.phase_fold(time=time,epoch_info=epoch_info_use,p_test=period,outpath=None,bin=50)

    lc=hawk.get_hist(time,len_bin=bin_len)
    T_tot=epoch_info_use[:,1][-1]-epoch_info_use[:,0][0]
    freq = np.arange(1 / T_tot, 0.5 / bin_len, 1 / (5* T_tot))
    freq = freq[np.where(freq > 1 / 50000.)]
    (FP, out_period, max_NormLSP)=hawk.get_LS(lc.time,lc.counts,freq=freq)
    print('Period=',format(out_period))

if __name__=='__main__':
    main_process()