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
import linecache
from astropy.timeseries import LombScargle
import funcs_timing as funcs

def plot_LS():
    bin_len=100
    path='/Users/baotong/eSASS/data/raw_data/47_Tuc/txt/txt_merge_0.5_5/'
    figurepath='/Users/baotong/eSASS/data/raw_data/47_Tuc/fig_LS/LS_bin{0}/'.format(bin_len)
    srcID=np.arange(1,889,1)
    FP_all=[];out_period_all=[]
    for i in range(len(srcID)):
        srcevt=np.loadtxt(path+'{0}.txt'.format(srcID[i]))
        if len(srcevt)==0:continue
        time=srcevt[:,0];energy=srcevt[:,1]
        (TSTART, TSTOP, OBSID, exptime) = funcs.read_epoch(path + 'epoch_src_{0}.txt'.format(srcID[i]))
        lc=funcs.get_hist(time,len_bin=bin_len,tstart=TSTART[0], tstop=TSTOP[-1])
        x = lc.time;flux = lc.counts

        print(i)
        T_tot = lc.time[-1] - lc.time[0]
        freq = np.arange(1 / T_tot, 0.5 / bin_len, 1 / (5 * T_tot))
        freq=freq[np.where(freq > 1 / 10000.)]
        [FP, out_period,max_NormLSP] = funcs.get_LS(x, flux, freq, outpath=figurepath,outname=str(srcID[i]))
        FP_all.append(FP);out_period_all.append(out_period)
    FP_all=np.array(FP_all);out_period_all=np.array(out_period_all)
    LS_info=np.column_stack((srcID,FP_all,out_period_all))
    np.savetxt(figurepath+'LS_info_bin{0}.txt'.format(bin_len),LS_info,fmt='%10d %10.5f %10.5f')
if __name__=='__main__':
    plot_LS()