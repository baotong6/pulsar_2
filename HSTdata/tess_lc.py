#!/bin/bash
# -*- coding: utf-8 -*-
# written by Tong
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
from astropy.timeseries import LombScargle
font1 = {'family': 'Normal',
         'weight': 'normal',
         'size': 16, }
def get_LS(time, flux,error,freq):
    x = time
    y = flux
    dy=error
    FP=1.0
    LS = LombScargle(x, y,normalization = 'standard')
    # LS = LombScargle(x, y, normalization='psd')
    power = LS.power(freq)
    FP=LS.false_alarm_probability(power.max(),minimum_frequency = freq[0],maximum_frequency = freq[-1],method='baluev')
    FP_99 = LS.false_alarm_level(0.0027,minimum_frequency = freq[0], maximum_frequency = freq[-1],method='baluev')
    FP_95 = LS.false_alarm_level(0.05, minimum_frequency=freq[0],
                                 maximum_frequency=freq[-1], method='baluev')
    FP_68 = LS.false_alarm_level(0.32,minimum_frequency=freq[0],
                                 maximum_frequency=freq[-1], method='baluev')

    # if FP<0.01:print(dataname)
    plt.figure(1, (9, 6))
    plt.title('Lomb-Scargle Periodogram')
    plt.plot(freq, power)
    plt.semilogx()
    # print(1. / freq[np.where(power == np.max(power))])
    # print(np.where(power == np.max(power)))
    out_period=1./freq[np.where(power==np.max(power))]

    plt.plot([freq[0], freq[-1]], [FP_99, FP_99], '--')
    plt.plot([freq[0], freq[-1]], [FP_95, FP_95], '--')
    plt.plot([freq[0], freq[-1]], [FP_68, FP_68], '--')
    plt.text(freq[0], FP_99, '99%')
    plt.text(freq[0], FP_95, '95%')
    plt.text(freq[0], FP_68, '68%')

    plt.xlabel('Frequency (1/day)',font1)
    plt.ylabel('Normalized LSP',font1)
    plt.tick_params(labelsize=16)
    plt.show()
    return [FP,out_period]

def phase_fold(time,rate,error,period,binnumber=20,shift=0):
    ## 周期折叠，period为这个源的周期 ##
    turns=time*1/period-(time*1/period).astype('int')+shift
    fluxmean=np.zeros(binnumber);error_mean=np.zeros(binnumber)
    for i in range(len(fluxmean)):
        fluxmean[i]=np.mean(rate[np.where((turns<(i+1)/binnumber)&(turns>i/binnumber))])
        error_mean[i]=np.mean(error[np.where((turns<(i+1)/binnumber)&(turns>i/binnumber))])
    x=np.linspace(0,1,binnumber+1)
    x=x[:-1]

    x2=np.concatenate((x,x+1));y2=np.concatenate((fluxmean,fluxmean))
    y2_error=np.concatenate((error_mean,error_mean))
    ## 为了画出来漂亮一些，仅此而已
    ## 取平均值且取bin的画法
    plt.figure(1,(9,6))
    plt.errorbar(x2,y2,y2_error,fmt='ro')
    plt.xlabel('Phase',font1)
    plt.ylabel('Count rate',font1)
    plt.tick_params(labelsize=16)
    plt.show()

    ## 把所有采样点都放进去，这个源看起来有phase变化 ##
    plt.figure(1,(9,6))
    plt.errorbar(turns,rate,error,fmt='ro')
    plt.xlabel('Phase',font1)
    plt.ylabel('Count rate',font1)
    plt.tick_params(labelsize=16)
    plt.show()

def read_lc(path,filename):
    lc_file=fits.open(path+filename)
    time=lc_file[1].data['TIME']
    flux=lc_file[1].data['PDCSAP_FLUX']
    flux_err=lc_file[1].data['PDCSAP_FLUX_ERR']

    ##找出并去掉nan值##
    NAN_index=np.isnan(flux)
    valid_index=np.where(NAN_index==False)[0]
    time=time[valid_index];flux=flux[valid_index];flux_err=flux_err[valid_index]
    ##看起来最后两个点有问题##
    time=time[:-3];flux=flux[:-3];flux_err=flux_err[:-3]

    # print(time)
    plt.figure(1,(9,6))
    plt.title('Lightcurve',font1)
    plt.xlabel('Time (BTJD)',font1)
    plt.ylabel('PDCSAP_FLUX',font1)
    plt.tick_params(labelsize=16)
    plt.errorbar(time,flux,flux_err,fmt='bo')
    plt.show()
    T_tot=time[-1]-time[0]
    dt=time[1]-time[0]
    ##is that right?##
    freq=np.arange(1/T_tot,0.5/dt,1/(5*T_tot))
    [FP,period]=get_LS(time,flux,flux_err,freq)
    print(period)
    phase_fold(time,flux,flux_err,period,binnumber=100,shift=0.)

if __name__ == "__main__":
    read_lc('/Users/baotong/Desktop/HST/','hlsp_tess-spoc_tess_phot_0000000002771358-s0019_tess_v1_lc.fits')

