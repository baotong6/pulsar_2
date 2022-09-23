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
from scipy import optimize
from scipy import interpolate
from scipy.optimize import curve_fit
import pandas as pd
from astropy.stats import poisson_conf_interval
import scipy
from astropy.timeseries import LombScargle
from tkinter import _flatten
import rednoise as rednoise
def plot_lc(source_id,period=None):
    # path='/Users/baotong/Desktop/period_NGC6397/txt_all_obs_0.5_8/'
    path='/Users/baotong/Desktop/period_Tuc/txt_startover/txt_all_obs_p{0}/'.format(90)

    (evt_file,epoch_file)=rednoise.load_data(source_id,ecf=90,ifobsid=[953,955,2735,2736,2737,2738,16527,15747,16529,17420,15748,16528])
    obs_id=epoch_file[:,2]
    evt_id=evt_file[:,2]
    evt_list=evt_file[:,0]
    period =period
    bin_len=8000
    x_all = [];
    y_all = [];
    xerr_all = [];
    yerr_all = []
    for i in range(len(obs_id)):
        x=[];y=[];xerr=[];yerr=[]
        time=evt_list[np.where(evt_id==obs_id[i])]
        if len(time)<2:
            print('attention')
            continue

        k=0;
        while (time[0]+k*bin_len)<time[-1]:
            evt_temp=time[np.where((time>(time[0]+k*bin_len))&(time<(time[0]+(k+1)*bin_len)))]
            x.append(time[0]+(k+0.5)*bin_len)
            if (time[0]+(k+1)*bin_len)<time[-1]:
                exp_temp=bin_len
            else:
                exp_temp=time[-1]-(time[0]+k*bin_len)
            y.append(len(evt_temp)/exp_temp)
            xerr.append(0.5*bin_len)
            yerr.append(np.sqrt(len(evt_temp))/exp_temp)
            k+=1
        print(len(x))

        x_all.append(x);y_all.append(y);xerr_all.append(xerr);yerr_all.append(yerr)


    # print(x_all[16])
    return_time=x_all[0];return_flux=y_all[0]
    for i in range(1,len(x_all)):
        return_time=np.concatenate((return_time,x_all[i]))
        return_flux=np.concatenate((return_flux,y_all[i]))
    plt.figure(1)
    plt.errorbar(x_all[0], y_all[0], xerr=xerr_all[0], yerr=yerr_all[0])

    def fmax(x, a, b):
        return a * np.sin(x * 2 * np.pi / period) + b
    xsin=np.linspace(x_all[0][0],x_all[-1][-1],100000)
    ysin=0.001*np.sin(2*np.pi/period*xsin)+0.005
    # plt.plot(xsin,ysin)
    for i in range(1,len(x_all)):
        if (np.mod(x_all[i][0],period)/period)>(np.mod(x_all[i-1][-1],period)/period):
            gap=(np.mod(x_all[i][0],period)/period-np.mod(x_all[i-1][-1],period)/period)*period
        else:
            gap = (np.mod(x_all[i][0], period) / period - np.mod(x_all[i - 1][-1], period) / period) * period+period
        # print(gap)
        x_all[i] =x_all[i]- x_all[i][0] + x_all[i - 1][-1] + gap
        plt.errorbar(x_all[i],y_all[i],xerr=xerr_all[i],yerr=yerr_all[i])
        plt.fill_between([x_all[i-1][-1],x_all[i][0]],0.02,facecolor='yellow',alpha=0.2)
    x_all_list=x_all[0];y_all_list=y_all[0]
    for i in range(1,len(x_all)):
        x_all_list=np.concatenate((x_all_list,x_all[i]))
        y_all_list= np.concatenate((y_all_list, y_all[i]))

    fita, fitb = optimize.curve_fit(fmax, x_all_list, y_all_list, [0.5, 0.01])

    xperiod=x_all[0][0]+period*np.arange(0,int(500000/period),1)
    for i in range(len(xperiod)):
        plt.plot([xperiod[i],xperiod[i]],[0,0.02],'--',c='grey')
    plt.plot(xsin[np.where(xsin<x_all[-1][-1])], fmax(xsin[np.where(xsin<x_all[-1][-1])], fita[0], fita[1]),'--',c='r')
    plt.xlabel('time')
    plt.ylabel('counts rate')
    plt.show()
    return [return_time,return_flux]
    # plt.errorbar(x,y,xerr=xerr,yerr=yerr)
    # plt.show()
# plot_lc('212')

def get_LS(time, flux,freq):
    x = time
    y = flux
    LS = LombScargle(x, y,normalization = 'standard')
    # LS = LombScargle(x, y, dy, normalization='psd')
    power = LS.power(freq)
    FP=LS.false_alarm_probability(power.max(),samples_per_peak=5, nyquist_factor=5,minimum_frequency = freq[0], maximum_frequency = freq[-1],method='baluev')
    FP_99 = LS.false_alarm_level(0.01, samples_per_peak=10, nyquist_factor=5,minimum_frequency = freq[0], maximum_frequency = freq[-1],method='naive')
    FP_90 = LS.false_alarm_level(0.1, samples_per_peak=10, nyquist_factor=5, minimum_frequency=freq[0],
                                 maximum_frequency=freq[-1], method='naive')
    FP_68 = LS.false_alarm_level(0.32, samples_per_peak=10, nyquist_factor=5, minimum_frequency=freq[0],
                                 maximum_frequency=freq[-1], method='naive')


    plt.plot(freq, power)
    print(1./freq[np.where(power==np.max(power))])
    plt.plot([freq[0], freq[-1]], [FP_99, FP_99], '--')
    plt.plot([freq[0], freq[-1]], [FP_90, FP_90], '--')
    plt.plot([freq[0], freq[-1]], [FP_68, FP_68], '--')
    plt.show()
    # plt.savefig('/Users/baotong/Desktop/CDFS/fig_LS_ep{0}/{1}.eps'.format(k,dataname))
    # plt.close()

[time,flux]=plot_lc('312',period=48780.49)

border = 5000
vary = np.array([i for i in range(0, border)])
freq = 1 /20000. + vary * 1.e-8
# get_LS(time, flux,freq)