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
# import scipy.signal.lombscargle as LombScargle
from astropy.stats import poisson_conf_interval
import stingray as sr
from stingray.events import EventList
from stingray.lightcurve import Lightcurve
from stingray import Lightcurve, Crossspectrum, sampledata,Powerspectrum,AveragedPowerspectrum
from stingray.simulator import simulator, models

def get_hist(t, len_bin):
    ###将输入的time信息，按照len_bin的长度输出为lc
    t_test = t-t[0];dt=len_bin
    ev = EventList()
    ev.time = t_test
    lc_new = ev.to_lc(dt=dt, tstart=ev.time[0] - 0.5 * dt, tseg=ev.time[-1] - ev.time[0])
    return lc_new
def get_hist_withbkg(t,t_bkg, len_bin):
    ###将输入的time信息，按照len_bin的长度输出为lc
    t_test = t-t[0];t_bkg_test=t_bkg-t[0];dt=len_bin
    t_bkg_test=np.delete(t_bkg_test,t_bkg_test<0)
    ev = EventList();ev_bkg=EventList()
    ev.time=t_test;ev_bkg.time=t_bkg_test
    lc_new = ev.to_lc(dt=dt, tstart=ev.time[0] - 0.5 * dt, tseg=ev.time[-1] - ev.time[0])
    lc_bkg = ev_bkg.to_lc(dt=dt, tstart=ev.time[0] - 0.5 * dt, tseg=ev.time[-1] - ev.time[0])
    lc_out=lc_new
    lc_out.counts=lc_new.counts-(1/12.)*lc_bkg.counts
    return lc_out
# def get_LS_myself(time,flux,freq):

def get_LS(time, flux,freq,dataname,k):
    x = time
    y = flux
    # dy=np.sqrt(y)
    # plt.scatter(x,y)
    # plt.show()

    # LS = LombScargle(x, y, dy = 1, normalization = 'standard', fit_mean = True,
    #                  center_data = True).power(freq, method = 'cython')
    LS = LombScargle(x, y,normalization = 'standard')
    # LS = LombScargle(x, y, dy, normalization='psd')
    power = LS.power(freq)
    FP=LS.false_alarm_probability(power.max(),minimum_frequency = freq[0],maximum_frequency = freq[-1],method='baluev')
    FP_99 = LS.false_alarm_level(0.0027,minimum_frequency = freq[0], maximum_frequency = freq[-1],method='baluev')
    FP_90 = LS.false_alarm_level(0.05, minimum_frequency=freq[0],
                                 maximum_frequency=freq[-1], method='baluev')
    FP_68 = LS.false_alarm_level(0.32,minimum_frequency=freq[0],
                                 maximum_frequency=freq[-1], method='baluev')

    if FP<0.01:print(dataname)
    plt.title('{0},FP={1}'.format(dataname,FP))
    # plt.semilogx()
    plt.plot(freq, power)
    plt.semilogx()
    print(1. / freq[np.where(power == np.max(power))])
    print(np.where(power == np.max(power)))
    if FP<0.01:print(1./freq[np.where(power==np.max(power))]);print(np.where(power==np.max(power)))
    out_period=1./freq[np.where(power==np.max(power))]
    plt.plot([freq[0], freq[-1]], [FP_99, FP_99], '--')
    plt.plot([freq[0], freq[-1]], [FP_90, FP_90], '--')
    plt.plot([freq[0], freq[-1]], [FP_68, FP_68], '--')
    plt.show()
    # plt.savefig('/Users/baotong/Desktop/CDFS/fig_LS_ep{0}_ovsamp_5_baluev/{1}.eps'.format(k,dataname))
    # plt.close()
    return [FP,out_period]
# ### uneven sampled freq ###
def get_freq_unsamp(exptime):
    p_unsamp = []
    freq_unsamp = []
    p_unsamp.append(1e-6 * exptime)
    while p_unsamp[-1] < exptime:
        if p_unsamp[-1] < 0.3 * exptime:
            p_unsamp.append(p_unsamp[-1] + p_unsamp[-1] ** 2 / (2 * exptime * 2))
        elif 0.3 * exptime < p_unsamp[-1] < 0.5 * exptime:
            p_unsamp.append(p_unsamp[-1] + p_unsamp[-1] ** 2 / (3 * exptime * 2))
        else:
            p_unsamp.append(p_unsamp[-1] + p_unsamp[-1] ** 2 / (4 * exptime * 2))
    p_unsamp = np.array(p_unsamp)
    p_unsamp=p_unsamp[np.where((p_unsamp>200)&(p_unsamp<10000))]
    freq_unsamp = 1.0 / p_unsamp
    freq_unsamp=np.sort(freq_unsamp)
    return freq_unsamp
    ### uneven sampled freq ###
# path='/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep1/'
# border = 500000

bin_len=100.0
# vary = np.array([i for i in range(0, border)])
# source_id=np.linspace(1+599,1055,1055-599)
# source_id=source_id.astype(int)

# highp_id_ep1=[89,168,643]
# highp_id_ep2=[51,141]
# highp_id_ep3=[89,106,127,810,867]
# highp_id_ep4=[37,220]
highp_id_ep3=[19]
source_id=highp_id_ep3
# exptime = time[-1] - time[0]
# freq = 1 /(5*86400.) + vary * 1.e-8
# flux=get_hist(time,bin_len)
# x=np.arange(bin_len/2.,(time[-1]-time[0])+bin_len/2.,bin_len)
# get_LS(x,flux,freq,str(source_id[0]))
def plot_CDFS_ep_LS(k_num):
# k_num=[1,2,3,4]
    ##写入txt文本方便plot
    for k in k_num:
        FP_print = [];
        out_period_print = [];
        source_name = [];
        src_cts = [];
        bkg_cts = [];
        cts_rate = []
        path = '/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep{0}/'.format(k)
        for i in range(len(source_id)):
            time=np.loadtxt(path+'{0}.txt'.format(source_id[i]))
            bkg_time=np.loadtxt(path+'{0}_bkg.txt'.format(source_id[i]))
            epoch = np.loadtxt(path + 'epoch_src_{0}.txt'.format(source_id[i]))
            if len(epoch)==0:continue
            if len(np.shape(epoch)) == 1:
                exptime = epoch[3]
            else:
                exptime = np.sum(epoch[:, 3])

            if len(time)<4:continue

            if len(bkg_time)<2:bkg_time=[]

            else:
                time=time[:,0]
                bkg_time=bkg_time[:,0]
                src_cts.append(len(time));bkg_cts.append(len(bkg_time))
                T_tot=time[-1]-time[0]
                # freq = get_freq_unsamp(exptime)
                freq=np.arange(1/T_tot,0.5/bin_len,1/(5*T_tot))
                freq=freq[np.where(freq > 1 / 20000.)]
                # print(freq[100]-freq[99])
                lc=get_hist_withbkg(time,bkg_time,bin_len)
                counts=np.sum(lc.counts)
                cts_rate.append(counts/exptime)
                # flux=get_hist(time,bin_len)
                x=lc.time;flux=lc.counts
                (FP,out_period)=get_LS(x,flux,freq,str(source_id[i]),k)
                FP_print.append(FP);out_period_print.append(out_period);source_name.append(source_id[i])

        result = np.column_stack((source_name,FP_print, out_period_print,src_cts,bkg_cts,cts_rate))
        # np.savetxt('/Users/baotong/Desktop/CDFS/fig_LS_ep{0}_ovsamp_5_baluev/LS_result_{0}_tail.txt'.format(k),result,fmt='%10d %15.10f %15.5f %10d %10d %15.10f')
plot_CDFS_ep_LS([3])
