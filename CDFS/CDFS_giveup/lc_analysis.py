#!/bin/bash
# -*- coding: utf-8 -*-
# written by Tong
# plot the phase-folded light curve from txt file (version for xmm)
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
from astropy.timeseries import LombScargle
import stingray
from stingray import Lightcurve, Powerspectrum, AveragedPowerspectrum

font1 = {'family': 'Normal',
         'weight': 'normal',
         'size': 16, }

def filter_energy(time,energy,band):
    T=time
    E=energy
    i=0
    if len(time)!=len(energy):
        print('error')
        return None
    else:
        while i <len(E):
            if E[i]<band[0] or E[i]>band[1]:
                E=np.delete(E,i)
                T=np.delete(T,i)
            else:
                i+=1
    return T
def get_LS(time, flux,freq):
    x = time
    y = flux

    # LS = LombScargle(x, y, dy = 1, normalization = 'standard', fit_mean = True,
    #                  center_data = True).power(freq, method = 'cython')
    LS = LombScargle(x, y,normalization = 'standard')
    FP=0
    power = LS.power(freq)
    FP=LS.false_alarm_probability(power.max(),minimum_frequency = freq[0], maximum_frequency = freq[-1],method='baluev')
    FP_99 = LS.false_alarm_level(0.0027, minimum_frequency = freq[0], maximum_frequency = freq[-1],method='baluev')
    FP_90 = LS.false_alarm_level(0.05,  minimum_frequency=freq[0],
                                 maximum_frequency=freq[-1], method='baluev')
    FP_68 = LS.false_alarm_level(0.32, minimum_frequency=freq[0],
                                 maximum_frequency=freq[-1], method='baluev')
    plt.figure(1,(9,6))
    plt.plot([freq[0], freq[-1]], [FP_99, FP_99], '--')
    plt.plot([freq[0], freq[-1]], [FP_90, FP_90], '--')
    plt.plot([freq[0], freq[-1]], [FP_68, FP_68], '--')
    # plt.plot([1 /950.7, 1 / 950.7], [0, np.max(power)], '--', linewidth=1)
    plt.title('FP={0}'.format(FP))
    plt.semilogx()
    # plt.xlim(1000.,1500.)
    plt.plot(freq, power)
    print(1./freq[np.where(power==np.max(power))])

    plt.show()
    res=1e5*power
    res=np.round(res,2)
    return [FP, 1. / freq[np.where(power == np.max(power))],np.max(power),res]

def plot_pds(time,flux):
    lc = Lightcurve(time, flux)
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    ax.plot(lc.time, lc.counts, lw=2, color='blue')
    ax.set_xlabel("Time (s)", fontproperties=font1)
    ax.set_ylabel("Counts (cts)", fontproperties=font1)
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    ax.tick_params(which='major', width=1.5, length=7)
    ax.tick_params(which='minor', width=1.5, length=4)
    plt.show()

    ps = Powerspectrum(lc,norm='leahy')
    fig, ax1 = plt.subplots(1, 1, figsize=(9, 6), sharex=True)
    ax1.loglog()
    ax1.step(ps.freq, ps.power, lw=2, color='blue')
    ax1.plot([1 / 950.7, 1 / 950.7], [0, np.max(ps.power)], '--', linewidth=1)
    ax1.set_ylabel("Frequency (Hz)", fontproperties=font1)
    ax1.set_ylabel("Power (raw)", fontproperties=font1)
    ax1.set_yscale('log')
    ax1.tick_params(axis='x', labelsize=16)
    ax1.tick_params(axis='y', labelsize=16)
    ax1.tick_params(which='major', width=1.5, length=7)
    ax1.tick_params(which='minor', width=1.5, length=4)
    for axis in ['top', 'bottom', 'left', 'right']:
        ax1.spines[axis].set_linewidth(1.5)
    plt.show()

    # avg_ps = AveragedPowerspectrum(lc, 500,dt=lc.time[1]-lc.time[0],norm='leahy')
    # print("Number of segments: %d" % avg_ps.m)
    # fig, ax1 = plt.subplots(1, 1, figsize=(9, 6))
    # ax1.loglog()
    # ax1.step(avg_ps.freq, avg_ps.power, lw=2, color='blue')
    # ax1.set_xlabel("Frequency (Hz)", fontproperties=font1)
    # ax1.set_ylabel("Power (raw)", fontproperties=font1)
    # ax1.set_yscale('log')
    # ax1.tick_params(axis='x', labelsize=16)
    # ax1.tick_params(axis='y', labelsize=16)
    # ax1.tick_params(which='major', width=1.5, length=7)
    # ax1.tick_params(which='minor', width=1.5, length=4)
    # for axis in ['top', 'bottom', 'left', 'right']:
    #     ax1.spines[axis].set_linewidth(1.5)
    # plt.show()


def read_SAS_lc(srcName,obsID,dt,useinst='pn',obsmode='single'):
    # 1,2,3分别代表mos1,mos2,pn的light curve，也可以加起来用，记为_all;
    # 根据实际情况来决定lomb-scargle的输入
    def read_lc_oneobs(obsID,dt):
        path = '/Volumes/pulsar/xmmCDFS/{0}/cal/'.format(obsID)
        print(obsID)
        os.chdir(path)
        mode = ['mos1', 'mos2', 'pn']
        filename1 = mode[0] + '_' + srcName + '_lccorr_bin{0}_0.5_8.lc'.format(int(dt));
        filename2 = mode[1] + '_' + srcName + '_lccorr_bin{0}_0.5_8.lc'.format(int(dt));
        filename3 = mode[2] + '_' + srcName + '_lccorr_bin{0}_0.5_8.lc'.format(int(dt))

        # filename1 = mode[0] + '_' + srcName + '_src_lc_bin{0}_1_8.lc'.format(int(dt));
        # filename2 = mode[1] + '_' + srcName + '_src_lc_bin{0}_1_8.lc'.format(int(dt));
        # filename3 = mode[2] + '_' + srcName + '_src_lc_bin{0}_1_8.lc'.format(int(dt))

        lc1 = fits.open(filename1);
        lc2 = fits.open(filename2);
        lc3 = fits.open(filename3)
        time1 = lc1[1].data['TIME'];
        rate1 = lc1[1].data['RATE']
        time2 = lc2[1].data['TIME'];
        rate2 = lc2[1].data['RATE']
        time3 = lc3[1].data['TIME'];
        rate3 = lc3[1].data['RATE']
        rate1 = np.nan_to_num(rate1);
        rate2 = np.nan_to_num(rate2);
        rate3 = np.nan_to_num(rate3)
        rate1[np.where(rate1 < 0)] = 0;
        rate2[np.where(rate2 < 0)] = 0;
        rate3[np.where(rate3 < 0)] = 0
        rate1[np.where(rate1 > 0.5)] = 0
        rate2[np.where(rate2 > 0.5)] = 0
        rate3[np.where(rate3 > 0.5)] = 0
        if useinst == 'pn':
            time_all = time3;
            rate_all = rate3
        if useinst == 'mos1+mos2+pn':
            time_all = time3;
            rate_all = rate1 + rate2 + rate3
        return [time_all, rate_all]
    if obsmode=='single':
        (time_all,rate_all)=read_lc_oneobs(obsID,dt)
        print(np.sum(rate_all*dt))
        T_tot = time_all[-1] - time_all[0]
        freq = np.arange(1. / T_tot, 0.5 / dt, 1. / (10 * T_tot))
        freq = freq[np.where(freq > 1 / 10000.)]
        get_LS(time_all, rate_all, freq)
        plot_pds(time_all, rate_all)
    if obsmode=='multiple':
        time_all=[];rate_all=[]
        for ID in obsID:
            (time, rate) = read_lc_oneobs(ID, dt)
            time_all=np.concatenate((time_all,time))
            rate_all=np.concatenate((rate_all,rate))
        print('counts={0}'.format(np.sum(rate_all*dt)))
        T_tot = time_all[-1] - time_all[0]
        freq = np.arange(1. / T_tot, 0.5 / dt, 1. / (10 * T_tot))
        freq=freq[np.where(freq>1/10000.)]
        get_LS(time_all, rate_all, freq)
        plot_pds(time_all, rate_all)

obsList=["0108060401","0108060501","0108060601","0108060701","0108061801","0108061901","0108062101",
         "0108062301","0555780101","0555780201","0555780301","0555780401","0555780501","0555780601",
         "0555780701","0555780801","0555780901","0555781001","0555782301","0604960101","0604960201",
         "0604960301","0604960401","0604961101","0604961201","0604960701","0604960501","0604961301",
         "0604960601","0604960801","0604960901","0604961001","0604961801"]
# read_SAS_lc('XID19','0604960601',dt=50,useinst='pn',obsmode='single')
if __name__=='__main__':
    read_SAS_lc('XID19',obsList[0:2],  dt=100,useinst='mos1+mos2+pn',obsmode='multiple')
    read_SAS_lc('XID19',obsList[2:8],  dt=100,useinst='mos1+mos2+pn',obsmode='multiple')
    read_SAS_lc('XID19',obsList[8:12], dt=100,useinst='mos1+mos2+pn',obsmode='multiple')
    read_SAS_lc('XID19',obsList[12:19],dt=100,useinst='mos1+mos2+pn',obsmode='multiple')
    read_SAS_lc('XID19',obsList[19:23],dt=100,useinst='mos1+mos2+pn',obsmode='multiple')
# read_SAS_lc('XID19',obsList[28],dt=100,useinst='mos1+mos2+pn',obsmode='single')
# read_SAS_lc('XID19','0108060501',dt=50,useinst='mos1+mos2+pn',obsmode='single')
# for i in range(len(obsList)):
#     read_SAS_lc('XID89',obsList[i],dt=50,useinst='mos1+mos2+pn',obsmode='single')
