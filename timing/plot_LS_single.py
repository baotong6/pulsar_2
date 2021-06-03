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

font1 = {'family': 'Normal',
         'weight': 'normal',
         'size': 16, }

def get_hist(t, len_bin,tstart=0,tstop=0):
    ###将输入的time信息，按照len_bin的长度输出为lc
    if tstart==0 and tstop==0:
        tstart=t[0]
        tstop=t[-1]
        tseg=tstop-tstart
    else:tseg=tstop-tstart

    t_test = t;dt=len_bin
    ev = EventList()
    ev.time = t_test
    lc_new = ev.to_lc(dt=dt, tstart=tstart-0.5*dt, tseg=tseg+0.5*dt)
    return lc_new

def get_hist_withbkg(t,t_bkg, len_bin,tstart,tstop):
    ###将输入的time信息，按照len_bin的长度输出为lc
    t_test = t;t_bkg_test=t_bkg;dt=len_bin
    t_bkg_test=np.delete(t_bkg_test,t_bkg_test<0)
    ev = EventList();ev_bkg=EventList()
    ev.time=t_test;ev_bkg.time=t_bkg_test
    lc_new = ev.to_lc(dt=dt, tstart=tstart , tseg=tstop - tstart)
    lc_bkg = ev_bkg.to_lc(dt=dt, tstart=tstart , tseg=tstop - tstart)
    lc_out=lc_new
    lc_out.counts=lc_new.counts-(1/12.)*lc_bkg.counts
    return lc_out
# def get_LS_myself(time,flux,freq):

def get_LS(time, flux,freq):
    x = time
    y = flux
    LS = LombScargle(x, y,normalization = 'standard')
    # LS = LombScargle(x, y, dy, normalization='psd')
    power = LS.power(freq)
    FP=LS.false_alarm_probability(power.max(),minimum_frequency = freq[0],maximum_frequency = freq[-1],method='baluev')
    FP_99 = LS.false_alarm_level(0.0027,minimum_frequency = freq[0], maximum_frequency = freq[-1],method='baluev')
    FP_90 = LS.false_alarm_level(0.05, minimum_frequency=freq[0],
                                 maximum_frequency=freq[-1], method='baluev')
    FP_68 = LS.false_alarm_level(0.32,minimum_frequency=freq[0],
                                 maximum_frequency=freq[-1], method='baluev')

    # if FP<0.01:print(dataname)
    # plt.title('{0},FP={1}'.format(dataname,FP))
    # plt.semilogx()
    plt.title('XID 19')
    plt.plot(freq, power)
    plt.semilogx()
    print(1. / freq[np.where(power == np.max(power))])
    print(np.where(power == np.max(power)))
    # if FP<0.01:print(1./freq[np.where(power==np.max(power))]);print(np.where(power==np.max(power)))
    # out_period=1./freq[np.where(power==np.max(power))]
    plt.plot([freq[0], freq[-1]], [FP_99, FP_99], '--')
    plt.plot([freq[0], freq[-1]], [FP_90, FP_90], '--')
    plt.plot([freq[0], freq[-1]], [FP_68, FP_68], '--')
    plt.show()
    # plt.savefig('/Users/baotong/Desktop/CDFS/fig_LS_ep{0}_ovsamp_5_baluev/{1}.eps'.format(k,dataname))
    # plt.close()
    return [FP,out_period]

def read_txt():
    # path='/Volumes/pulsar/CDFS/merge_data/timing/txt_all_obs_0.5_8/'
    # file=np.loadtxt(path+'XID19.txt')
    # epoch_file=np.loadtxt(path+'epoch_src_XID19.txt')

    path='/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep3/'
    file=np.loadtxt(path+'210.txt')
    epoch_file=np.loadtxt(path+'epoch_src_210.txt')

    time=file[:,0]
    tstart=epoch_file[:,0];tstop=epoch_file[:,1]
    # tstart=[epoch_file[0]];tstop=[epoch_file[1]]

    T_TOT = 2004251.6621800959;dt=100
    lc=get_hist(time,dt,tstart[0],tstop[-1])
    # lc_cut0=lc.truncate(tstart[0],tstop[0],method='time')
    lc_all=lc
    # i=0
    # while i < len(tstart):
    #     lc_temp=lc.truncate(tstart[i],tstop[i],method='time')
    #     lc_all=lc_all.join(lc_temp)
    #     i+=1
    def make_period_range(pmin, pmax, expT):
        P = [pmin]
        while P[-1] < pmax:
            dP = 0.5* P[-1] ** 2 / (expT - P[-1])
            P.append(P[-1] + dP)
        return np.array(P)

    freq=np.arange(1/T_TOT,1/200,1/(5*T_TOT))
    freq=freq[np.where(freq>1/20000.)]
    print(np.sum(lc_all.counts))

    # ps = Powerspectrum(lc_all,norm='leahy')
    # fig, ax1 = plt.subplots(1, 1, figsize=(9, 6), sharex=True)
    # ax1.loglog()
    # ax1.step(ps.freq, ps.power, lw=2, color='blue')
    # # ax1.plot([1 / 950.7, 1 / 950.7], [0, np.max(ps.power)], '--', linewidth=1)
    # ax1.set_ylabel("Frequency (Hz)", fontproperties=font1)
    # ax1.set_ylabel("Power (raw)", fontproperties=font1)
    # ax1.set_yscale('log')
    # ax1.tick_params(axis='x', labelsize=16)
    # ax1.tick_params(axis='y', labelsize=16)
    # ax1.tick_params(which='major', width=1.5, length=7)
    # ax1.tick_params(which='minor', width=1.5, length=4)
    # for axis in ['top', 'bottom', 'left', 'right']:
    #     ax1.spines[axis].set_linewidth(1.5)
    # plt.show()

    # avg_ps = AveragedPowerspectrum(lc_all, 60000,dt=lc_all.time[1]-lc_all.time[0],norm='leahy')
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




    get_LS(lc_all.time,lc_all.counts,freq)
if __name__ == '__main__':
    read_txt()