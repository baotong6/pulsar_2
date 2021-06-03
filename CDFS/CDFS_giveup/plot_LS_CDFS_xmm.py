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
    return (T,E)

def filter_obs(src_evt,bkg_evt,useid):
    src_evt_use = src_evt[np.where(src_evt[:-1] == useid[0])[0]]
    bkg_evt_use = bkg_evt[np.where(bkg_evt[:-1] == useid[0])[0]]
    i=1
    while i < len(useid):
        id=useid[i]
        src_evt_use_temp=src_evt[np.where(src_evt[:-1]==id)[0]]
        bkg_evt_use_temp=bkg_evt[np.where(bkg_evt[:-1]==id)[0]]
        src_evt_use = np.concatenate((src_evt_use, src_evt_use_temp))
        bkg_evt_use = np.concatenate((bkg_evt_use, bkg_evt_use_temp))
        i+=1
    return (src_evt_use,bkg_evt_use)
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
def get_hist_withbkg(t,t_bkg, len_bin,tstart=0,tstop=0):
    ###将输入的time信息，按照len_bin的长度输出为lc
    if tstart==0 and tstop==0:
        tstart=t[0]
        tstop=t[-1]
        tseg=tstop-tstart
    else:
        tseg = tstop - tstart
    t_test = t;t_bkg_test=t_bkg;dt=len_bin
    t_bkg_test=np.delete(t_bkg_test,t_bkg_test<0)
    ev = EventList();ev_bkg=EventList()
    ev.time=t_test;ev_bkg.time=t_bkg_test
    lc_new = ev.to_lc(dt=dt, tstart=tstart-0.5*dt, tseg=tseg+0.5*dt)
    lc_bkg = ev_bkg.to_lc(dt=dt, tstart=tstart-0.5*dt, tseg=tseg+0.5*dt)
    lc_bkg_level=np.mean(lc_bkg.counts)
    lc_out=lc_new
    lc_out.counts=lc_new.counts-lc_bkg.counts
    # lc_out.counts = lc_new.counts
    # print(lc_bkg.counts)
    return lc_out

# def get_LS_myself(time,flux,freq):

def get_LS(time, flux,freq,dataname='default',outpath='/Users/baotong/Desktop/'):
    x = time
    y = flux
    # dy=np.sqrt(y)
    # plt.scatter(x,y)
    # plt.show()
    FP=1.0
    # LS = LombScargle(x, y, dy = 1, normalization = 'standard', fit_mean = True,
    #                  center_data = True).power(freq, method = 'cython')
    LS = LombScargle(x, y,normalization = 'standard')
    # LS = LombScargle(x, y, normalization='psd')
    power = LS.power(freq)
    FP=LS.false_alarm_probability(power.max(),minimum_frequency = freq[0],maximum_frequency = freq[-1],method='baluev')
    FP_99 = LS.false_alarm_level(0.01,minimum_frequency = freq[0], maximum_frequency = freq[-1],method='baluev')
    FP_95 = LS.false_alarm_level(0.05, minimum_frequency=freq[0],
                                 maximum_frequency=freq[-1], method='baluev')
    FP_68 = LS.false_alarm_level(0.32,minimum_frequency=freq[0],
                                 maximum_frequency=freq[-1], method='baluev')
    #
    # if FP<0.01:print(dataname)
    plt.title('{0},FP={1}'.format(dataname,FP))
    # plt.semilogx()
    plt.plot(freq, power)
    plt.semilogx()
    # print(1. / freq[np.where(power == np.max(power))])
    # print(np.where(power == np.max(power)))
    # if FP<0.01:print(1./freq[np.where(power==np.max(power))]);print(np.where(power==np.max(power)))
    out_period=1./freq[np.where(power==np.max(power))[0]]
    plt.plot([freq[0], freq[-1]], [FP_99, FP_99], '--')
    plt.plot([freq[0], freq[-1]], [FP_95, FP_95], '--')
    plt.plot([freq[0], freq[-1]], [FP_68, FP_68], '--')
    plt.text(freq[0],FP_99,'99%')
    plt.text(freq[0], FP_95, '95%')
    plt.text(freq[0], FP_68, '68%')
    plt.show()
    # plt.savefig(outpath+'{0}_15sec.eps'.format(dataname))
    # plt.close()
    return [FP,out_period]

bin_len=100
highp_id_ep3=[289]
# source_id=np.arange(1,340,1)
source_id=highp_id_ep3
def plot_CDFS_ep_LS(ep,det,band):
    FP_print = [];
    out_period_print = [];
    source_name = [];
    src_cts = [];
    bkg_cts = [];
    cts_rate = []

    path = '/Users/baotong/Desktop/CDFS/xmm_CDFS/xmm_txt/epoch{0}_txt/'.format(ep)
    outpath=path+'fig_LS/'
    for i in range(len(source_id)):
        srcevt = np.loadtxt(path + 'NID{0}_{1}_10sec.txt'.format(source_id[i], det))
        epoch = np.loadtxt(path + 'epoch_NID{0}_{1}_10sec.txt'.format(source_id[i], det))
        if len(srcevt)<5 or len(epoch)==0:continue
        if epoch.ndim == 1: epoch = np.array([epoch])
        time=srcevt[:,0];energy=srcevt[:,1]
        (time, energy) = filter_energy(time, energy, band)
        src_cts.append(len(time))
        source_name.append(int(source_id[i]))
        cts_rate.append(len(time)/np.sum(epoch[:, 3]))
        lc = get_hist(time, bin_len, tstart=epoch[:,0][0], tstop=epoch[:,1][-1])
        T_tot = lc.time[-1] - lc.time[0]
        freq=np.arange(1/T_tot,0.5/bin_len-0.0002,1/(5*T_tot))
        freq=freq[np.where(freq > 1 / 20000.)]
        x=lc.time;flux=lc.counts
        (FP, out_period) = get_LS(x, flux, freq, dataname=source_id[i],outpath=outpath)
        FP_print.append(FP)
        out_period_print.append(out_period)

    result = np.column_stack((source_name, FP_print, out_period_print, src_cts, cts_rate))
    print(out_period_print)
    # np.savetxt(path+'LS_result_mos_ep{0}_only_srccts_15sec.txt'.format(ep), result,
    #            fmt='%10d %15.10f %15.5f %10d %15.10f')

plot_CDFS_ep_LS(1,'mos',[500,8000])
plot_CDFS_ep_LS(2,'mos',[500,8000])
plot_CDFS_ep_LS(3,'mos',[500,8000])
plot_CDFS_ep_LS(4,'mos',[500,8000])
plot_CDFS_ep_LS(5,'mos',[500,8000])
plot_CDFS_ep_LS(6,'mos',[500,8000])

#             return None
# def plot_CDFS_ep_LS(source_id,mode):
#     FP_print = [];
#     out_period_print = [];
#     source_name = [];
#     src_cts = [];
#     bkg_cts = [];
#     cts_rate = []
#     # path = '/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep{0}/'.format(k)
#     path='/Users/baotong/Desktop/CDFS/xmm_txt/'
#     def get_lc_mode(det,band):
#         srcevt=np.loadtxt(path+'{0}_{1}_all_obs_10sec.txt'.format(source_id,det))
#         # bkgevt=np.loadtxt(path+'bkg_XID{0}_{1}_all_obs_10sec.txt'.format(source_id,det))
#         epoch = np.loadtxt(path + 'epoch_XID{0}_{1}_all_obs_10sec.txt'.format(source_id,det))
#         # useid=epoch[:,2][8:12]
#         # (srcevt,bkgevt)=filter_obs(srcevt,bkgevt,useid)
#         if len(epoch)==0:
#             return None
#         if len(np.shape(epoch)) == 1:
#             exptime = epoch[3]
#         else:
#             exptime = np.sum(epoch[:, 3])
#         if len(srcevt)<4:
#             return None
#         if len(bkgevt)<2:bkgevt=[]
#         else:
#             time=[];bkg_time=[];energy=[];bkg_energy=[]
#             srcevt[:, -1] = srcevt[:, -1].astype('int');
#             # bkgevt[:, -1] = bkgevt[:, -1].astype('int')
#             for k in range(len(useid)):
#                 time =np.concatenate((time, srcevt[:, 0][np.where(srcevt[:, -1] == useid[k])]))
#                 energy=np.concatenate((energy, srcevt[:, 1][np.where(srcevt[:, -1] == useid[k])]))
#                 # bkg_time = np.concatenate((bkg_time,bkgevt[:, 0][np.where(bkgevt[:, -1] == useid[k])]))
#                 # bkg_energy=np.concatenate((bkg_energy,bkgevt[:, 1][np.where(bkgevt[:, -1] == useid[k])]))
#             ## 看一下不同波段的lc如何 ##
#             (time,energy)=filter_energy(time,energy,band)
#             # (bkg_time,bkg_energy)=filter_energy(bkg_time,bkg_energy,band)
#             src_cts=len(time)
#             # bkg_cts=len(bkg_time)
#             print('src_counts={0}'.format(src_cts))
#             # print('bkg_counts={0}'.format(bkg_cts))
#             # lc=get_hist_withbkg(time,bkg_time,bin_len)
#             # lc = get_hist(time, bin_len,tstart=epoch[:,0][8],tstop=epoch[:,1][11])
#             # print(epoch[:,0][8],epoch[:,1][11])
#             # print(time[0],time[-1])
#             lc = get_hist(time, bin_len, tstart=time[0], tstop=time[-1])
#             # lc = get_hist_withbkg(time, bkg_time, bin_len,tstart=epoch[:,0][8],tstop=epoch[:,1][11])
#             return (lc,exptime)
#     (lc0,exptime0)=get_lc_mode(mode[0],[500,8000])
#     # (lc1, exptime1) = get_lc_mode(mode[1], [1000, 12000])
#     # (lc2, exptime2) = get_lc_mode(mode[2], [1000, 12000])
#     lc=lc0;exptime=exptime0
#     # lc_all=lc2;exptime=exptime1
#     # print(len(lc0.time))
#     # print(len(lc1.time))
#     # print(len(lc2.time))
#
#     # lc_all.counts+=lc0.counts
#     # lc_all.counts+=lc1.counts
#
#     T_tot=lc.time[-1]-lc.time[0]
#     freq=np.arange(1/T_tot,0.5/bin_len,1/(5*T_tot))
#     freq=freq[np.where(freq > 1 / 20000.)]
#     counts=np.sum(lc.counts)
#     cts_rate.append(counts/exptime)
#     x=lc.time;flux=lc.counts
#     (FP, out_period) = get_LS(x, flux, freq, source_id)
#         # FP_print.append(FP);out_period_print.append(out_period);source_name.append(source_id)
#
#     # result = np.column_stack((source_name,FP_print, out_period_print,src_cts,bkg_cts,cts_rate))
#     # np.savetxt('/Users/baotong/Desktop/CDFS/fig_LS_ep{0}_ovsamp_5_baluev/LS_result_{0}_tail.txt'.format(k),result,fmt='%10d %15.10f %15.5f %10d %10d %15.10f')
# plot_CDFS_ep_LS(source_id,['pn'])
