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
    lc_out=lc_new
    lc_out.counts=lc_new.counts-(1/12.)*lc_bkg.counts
    # lc_out.counts = lc_new.counts
    # print(lc_bkg.counts)
    return lc_out
# def get_LS_myself(time,flux,freq):


def get_LS(time, flux,freq,dataname,k):
    x = time
    y = flux
    # dy=np.sqrt(y)
    # plt.scatter(x,y)
    # plt.show()
    FP=1.0
    LS = LombScargle(x, y,dy=None,normalization = 'standard')
    # LS = LombScargle(x, y, normalization='psd')
    power = LS.power(freq)
    # print(power.max())
    FP=LS.false_alarm_probability(power.max(),minimum_frequency = freq[0],maximum_frequency = freq[-1],method='baluev')
    FP_99 = LS.false_alarm_level(0.01,minimum_frequency = freq[0], maximum_frequency = freq[-1],method='baluev')
    FP_95 = LS.false_alarm_level(0.05, minimum_frequency=freq[0],
                                 maximum_frequency=freq[-1], method='baluev')
    FP_68 = LS.false_alarm_level(0.32,minimum_frequency=freq[0],
                                 maximum_frequency=freq[-1], method='baluev')

    # if FP<0.01:print(dataname)
    plt.title('Epoch {2}: XID={0},FAP={1}'.format(dataname,np.round(FP,4),k),font1)
    # plt.semilogx()
    plt.plot(freq, power)
    # plt.plot([1/2492.73,1/2492.73],[0,np.max(power)],'--',linewidth=1)
    plt.semilogx()
    print(1. / freq[np.where(power == np.max(power))])
    # print(np.where(power == np.max(power)))
    # if FP<0.01:print(1./freq[np.where(power==np.max(power))]);print(np.where(power==np.max(power)))
    out_period=1./freq[np.where(power==np.max(power))]
    plt.plot([freq[0], freq[-1]], [FP_99, FP_99], '--')
    plt.plot([freq[0], freq[-1]], [FP_95, FP_95], '--')
    plt.plot([freq[0], freq[-1]], [FP_68, FP_68], '--')
    plt.text(freq[0], FP_99, 'FAP 99%',font1)
    plt.text(freq[0], FP_95, '95%',font1)
    plt.text(freq[0], FP_68, '68%',font1)
    plt.xlabel('Frequency (Hz)',font1)
    plt.ylabel('Normalized LS Periodogram',font1)
    plt.tick_params(labelsize=16)
    plt.show()
    # plt.savefig('/Users/baotong/Desktop/CDFS/fig_LS_ep{0}_ovsamp_5_baluev/{1}.eps'.format(k,dataname))
    # plt.savefig('/Users/baotong/Desktop/CDFS/fig_LS_ep{0}_samp_1_baluev/{1}.eps'.format(k,dataname))
    # plt.close()
    return [FP,out_period]

def plot_CDFS_ep_LS(k_num):
    ##写入txt文本方便plot
    for k in k_num:
        FP_print = [];
        out_period_print = [];
        source_name = [];
        src_cts = [];
        bkg_cts = [];
        cts_rate = []
        path = '/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep{0}/'.format(k)
        print(source_id)
        # path='/Users/baotong/Desktop/CDFS/xmm_txt/'
        for i in range(len(source_id)):
            srcevt=np.loadtxt(path+'{0}.txt'.format(source_id[i]))
            bkgevt=np.loadtxt(path+'{0}_bkg.txt'.format(source_id[i]))
            epoch = np.loadtxt(path + 'epoch_src_{0}.txt'.format(source_id[i]))
            if len(srcevt) < 4 or len(epoch) == 0: continue
            if epoch.ndim == 1: epoch = np.array([epoch])
            #
            # srcevt=np.loadtxt(path+'XID{0}_pn_all_obs.txt'.format(source_id[i]))
            # bkgevt=np.loadtxt(path+'bkg_XID{0}_pn_all_obs.txt'.format(source_id[i]))
            # epoch = np.loadtxt(path + 'epoch_XID{0}_pn_all_obs.txt'.format(source_id[i]))
            useid=epoch[:,2][13:23]
            # (srcevt,bkgevt)=filter_obs(srcevt,bkgevt,useid)
            if len(epoch)==0:continue
            if len(np.shape(epoch)) == 1:
                exptime = epoch[3]
            else:
                exptime = np.sum(epoch[:, 3])

            if len(srcevt)<4:continue

            if len(bkgevt)<2:bkgevt=[]

            else:
                time=[];bkg_time=[];energy=[];bkg_energy=[]
                srcevt[:, -1] = srcevt[:, -1].astype('int');
                bkgevt[:, -1] = bkgevt[:, -1].astype('int')

                for id_k in range(len(useid)):
                    time =np.concatenate((time, srcevt[:, 0][np.where(srcevt[:, -1] == useid[id_k])]))
                    energy=np.concatenate((energy, srcevt[:, 1][np.where(srcevt[:, -1] == useid[id_k])]))
                    bkg_time = np.concatenate((bkg_time,bkgevt[:, 0][np.where(bkgevt[:, -1] == useid[id_k])]))
                    bkg_energy=np.concatenate((bkg_energy,bkgevt[:, 1][np.where(bkgevt[:, -1] == useid[id_k])]))
                ## 看一下不同波段的lc如何 ##
                (time,energy)=filter_energy(time,energy,[500,8000])
                (bkg_time,bkg_energy)=filter_energy(bkg_time,bkg_energy,[500,8000])

                src_cts.append(len(time));bkg_cts.append(len(bkg_time))
                # print('src_counts={0}'.format(src_cts))
                # print('bkg_counts={0}'.format(bkg_cts))
                # print(freq[100]-freq[99])
                lc = get_hist_withbkg(time, bkg_time, bin_len, time[0], time[-1])
                counts=np.sum(lc.counts)
                cts_rate.append(counts/exptime)
                # flux=get_hist(time,bin_len)
                x=lc.time;flux=lc.counts
                print('k={0}'.format(k))

                T_tot=lc.time[-1]-lc.time[0]
                # freq = get_freq_unsamp(exptime)
                freq=np.arange(1/T_tot,0.5/bin_len-0.0002,1/(5*T_tot))
                freq=freq[np.where(freq > 1 / 20000.)]
                (FP, out_period) = get_LS(x, flux, freq, str(source_id[i]), k)
                FP_print.append(FP);out_period_print.append(out_period);source_name.append(source_id[i])
        result = np.column_stack((source_name,FP_print, out_period_print,src_cts,bkg_cts,cts_rate))
        # np.savetxt('/Users/baotong/Desktop/CDFS/fig_LS_ep{0}_ovsamp_5_baluev/LS_result_{0}_new_psf90.txt'.format(k),result,fmt='%10d %15.10f %15.5f %10d %10d %15.10f')
        # np.savetxt('/Users/baotong/Desktop/CDFS/fig_LS_ep{0}_samp_1_baluev/LS_result_{0}_new_psf90.txt'.format(k),
        #            result, fmt='%10d %15.10f %15.5f %10d %10d %15.10f')

if __name__=='__main__':
    bin_len = 100
    high_id_ep3 = ['643']
    # source_id=np.arange(1,500,1)
    source_id = high_id_ep3

    figurepath = '/Users/baotong/Desktop/aas/AGN_CDFS/figure/'
    plot_CDFS_ep_LS([4])
    # plt.figure(1, (8, 8))
    # plt.subplot(211)
    # plot_CDFS_ep_LS([2])
    # plt.subplot(212)
    # plot_CDFS_ep_LS([3])
    # plt.savefig(figurepath+'XID19.eps',bbox_inches='tight',pad_inches=0.0)
    plt.show()