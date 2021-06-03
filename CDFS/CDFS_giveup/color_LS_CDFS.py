#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib import ticker, cm
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
def get_hist(t, len_bin):
    ###将输入的time信息，按照len_bin的长度输出为lc
    t_test = t;dt=len_bin
    ev = EventList()
    ev.time = t_test
    lc_new = ev.to_lc(dt=dt, tstart=ev.time[0] - 0.5 * dt, tseg=ev.time[-1] - ev.time[0])
    return lc_new
def get_hist_withbkg(t,t_bkg, len_bin):
    ###将输入的time信息，按照len_bin的长度输出为lc
    t_test = t;t_bkg_test=t_bkg;dt=len_bin
    t_bkg_test=np.delete(t_bkg_test,t_bkg_test<0)
    ev = EventList();ev_bkg=EventList()
    ev.time=t_test;ev_bkg.time=t_bkg_test
    lc_new = ev.to_lc(dt=dt, tstart=ev.time[0] - 0.5 * dt, tseg=ev.time[-1] - ev.time[0])
    lc_bkg = ev_bkg.to_lc(dt=dt, tstart=ev.time[0] - 0.5 * dt, tseg=ev.time[-1] - ev.time[0])
    lc_out=lc_new
    lc_out.counts=lc_new.counts-(1/12.)*lc_bkg.counts
    print(len(lc_out.time))
    return lc_out
# def get_LS_myself(time,flux,freq):

def get_LS(time, flux,freq,dataname):
    x = time
    y = flux
    # dy=np.sqrt(y)
    # plt.scatter(x,y)
    # plt.show()

    # LS = LombScargle(x, y, dy = 1, normalization = 'standard', fit_mean = True,
    #                  center_data = True).power(freq, method = 'cython')
    LS = LombScargle(x, y,normalization = 'standard')
    # LS = LombScargle(x, y, normalization='psd')
    power = LS.power(freq)
    power_freq=np.max(power[np.where((1/freq-2493.76)<20)])
    FP=LS.false_alarm_probability(power_freq,minimum_frequency = freq[0],maximum_frequency = freq[-1],method='baluev')
    # FP_99 = LS.false_alarm_level(0.0027,minimum_frequency = freq[0], maximum_frequency = freq[-1],method='baluev')
    # FP_90 = LS.false_alarm_level(0.05, minimum_frequency=freq[0],
    #                              maximum_frequency=freq[-1], method='baluev')
    # FP_68 = LS.false_alarm_level(0.32,minimum_frequency=freq[0],
    #                              maximum_frequency=freq[-1], method='baluev')
    #
    # if FP<0.01:print(dataname)
    plt.title('{0},FP={1}'.format(dataname,FP))
    # plt.semilogx()
    plt.plot(freq, power)
    # plt.plot([1/2492,1/2492.],[0,np.max(power)],'--',linewidth=1)
    plt.semilogx()
    print(1. / freq[np.where(power == np.max(power))])
    # if FP<0.01:print(1./freq[np.where(power==np.max(power))]);print(np.where(power==np.max(power)))
    out_period=1./freq[np.where(power==np.max(power))]
    # if np.abs(out_period-950.7317)>30:FP=1
    # plt.savefig('/Users/baotong/Desktop/CDFS/fig_LS_ep{0}_ovsamp_5_baluev/{1}.eps'.format(k,dataname))
    plt.close()
    return [FP,out_period]
def LSP_band(source_id,k_num,band,epochnum=[0,31]):
    bin_len=100
    # k_num=[1,2,3,4]
        ##写入txt文本方便plot
    path = '/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep{0}/'.format(k_num)
    srcevt=np.loadtxt(path+'{0}.txt'.format(source_id))
    bkgevt=np.loadtxt(path+'{0}_bkg.txt'.format(source_id))
    epoch = np.loadtxt(path + 'epoch_src_{0}.txt'.format(source_id))
    useid=epoch[:,2][epochnum[0]:epochnum[1]]
    if len(epoch)==0:
        return (1,0)
    if len(np.shape(epoch)) == 1:
        exptime = epoch[3]
    else:
        exptime = np.sum(epoch[:, 3])
    if len(srcevt)<4:
        return (1,0)
    if len(bkgevt)<2:bkgevt=[]

    else:
        time=[];bkg_time=[];energy=[];bkg_energy=[]
        srcevt[:, -1] = srcevt[:, -1].astype('int');
        bkgevt[:, -1] = bkgevt[:, -1].astype('int')
        for k in range(len(useid)):
            time =np.concatenate((time, srcevt[:, 0][np.where(srcevt[:, -1] == useid[k])]))
            energy=np.concatenate((energy, srcevt[:, 1][np.where(srcevt[:, -1] == useid[k])]))
            bkg_time = np.concatenate((bkg_time,bkgevt[:, 0][np.where(bkgevt[:, -1] == useid[k])]))
            bkg_energy=np.concatenate((bkg_energy,bkgevt[:, 1][np.where(bkgevt[:, -1] == useid[k])]))
        ## 看一下不同波段的lc如何 ##
        (time,energy)=filter_energy(time,energy,band)
        (bkg_time,bkg_energy)=filter_energy(bkg_time,bkg_energy,band)
        src_cts=len(time);bkg_cts=len(bkg_time)
        if len(time) < 4:
            return (1, 0)
        if len(bkg_time) < 2: bkg_time = []

        T_tot=time[-1]-time[0]
        freq=np.arange(1/T_tot,0.5/bin_len,1/(5*T_tot))
        freq=freq[np.where(freq > 1 / 20000.)]
        lc=get_hist_withbkg(time,bkg_time,bin_len)
        counts=np.sum(lc.counts)
        cts_rate=(counts/exptime)
        x=lc.time;flux=lc.counts
        (FP, out_period) = get_LS(x, flux, freq, source_id)
    return (FP,out_period)
# LSP_band('19','3',[500,8000])

def plot_color_band_LSP(source_id,k_num,binsize=5):
    # FP=np.zeros((binsize,binsize))
    # out_period=np.zeros((binsize,binsize))
    # band0=np.linspace(500,8000,binsize)
    # band1=np.linspace(500,8000,binsize)
    # for i in range(binsize):
    #     for j in range(binsize):
    #         (FP[i][j],out_period[i][j])=LSP_band(source_id,k_num,[band0[i],band1[j]])
    # print(FP)
    # im=plt.contourf(band0,band1,FP.T,locator=ticker.LogLocator(),levels=np.logspace(-3,0,10),cmap="Oranges_r")
    # plt.colorbar(im)
    # plt.show()

    epoch_id_num = 12
    # FP=np.zeros((epoch_id_num,epoch_id_num))+1
    # out_period=np.zeros((epoch_id_num,epoch_id_num))
    # for i in range(epoch_id_num):
    #     for j in range(i,epoch_id_num):
    #         (FP[i][j],out_period[i][j])=LSP_band(source_id,k_num,[500,8000],epochnum=[i,j])
    # np.savetxt('/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep3/FP_obs_color_{0}.txt'.format(source_id),FP)
    # np.savetxt('/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep3/period_obs_color_{0}.txt'.format(source_id),out_period)

    FP=np.loadtxt('/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep3/FP_obs_color_{0}.txt'.format(source_id))

    im=plt.contourf(np.linspace(1,epoch_id_num,epoch_id_num),np.linspace(1,epoch_id_num,epoch_id_num),FP.T,levels=np.linspace(0,1,11),cmap="OrRd_r")
    plt.colorbar(im)
    plt.show()
if __name__=='__main__':
    plot_color_band_LSP('19','3')