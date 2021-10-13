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
import stingray as sr
from stingray.events import EventList
from stingray.lightcurve import Lightcurve
from stingray import Lightcurve, Crossspectrum, sampledata,Powerspectrum,AveragedPowerspectrum
from stingray.simulator import simulator, models
from astropy.stats import poisson_conf_interval

def plot_LS_all_single_obs():
    bin_len=100
    # obsid='700011'
    obsIDlist=[700011,700013,700014,700163,700173,700174,700175]
    for obsid in obsIDlist:
        obsid=str(obsid)
    # path='/Users/baotong/eSASS/data/raw_data/47_Tuc/txt/txt_merge_0.5_5/'
        path='/Users/baotong/eSASS/data/raw_data/47_Tuc/txt/txt_psf75_{0}/'.format(obsid)
        figurepath='/Users/baotong/eSASS/data/raw_data/47_Tuc/fig_LS/LS_psf75_{0}_bin{1}/'.format(obsid,bin_len)
        if not os.path.exists(figurepath):os.mkdir(figurepath)
        srcID=np.arange(1,889,1)
        FP_all=[];out_period_all=[];srcID_output=[];counts_all=[]
        for i in range(len(srcID)):
            srcevt=np.loadtxt(path+'{0}_{1}.txt'.format(srcID[i],obsid))
            if len(srcevt)==0:continue
            if srcevt.ndim==1:continue
            if len(srcevt)<10:continue ##去掉10个光子以下的
            counts=len(srcevt)
            time=srcevt[:,0];energy=srcevt[:,1]
            # (TSTART, TSTOP, OBSID, exptime) = funcs.read_epoch(path + 'epoch_src_{0}.txt'.format(srcID[i]))
            # lc=funcs.get_hist(time,len_bin=bin_len,tstart=TSTART[0], tstop=TSTOP[-1])
            lc=funcs.get_hist(time,len_bin=bin_len)
            x = lc.time;flux = lc.counts

            print(i)
            T_tot = lc.time[-1] - lc.time[0]
            freq = np.arange(1 / T_tot, 0.5 / bin_len, 1 / (5 * T_tot))
            freq=freq[np.where(freq > 1 / 10000.)]
            [FP, out_period,max_NormLSP] = funcs.get_LS(x, flux, freq, outpath=figurepath,outname=str(srcID[i]),show=False)
            FP_all.append(FP);out_period_all.append(out_period);srcID_output.append(srcID[i]);counts_all.append(counts)

        FP_all=np.array(FP_all);out_period_all=np.array(out_period_all);srcID_output=np.array(srcID_output);counts_all=np.array(counts_all)
        LS_info=np.column_stack((srcID_output,FP_all,out_period_all,counts_all))
        np.savetxt(figurepath+'LS_info_bin{0}.txt'.format(bin_len),LS_info,fmt='%10d %10.5f %10.5f %10d')

def filter_obs(src_evt,useid):
    src_evt_use = src_evt[np.where(src_evt[:-1] == useid[0])[0]]
    i=1
    while i < len(useid):
        id=useid[i]
        src_evt_use_temp=src_evt[np.where(src_evt[:-1]==id)[0]]
        src_evt_use = np.concatenate((src_evt_use, src_evt_use_temp))
        i+=1
    return src_evt_use

def get_netlc_merge(srcID,obsIDlist,ecf=75,bin_len=100):
    record=0
    src_cts=0;bkg_cts=0;net_cts=0;bkgscale_meanlist=[]
    for i in range(len(obsIDlist)):
        obsid=obsIDlist[i]
        path='/Users/baotong/eSASS/data/raw_data/47_Tuc/txt/txt_psf{0}_{1}/'.format(ecf,obsid)
        srcevt = np.loadtxt(path + '{0}_{1}.txt'.format(srcID,obsid))
        bkgevt= np.loadtxt(path + '{0}_bkg_{1}.txt'.format(srcID,obsid))
        if srcevt.ndim<2 or bkgevt.ndim<2 or len(srcevt)<10:
            continue
        src_area=np.loadtxt(path+'src_area.txt')
        bkg_area=np.loadtxt(path+'bkg_area.txt')
        bkgscale=(src_area[:,1][np.where(src_area[:,0]==srcID)])/(bkg_area[:,1][np.where(bkg_area[:,0]==srcID)])
        bkgscale=bkgscale[0]
        time = srcevt[:, 0];
        energy_src = srcevt[:, 1]
        time_bkg = bkgevt[:, 0];
        energy_bkg = bkgevt[:, 1]

        lc = funcs.get_hist_withbkg(time,time_bkg, len_bin=bin_len,bkgscale=bkgscale,tstart=0,tstop=0)
        src_cts+=len(time);bkg_cts+=len(time_bkg);net_cts+=np.sum(lc.counts)
        bkgscale_meanlist.append(bkgscale)
        if i==0:
            lc_all=lc;
            record=1
        elif record==0:
            lc_all=lc
        else:
            lc_all=lc_all.join(lc)
    return (lc_all,bkgscale_meanlist,src_cts,bkg_cts,net_cts)

def plot_LS_merge(srcID,obsIDlist,ecf=75,bin_len=100):
    (lc_all,bkgscale_meanlist,src_cts,bkg_cts,net_cts)=get_netlc_merge(srcID, obsIDlist, bin_len=bin_len)
    bkgscale_mean = np.mean(np.array(bkgscale_meanlist))
    x=lc_all.time;flux=lc_all.counts
    T_tot = x[-1] - x[0]
    freq = np.arange(1 / T_tot, 0.5 / bin_len, 1 / (10* T_tot))
    freq = freq[np.where(freq > 1 / 10000.)]
    figurepath='/Users/baotong/eSASS/data/raw_data/47_Tuc/fig_LS/LS_net_psf{0}_bin{1}/'.format(ecf,bin_len)
    [FP, out_period, max_NormLSP] = funcs.get_LS(x, flux, freq, outpath=figurepath, outname=str(srcID),save=True,show=False)
    return [srcID,1-FP,out_period,bkgscale_mean,src_cts,bkg_cts,net_cts]


def pfold_fromlc(lc,epoch_info,p_test,bin,shift,path_out,label='test'):
    x=lc.time
    flux=lc.counts
    plt.step(x,flux)
    plt.show()
    net_counts=np.sum(flux)
    def trans(t,p_test,shift):
        ti =t
        v = 1.0 /p_test
        turns = v * ti
        turns += shift
        # 初始相位
        for i in range(len(turns)):
            turns[i] = turns[i] - int(turns[i])
        return turns

    turns=trans(x,p_test,shift)
    loc=np.zeros(bin)
    for i in range(len(turns)):
        loc[int(turns[i]*bin)] += flux[i]
    x = np.array([(i / bin + 0.5 / bin) for i in range(bin)])
    x2=np.concatenate((x,x+1))
    y2=np.concatenate((loc,loc))
    T_in_perbin = funcs.get_T_in_mbins(epoch_info, 2 * np.pi / p_test, bin, shift * 2 * np.pi)
    correct_gap = T_in_perbin / (sum(T_in_perbin) / len(T_in_perbin))
    y2 /= np.concatenate((correct_gap, correct_gap))
    y2_err=np.array(poisson_conf_interval(y2,interval='frequentist-confidence'))
    y2_err[0]=y2-y2_err[0]
    y2_err[1]=y2_err[1]-y2

    plt.title("#{0} P={1:.2f},net_C={2:.2f}".format(label,p_test,net_counts), fontsize = 18)
    plt.xlabel('Phase',funcs.font1)
    plt.ylabel('Counts/bin',funcs.font1)
    plt.tick_params(labelsize = 18)
    plt.ylim(0,(np.max(y2)+np.max(y2)**0.5)*1.05)
    plt.step(x2,y2,color='red')
    plt.errorbar(x2 - 0.5 / bin, y2, yerr=y2_err, fmt='.', capsize=1, elinewidth=1, ecolor='red')
    plt.savefig(path_out + 'pfold_lc_{0}.eps'.format(label))
    plt.show()


if __name__=='__main__':
    figurepath = '/Users/baotong/eSASS/data/raw_data/47_Tuc/fig_LS/LS_net_psf75_bin100/'
    obsIDlist=[700011, 700163, 700013, 700014, 700173, 700174, 700175]
    # obsIDlist = [700011]
    # plot_LS_merge(srcID=5,obsIDlist=obsIDlist)

    # srcIDlist=np.arange(1,889,1)
    # src_LS_info=[]
    # for i in range(len(srcIDlist)):
    #     info_temp=plot_LS_merge(srcID=srcIDlist[i],obsIDlist=obsIDlist)
    #     src_LS_info.append(info_temp)
    # src_LS_info=np.array(src_LS_info)
    # np.savetxt(figurepath+'LS_info.txt',src_LS_info,fmt='%10d %10.5f %10.5f %10.5f %10d %10d %10.2f')


    ##输出net_counts，运行一次即可;更新后已废
    # net_counts=[]
    # for srcID in srcIDlist:
    #     (lc_all,bkgscale_meanlist,src_cts,bkg_cts) = get_netlc_merge(srcID, obsIDlist, bin_len=100)
    #     net_counts.append(np.sum(lc_all.counts))
    # np.savetxt(figurepath+'net_counts.txt',np.array(net_counts),fmt='%10d')

    ##pfold
    srcID=735;
    # [srcID,FP,out_period,bkgscale_mean,src_cts,bkg_cts,net_cts]=plot_LS_merge(srcID, obsIDlist, bin_len=100)
    # print(net_cts)

    (lc_all,bkgscale_meanlist,src_cts,bkg_cts,net_cts) = get_netlc_merge(srcID, obsIDlist, ecf=75,bin_len=100)
    # # print(np.sum(lc_all))
    epoch_info = np.loadtxt('/Users/baotong/eSASS/data/raw_data/47_Tuc/txt/txt_merge_psf75_0.2_5/epoch_src_{0}.txt'.format(srcID))
    # # #
    pfold_fromlc(lc=lc_all,epoch_info=epoch_info,p_test=19872.06265,bin=20,shift=0.,path_out=figurepath,label=str(srcID))