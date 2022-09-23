#!/bin/bash
# -*- coding: utf-8 -*-
# written by Tong
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
from scipy.optimize import curve_fit
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.timeseries import LombScargle
import hawkeye as hawk

pos_all={'Tuc':[6.0236250, -72.0812833, 3.17 * 60, 3.17 / 8.8 * 60],
         'terzan5':[267.0202083, -24.7790556, 0.72 * 60, 0.72 / 3.4 * 60],
         'M28':[276.1363750, -24.8702972, 1.97 * 60, 1.97 / 8.2 * 60],
         'omg':[201.69700, -47.47947, 5 * 60, 5 / 2.1 * 60],
         'NGC6397':[265.17539, -53.67433, 2.9 * 60, 2.9 / 58 * 60],
         'NGC6752':[287.71713, -59.98455, 1.91* 60, 1.91 / 11.24 * 60],
         'NGC6266':[255.303333,-30.113722,0.92* 60,0.92/4.2*60],
         'M30':[325.092167,-23.179861,1.03* 60,1.03/17.2*60]}

gcname=['Tuc','terzan5','omg','M28','NGC6397','NGC6752','NGC6266','M30']
srcnum=[592,489,300,502,376,244,146,84]
catname=['xray_properties-592.fits','cheng2019_terzan.fit','cheng2020_omg.fit',
         'cheng2020_M28.fit','ngc6397_catalog.fits','ngc6752_catalog.fits',
         'NGC6266_p50_i5_src_1_2_4_8.fits','M30_p50_i5_src_1_2_4_8.fits']

def plot_profile():
    # for i in range(len(gcname)):
    #     path = '/Users/baotong/Desktop/period_' + gcname[i] + '/'
    #     src_info=np.loadtxt(path+'src_info.txt')
    #     counts_all=src_info[:,3];exptime_all=src_info[:,4];VI_all=src_info[:,5]
    #     print('bright_src=',len(np.where(counts_all>100)[0]))
    #     bins_counts=np.logspace(0.5,4,30)
    #     plt.hist(counts_all,bins=bins_counts,histtype='step',lw=2,linestyle='-')
    #     plt.semilogx()
    #     plt.show()
    #     plt.scatter(exptime_all,VI_all)
    #     plt.semilogy()
    #     plt.show()
    label=['x','^','v','o','D','*','s']
    color_list=['r','g','b','k','orange','purple','magenta']
    result_all=pd.read_excel('/Users/baotong/Desktop/period_terzan5/candidate_allGC.xlsx','allbutTuc')
    ra = np.array(result_all['ra'])
    dec = np.array(result_all['dec'])
    seq = np.array(result_all['seq'])
    period=np.array(result_all['period_all'])
    type = np.array(result_all['GC'])
    dist=np.array(result_all['proj_dist'])
    counts=np.array(result_all['counts'])
    exptime=np.array(result_all['expT'])
    CR=counts/exptime
    # plt.scatter(period/3600,CR)
    # plt.loglog()

    bins_p=np.logspace(np.log10(1.0), np.log10(15), 10)
    P_min = 7./6.
    P_gap = [7740.0 / 3600., 11448.0 / 3600.]
    fig = plt.figure(1, figsize=(9, 6))
    ax1 = fig.add_subplot(111)
    ax1.plot()
    ax1.hist(period/3600, bins=bins_p, histtype='step', lw=1.5, color='blue', linestyle='-')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.plot([P_min, P_min], [0, 10], '--',color='grey')
    ax1.plot([P_gap[0], P_gap[0]], [0, 10], '-',lw=2.,color='orange')
    ax1.plot([P_gap[1], P_gap[1]], [0, 10], '-',lw=2.,color='orange')
    ax1.text(P_gap[0]+0.3, 10, 'gap',fontsize=14)
    ax1.text(P_min-0.2, 10, 'minimum',fontsize=14)
    plt.show()

    # for i in range(len(gcname)-1):
    #     gc_srcid=np.where(type==gcname[i])[0]
    #     plt.scatter(period[gc_srcid]/3600,dist[gc_srcid]/pos_all[gcname[i]][3],marker=label[i],color=color_list[i],label=gcname[i],s=50)
    # plt.plot([P_min, P_min], [0, 10], '--',color='grey')
    # plt.plot([P_gap[0], P_gap[0]], [0, 10], '-',lw=2.,color='orange')
    # plt.plot([P_gap[1], P_gap[1]], [0, 10], '-',lw=2.,color='orange')
    # plt.xlabel('Period (h)',hawk.font1)
    # plt.ylabel(r'R/$r_{c}$',hawk.font1)
    # plt.tick_params(labelsize=18)
    # plt.loglog()
    # plt.legend()
    # plt.show()
    return None
def plot_P_L(save=0,show=1):
    label=['x','^','v','o','D','*','s']
    color_list=['r','g','b','k','orange','purple','magenta']
    result_all=pd.read_excel('/Users/baotong/Desktop/period_terzan5/candidate_allGC.xlsx','all')
    idex=np.where(result_all['judge']=='CV')[0]
    ra = np.array(result_all['ra'])[idex]
    dec = np.array(result_all['dec'])[idex]
    seq = np.array(result_all['seq'])[idex]
    period=np.array(result_all['period_all'])[idex]
    type = np.array(result_all['GC'])[idex]
    dist=np.array(result_all['proj_dist'])[idex]
    counts=np.array(result_all['counts'])[idex]
    exptime=np.array(result_all['expT'])[idex]
    L=np.array(result_all['L'])[idex]
    P_min = 7./6.
    P_gap = [7740.0 / 3600., 11448.0 / 3600.]
    print(len(L))
    fig, (ax1, ax2) = plt.subplots(ncols=1, nrows=2, sharex='all',gridspec_kw={'height_ratios': [1, 1]}, figsize=(10, 15))
    # fig=plt.figure(1,figsize=(10,20))
    # ax1 = fig.add_subplot(211)
    # ax2 = fig.add_subplot(212)
    for i in range(len(gcname)-1):
        gc_srcid=np.where(type==gcname[i])[0]
        ax1.scatter(period[gc_srcid]/3600,L[gc_srcid],marker=label[i],color=color_list[i],label=gcname[i],s=50)
        ax2.scatter(period[gc_srcid] / 3600, dist[gc_srcid] / pos_all[gcname[i]][2], marker=label[i],
                    color=color_list[i], label=gcname[i], s=50)
    ax1.plot([P_min, P_min], [0, 1e33], '--',color='grey')
    ax1.plot([P_gap[0], P_gap[0]], [0, 1e33], '-',lw=2.,color='orange')
    ax1.plot([P_gap[1], P_gap[1]], [0, 1e33], '-',lw=2.,color='orange')
    ax1.loglog()
    ax1.set_ylabel(r'Luminosity ($\rm erg~s^{-1}$)',hawk.font1)
    # ax1.set_xlabel('Period (h)',hawk.font1)
    ax1.set_xlim(1,28)
    ax1.tick_params(labelsize=18)
    # ax1.legend()
    ax1.sharex(ax2)
    ax2.plot([P_min, P_min], [0, 10], '--',color='grey')
    ax2.plot([P_gap[0], P_gap[0]], [0, 10], '-',lw=2.,color='orange')
    ax2.plot([P_gap[1], P_gap[1]], [0, 10], '-',lw=2.,color='orange')
    ax2.set_xlabel('Period (h)',hawk.font1)
    ax2.set_ylabel(r'R/$r_{h}$',hawk.font1)
    ax2.set_xlim(1,28)
    ax2.set_ylim(1e-2,20)
    ax2.tick_params(labelsize=18)
    ax2.loglog()
    ax2.legend(loc='upper center',ncol=7,handletextpad=0.1,columnspacing=0.1,handlelength=1.5,edgecolor='grey')
    plt.subplots_adjust(wspace=0,hspace = 0.05)

    # plt.figure(2)
    # plt.hist(np.log10(L),bins=np.linspace(31,33.5,15),histtype='step')
    # plt.semilogx()
    # plt.show()
    path_out='/Users/baotong/Desktop/aas/GCall/figure/'
    if save:
        plt.savefig(path_out + 'GC_CV_profile.pdf', bbox_inches='tight', pad_inches=0.05)
    if show:
        plt.show()

def plot_spec():
    label=['x','^','v','o','D','*','s']
    color_list=['r','g','b','k','orange','purple','magenta']
    result_all=pd.read_excel('/Users/baotong/Desktop/period_terzan5/candidate_allGC.xlsx','allbutTuc')
    ra = np.array(result_all['ra'])
    dec = np.array(result_all['dec'])
    seq = np.array(result_all['seq'])
    period=np.array(result_all['period_all'])
    type = np.array(result_all['GC'])
    dist=np.array(result_all['proj_dist'])
    counts=np.array(result_all['counts'])
    exptime=np.array(result_all['expT'])
    L=np.array(result_all['L'])
    kT=np.array(result_all['kT'])
    plt.figure(1)
    plt.scatter(kT,L)
    plt.loglog()
    plt.show()
if __name__=='__main__':
    plot_P_L(save=1,show=1)
    # plot_spec()