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
import NGC104.plot_pXS as plot_pXS
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
label = ['.', '^', 'v', 'o', 'D', '*', 's']
color_list = ['grey', 'g', 'b', 'k', 'orange', 'purple', 'magenta']
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

def plot_Phist(save=0,show=1):
    result_all=pd.read_excel('/Users/baotong/Desktop/period_terzan5/candidate_allGC.xlsx','all')
    result_allbutTuc=pd.read_excel('/Users/baotong/Desktop/period_terzan5/candidate_allGC.xlsx','allbutTuc')
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

    idex2=np.where(result_allbutTuc['judge']=='CV')[0]
    period_allbutTuc=np.array(result_allbutTuc['period_all'])[idex2]
    period_allbutTuc/=3600.
    period_GC = period / 3600.

    P_min = 7./6.
    P_gap = [7740.0 / 3600., 11448.0 / 3600.]

    bins = np.logspace(np.log10(0.5), np.log10(30), 41)
    bins_2 = np.logspace(np.log10(0.8), np.log10(100), 41)
    bins_spin = np.logspace(np.log10(3 / 36.), np.log10(2), 31)
    bins_p = np.logspace(np.log10(0.5), np.log10(12), 16)
    bins_p=np.concatenate((bins_p,[15,20.,25.,30.]))
    path_fits = '/Users/baotong/Desktop/period_LW/'
    RK = fits.open(path_fits + 'RK14.fit')
    orb = RK[1].data['Orb_Per']
    type1 = RK[1].data['Type1'];
    type2 = RK[1].data['Type2'];
    type3 = RK[1].data['Type3']
    M1 = RK[1].data['M1'];
    M2 = RK[1].data['M2'];
    M1_M2 = RK[1].data['M1_M2']
    orb = orb * 24;
    spin = RK[1].data['_3___Per']
    orb_DN = orb[np.where(type1 == 'DN')]
    orb_Polar = orb[np.where((type2 == 'AM') | (type3 == 'AM'))]
    orb_IP = orb[np.union1d(np.where((type2 == 'IP') | (type2 == 'DQ')), np.where((type3 == 'IP') | (type3 == 'DQ')))]
    spin_IP = spin[np.union1d(np.where((type2 == 'IP') | (type2 == 'DQ')), np.where((type3 == 'IP') | (type3 == 'DQ')))]
    spin_IP /= 3600.
    orb_all = np.concatenate((orb_DN, orb_IP, orb_Polar))
    orb_all = orb_all[np.where(orb_all > 0)]
    print('Too long period CV', len(orb_all[np.where(orb_all > 12.)]))
    print('CV number', len(orb_all))
    print('short period CV', len(orb_all[np.where(orb_all < 7700 / 3600.)]))
    print('long period CV', len(orb_all[np.where(orb_all > 11448.0 / 3600)]))
    (ra, dec, seq, period, L, Lmin, Lmax, type) =plot_pXS.load_LW('result_LW')
    period_LW = period[np.where((period > 3900) & (period < 40000))]
    period_LW /= 3600.
    # fig, axes = plt.subplots(2, 1, figsize=(15, 10), sharex='all')
    fig, (ax1, ax2) = plt.subplots(ncols=1, nrows=2, sharex='all', gridspec_kw={'height_ratios': [2, 1]},
                                   figsize=(9, 6))
    # ax1=fig.add_subplot(211)
    ax1.plot()
    ax1.hist(orb_all, bins=bins, histtype='step', lw=2, color='red', linestyle='--')
    ax1.hist(period_allbutTuc, bins=bins_p, histtype='step', lw=2, linestyle='-', facecolor='grey',
             hatch='/', edgecolor='k', fill=True)
    ax1.hist(period_GC, bins=bins_p, histtype='step', lw=3, color='c', linestyle='-')
    # ax1.hist(spin_IP, bins = bins_spin, histtype = 'step',lw=1.5, color = 'purple',linestyle='-')
    ax1.legend(['CVs in Solar Neighborhood', 'CVs in GCs(excluding 47 Tuc)', 'CVs in GCs'])
    P_min = 1.373333333
    P_gap = [7740.0 / 3600., 11448.0 / 3600.]
    ax1.set_ylim(8e-1, 360)
    ax1.set_xlim(0.5, 60)
    ax1.plot([P_min, P_min], [0, 220], '--',color='grey')
    ax1.plot([P_gap[0], P_gap[0]], [0, 220], '-',lw=2.,color='orange')
    ax1.plot([P_gap[1], P_gap[1]], [0, 220], '-',lw=2.,color='orange')
    ax1.text(P_gap[0]+0.3, 250, 'gap',fontsize=14)
    ax1.text(P_min-0.2, 250, 'minimum',fontsize=14)
    # ax1.set_xlabel('Period (hours)',font1)
    ax1.set_ylabel('Number of sources',plot_pXS.font1)
    ax1.tick_params(labelsize=16)
    ax1.set_xscale('log')
    ax1.set_yscale('log')

    # ax2=fig.add_subplot(212)
    ax2.hist(orb_all,bins=bins,histtype='step',lw=2,color='red',cumulative=1,density=1,linestyle='--')
    ax2.hist(period_allbutTuc, bins=bins_p, histtype='step', lw=2, linestyle='-',cumulative=1,density=1,facecolor='grey',
         hatch='/', edgecolor='k',fill=True)
    ax2.hist(period_GC, bins=bins_p, histtype='step', lw=3, color='c',cumulative=1,density=1, linestyle='-')
    # ax2.hist(spin_IP, bins = bins_spin, histtype = 'step',lw=1.5, color = 'purple',cumulative=1,density=1,linestyle='-')

    ax2.plot([P_min, P_min], [0, 1], '--',color='grey')
    ax2.plot([P_gap[0], P_gap[0]], [0, 1], '-',lw=2.,color='orange')
    ax2.plot([P_gap[1], P_gap[1]], [0, 1], '-',lw=2.,color='orange')

    ax2.set_xscale('log')
    ax2.set_xlabel('Period (hours)',plot_pXS.font1)
    ax2.set_ylabel('CDF',plot_pXS.font1)
    ax2.tick_params(labelsize=16)
    # ax2.set_yscale('log')
    path_out = '/Users/baotong/Desktop/aas/GCall/figure/'
    if save:
        plt.savefig(path_out + 'GC_NP.pdf', bbox_inches='tight', pad_inches=0.05)
    if show:
        plt.show()

def plot_P_L_profile(save=0,show=1):
    label=['x','^','v','o','D','*','s']
    color_list=['r','g','b','k','orange','purple','magenta']
    result_all=pd.read_excel('/Users/baotong/Desktop/period_terzan5/candidate_allGC.xlsx','allbutTuc')
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
        plt.savefig(path_out + 'GC_CV_profile_butTuc.pdf', bbox_inches='tight', pad_inches=0.05)
    if show:
        plt.show()
    return period
def plot_P_L(save=0,show=1):
    result_all = pd.read_excel('/Users/baotong/Desktop/period_terzan5/candidate_allGC.xlsx', 'all')
    idex = np.where(result_all['judge'] == 'CV')[0]
    ra = np.array(result_all['ra'])[idex]
    dec = np.array(result_all['dec'])[idex]
    seq = np.array(result_all['seq'])[idex]
    period = np.array(result_all['period_all'])[idex]
    type = np.array(result_all['GC'])[idex]
    dist = np.array(result_all['proj_dist'])[idex]
    counts = np.array(result_all['counts'])[idex]
    exptime = np.array(result_all['expT'])[idex]
    L = np.array(result_all['L'])[idex]
    P_min = 7. / 6.
    P_gap = [7740.0 / 3600., 11448.0 / 3600.]
    print(len(L))

    fig, ax1 = plt.subplots(ncols=1, nrows=1,figsize=(10, 10))
    # ax1 = fig.add_subplot(211)
    # ax2 = fig.add_subplot(212)
    for i in range(len(gcname) - 1):
        gc_srcid = np.where(type == gcname[i])[0]
        if len(gc_srcid)==0:continue
        else:
            if i == 0:
                s = 80
            else:
                s = 120
            ax1.scatter(period[gc_srcid] / 3600, L[gc_srcid], marker=label[i], color=color_list[i], label=gcname[i],
                    s=s)
    ax1.plot([P_min, P_min], [0, 1e33], '--', color='grey')
    ax1.plot([P_gap[0], P_gap[0]], [0, 1e33], '-', lw=2., color='orange')
    ax1.plot([P_gap[1], P_gap[1]], [0, 1e33], '-', lw=2., color='orange')
    ax1.loglog()
    # ax1.set_yscale('log')
    ax1.set_ylabel(r'Luminosity ($\rm erg~s^{-1}$)', hawk.font1)
    # ax1.set_xlabel('Period (h)',hawk.font1)
    ax1.set_xlim(0.8, 14)
    ax1.set_ylim(5e30, 1e34)
    ax1.tick_params(labelsize=18)
    ax1.set_xlabel('Period (h)', hawk.font1)
    ax1.tick_params(labelsize=18)
    ax1.legend(loc='lower center', ncol=len(gcname), handletextpad=0.1, columnspacing=0.1, handlelength=1.5, edgecolor='grey')
    plt.subplots_adjust(wspace=0, hspace=0.05)
    # plt.figure(2)
    # plt.hist(np.log10(L),bins=np.linspace(31,33.5,15),histtype='step')
    # plt.semilogx()
    # plt.show()
    path_out = '/Users/baotong/Desktop/aas/GCall/figure/'
    if save:
        plt.savefig(path_out + 'GC_CV_P_L.pdf', bbox_inches='tight', pad_inches=0.05)
    if show:
            plt.show()
    return period
def plot_P_L_type(save=0,show=1):
    label = ['x', '^', 'v']
    color_list = ['r', 'g', 'b', 'k', 'orange', 'purple', 'magenta']
    result_all = pd.read_excel('/Users/baotong/Desktop/period_terzan5/candidate_allGC.xlsx', 'all')
    idex = np.where(result_all['judge'] == 'CV')[0]
    ra = np.array(result_all['ra'])[idex]
    dec = np.array(result_all['dec'])[idex]
    seq = np.array(result_all['seq'])[idex]
    period = np.array(result_all['period_all'])[idex]
    type = np.array(result_all['GC'])[idex]
    dist = np.array(result_all['proj_dist'])[idex]
    counts = np.array(result_all['counts'])[idex]
    exptime = np.array(result_all['expT'])[idex]
    L = np.array(result_all['L'])[idex]
    P_min = 7. / 6.
    P_gap = [7740.0 / 3600., 11448.0 / 3600.]
    print(len(L))

    fig, ax1 = plt.subplots(ncols=1, nrows=1, figsize=(10, 6))
    # ax1 = fig.add_subplot(211)
    # ax2 = fig.add_subplot(212)
    for i in range(len(gcname) - 1):
        gc_srcid = np.where(type == gcname[i])[0]
        if len(gc_srcid) == 0:
            continue
        else:
            ax1.scatter(period[gc_srcid] / 3600, L[gc_srcid], marker=label[i], color=color_list[i], label=gcname[i],
                        s=50)
    ax1.plot([P_min, P_min], [0, 1e33], '--', color='grey')
    ax1.plot([P_gap[0], P_gap[0]], [0, 1e33], '-', lw=2., color='orange')
    ax1.plot([P_gap[1], P_gap[1]], [0, 1e33], '-', lw=2., color='orange')
    ax1.loglog()
    ax1.set_ylabel(r'Luminosity ($\rm erg~s^{-1}$)', hawk.font1)
    # ax1.set_xlabel('Period (h)',hawk.font1)
    ax1.set_xlim(1, 15)
    ax1.tick_params(labelsize=18)
    ax1.set_xlabel('Period (h)', hawk.font1)
    ax1.tick_params(labelsize=18)
    ax1.legend(loc='upper center', ncol=6, handletextpad=0.1, columnspacing=0.1, handlelength=1.5, edgecolor='grey')
    plt.subplots_adjust(wspace=0, hspace=0.05)
    # plt.figure(2)
    # plt.hist(np.log10(L),bins=np.linspace(31,33.5,15),histtype='step')
    # plt.semilogx()
    # plt.show()
    path_out = '/Users/baotong/Desktop/aas/GCall/figure/'
    if save:
        plt.savefig(path_out + 'GC_CV_P_L.pdf', bbox_inches='tight', pad_inches=0.05)
    if show:
        plt.show()
    return period

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


def plot_src_info(save=0,show=1):
    label = ['x', '^', 'v', 'o', 'D', '*', 's','.']
    color_list = ['r', 'g', 'b', 'k', 'orange', 'purple', 'magenta','cyan']
    for i in range(len(gcname)):
        path = '/Users/baotong/Desktop/period_' + gcname[i] + '/'
        src_info=np.loadtxt(path+'src_info.txt')
        counts_all=src_info[:,3];exptime_all=src_info[:,4];VI_all=src_info[:,5]
        print('bright_src=',len(np.where(counts_all>100)[0]))
        bins_counts=np.logspace(0.5,4,20)
        plt.hist(counts_all,bins=bins_counts,histtype='step',lw=2,linestyle='-',color=color_list[i],label=gcname[i])
        plt.semilogx()
        # plt.show()
        # plt.scatter(exptime_all,VI_all)
        # plt.semilogy()
        # plt.show()
    plt.legend()
    plt.show()
if __name__=='__main__':
    plot_P_L(save=0,show=1)
    # plot_spec()
    # plot_Phist(save=1,show=1)
    # plot_src_info(show=1)