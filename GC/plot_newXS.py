#!/bin/bash
# -*- coding: utf-8 -*-
# written by Tong
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
from matplotlib import ticker
from astropy.io import fits
from astropy.stats import poisson_conf_interval
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy import units as u
from scipy import optimize as op
import os
import hawkeye as hawk
import NGC104.plot_pXS as plot_pXS
import GC.localCVs as localCVs
import NGC104.CV_model as CV_model

pos_all = {'Tuc': [6.0236250, -72.0812833, 3.17 * 60, 3.17 / 8.8 * 60],
           'terzan5': [267.0202083, -24.7790556, 0.72 * 60, 0.72 / 3.4 * 60],
           'M28': [276.1363750, -24.8702972, 1.97 * 60, 1.97 / 8.2 * 60],
           'omg': [201.69700, -47.47947, 5 * 60, 5 / 2.1 * 60],
           'NGC6397': [265.17539, -53.67433, 2.9 * 60, 2.9 / 58 * 60],
           'NGC6752': [287.71713, -59.98455, 1.91 * 60, 1.91 / 11.24 * 60],
           'NGC6266': [255.303333, -30.113722, 0.92 * 60, 0.92 / 4.2 * 60],
           'M30': [325.092167, -23.179861, 1.03 * 60, 1.03 / 17.2 * 60]}

gcname = ['Tuc', 'terzan5', 'M28', 'omg', 'NGC6397', 'NGC6752', 'NGC6266']
catname = ['xray_properties-592.fits', 'cheng2019_terzan.fit', 'cheng2020_M28.fit',
           'cheng2020_omg.fit', 'ngc6397_catalog.fits', 'ngc6752_catalog.fits',
           'NGC6266_p50_i5_src_1_2_4_8.fits']
label = ['.', '^', 'v', 'o', 'D', '*', 's']
color_list = ['grey', 'g', 'b', 'k', 'orange', 'purple', 'magenta']
result_all = pd.read_excel('/Users/baotong/Desktop/period_terzan5/candidate_allGC.xlsx', 'all')
ra = np.array(result_all['ra'])
dec = np.array(result_all['dec'])
seq = np.array(result_all['seq'])
period = np.array(result_all['period_all'])
type = np.array(result_all['GC'])
dist = np.array(result_all['proj_dist'])
counts = np.array(result_all['counts'])
exptime = np.array(result_all['expT'])
judge=np.array(result_all['judge'])
CR = counts / exptime

def f(x,a,b):
    logS=a*x+b  #S~V^a
    return logS
def spectrafit(x,y,error):
    popt, pcov = op.curve_fit(f, np.log(x), np.log(y),absolute_sigma=True,sigma=np.log(error))
    # popt, pcov = op.curve_fit(f, np.log(x), np.log(y))
    perr = np.sqrt(np.diag(pcov))
    logydata=f(np.log(x),popt[0],popt[1])
    ydata=np.exp(logydata)
    return (popt,perr)

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
    # plt.scatter(period/3600,CR)
    # plt.loglog()
    bins_rc = np.logspace(np.log10(1.0), np.log10(120), 11)
    fig = plt.figure(1, figsize=(9, 6))
    ax1 = fig.add_subplot(111)
    ax1.plot()
    label_gc = ['47 Tuc', 'Terzan 5', r'$\omega~cen$', 'M 28', 'NGC 6397', 'NGC 6752', 'NGC 6266', 'M 30']
    dist_rc = []
    for i in range(len(gcname) - 1):
        gc_srcid = np.where(type == gcname[i])[0]
        dist_rc = np.concatenate((dist_rc, dist[gc_srcid] / pos_all[gcname[i]][3]))
    hist_drc = ax1.hist(dist_rc, bins=bins_rc, histtype='step', lw=1.5, color='blue', linestyle='-')
    print(np.sum(hist_drc[0]))
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel(r'R/$r_{h}$', hawk.font1)
    ax1.set_ylabel('Number of sources', hawk.font1)
    ax1.tick_params(labelsize=16)
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


def plot_Phist(save=0, show=1):
    result_all = pd.read_excel('/Users/baotong/Desktop/period_terzan5/candidate_allGC.xlsx', 'all')
    result_allbutTuc = pd.read_excel('/Users/baotong/Desktop/period_terzan5/candidate_allGC.xlsx', 'allbutTuc')
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

    idex2 = np.where(result_allbutTuc['judge'] == 'CV')[0]
    period_allbutTuc = np.array(result_allbutTuc['period_all'])[idex2]
    period_allbutTuc /= 3600.
    period_GC = period / 3600.

    P_min = 7. / 6.
    P_gap = [7740.0 / 3600., 11448.0 / 3600.]

    bins = np.logspace(np.log10(0.5), np.log10(30), 41)
    bins_2 = np.logspace(np.log10(0.5), np.log10(100), 41)
    bins_spin = np.logspace(np.log10(3 / 36.), np.log10(2), 31)
    bins_p = np.logspace(np.log10(0.8), np.log10(16), 13)
    bins_p = np.concatenate((bins_p, [20., 25., 30.]))
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
    (ra, dec, seq, period_LW, L, Lmin, Lmax, type) = plot_pXS.load_LW('result_LW')
    period_LW = period_LW[np.where((period_LW > 3900) & (period_LW < 40000))]
    period_LW /= 3600.
    # fig, axes = plt.subplots(2, 1, figsize=(15, 10), sharex='all')
    fig, (ax1, ax2) = plt.subplots(ncols=1, nrows=2, sharex='all', gridspec_kw={'height_ratios': [2, 1]},
                                   figsize=(15, 10))
    # ax1=fig.add_subplot(211)
    ax1.plot()
    ax1.hist(orb_all, bins=bins, histtype='step', lw=2, color='red', linestyle='--')
    ax1.hist(period_allbutTuc, bins=bins_p, histtype='step', lw=2, linestyle='-', facecolor='grey',
             hatch='/', edgecolor='k', fill=True)
    ax1.hist(period_GC, bins=bins_p, histtype='step', lw=4, color='c', linestyle='-')
    ax1.hist(period_LW, bins=bins_p, histtype='step', lw=3, color='blue', linestyle='-')
    # ax1.hist(spin_IP, bins = bins_spin, histtype = 'step',lw=1.5, color = 'purple',linestyle='-')
    ax1.legend(['CVs in Solar Neighborhood', 'CVs in GCs (excluding 47 Tuc)', 'CVs in GCs', 'CVs in Galactic Bulge'])
    P_min = 1.373333333
    P_gap = [7740.0 / 3600., 11448.0 / 3600.]
    ax1.set_ylim(8e-1, 360)
    ax1.set_xlim(0.6, 35)
    ax1.plot([P_min, P_min], [0, 220], '--', color='grey')
    ax1.plot([P_gap[0], P_gap[0]], [0, 220], '-', lw=2., color='orange')
    ax1.plot([P_gap[1], P_gap[1]], [0, 220], '-', lw=2., color='orange')
    ax1.text(P_gap[0] + 0.3, 250, 'gap', fontsize=14)
    ax1.text(P_min - 0.2, 250, 'minimum', fontsize=14)
    # ax1.set_xlabel('Period (hours)',font1)
    ax1.set_ylabel('Number of sources', plot_pXS.font1)
    ax1.tick_params(labelsize=16)
    ax1.set_xscale('log')
    ax1.set_yscale('log')

    # ax2=fig.add_subplot(212)
    ax2.hist(orb_all, bins=bins, histtype='step', lw=2, color='red', cumulative=1, density=1, linestyle='--')
    ax2.hist(period_allbutTuc, bins=bins_p, histtype='step', lw=2, linestyle='-', cumulative=1, density=1,
             facecolor='grey',
             hatch='/', edgecolor='k', fill=True)
    ax2.hist(period_GC, bins=bins_p, histtype='step', lw=4, color='c', cumulative=1, density=1, linestyle='-')
    ax2.hist(period_LW, bins=bins_p, histtype='step', lw=3, color='blue', cumulative=1, density=1, linestyle='-')
    # ax2.hist(spin_IP, bins = bins_spin, histtype = 'step',lw=1.5, color = 'purple',cumulative=1,density=1,linestyle='-')

    ax2.plot([P_min, P_min], [0, 1], '--', color='grey')
    ax2.plot([P_gap[0], P_gap[0]], [0, 1], '-', lw=2., color='orange')
    ax2.plot([P_gap[1], P_gap[1]], [0, 1], '-', lw=2., color='orange')

    ax2.set_xscale('log')
    ax2.set_xlabel('Period (hours)', plot_pXS.font1)
    ax2.set_ylabel('CDF', plot_pXS.font1)
    ax2.tick_params(labelsize=16)
    plt.subplots_adjust(wspace=0, hspace=0.05)
    # ax2.set_yscale('log')
    path_out = '/Users/baotong/Desktop/aas/GCall/figure/'
    if save:
        plt.savefig(path_out + 'GC_NP.pdf', bbox_inches='tight', pad_inches=0.05)
    if show:
        plt.show()


def plot_P_L_profile(save=0, show=1):
    label = ['.', '^', 'v', 'o', 'D', '*', 's']
    color_list = ['grey', 'g', 'b', 'k', 'orange', 'purple', 'magenta']
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

    (ra_LW, dec_LW, seq_LW, period_LW, L_LW, Lmin, Lmax, type_LW) = plot_pXS.load_LW('result_LW')
    LW_index = np.where((period_LW > 3900) & (period_LW < 40000))[0]
    period_LW = period_LW[LW_index]
    period_LW /= 3600.
    L_LW = L_LW[LW_index]
    fig, f1_axes = plt.subplots(ncols=2, nrows=2, sharey='row', sharex='col',
                                gridspec_kw={'height_ratios': [1, 1], 'width_ratios': [3, 1]},
                                figsize=(12, 12))
    ax1 = f1_axes[0, 0];
    ax2 = f1_axes[1, 0];
    ax3 = f1_axes[0, 1];
    ax4 = f1_axes[1, 1]
    label_gc = ['47 Tuc', 'Terzan 5', r'$\omega~cen$', 'M 28', 'NGC 6397', 'NGC 6752', 'NGC 6266', 'M 30']
    dist_rc = []
    bins_lx = np.logspace(np.log10(1e31), np.log10(1e33), 11)
    bins_rc = np.logspace(np.log10(0.15), np.log10(130), 11)
    for i in range(len(gcname) - 1):
        gc_srcid = np.where(type == gcname[i])[0]
        dist_rc = np.concatenate((dist_rc, dist[gc_srcid] / pos_all[gcname[i]][3]))
        ax1.scatter(period[gc_srcid] / 3600, L[gc_srcid], marker=label[i], color=color_list[i], label=label_gc[i], s=70)
        ax2.scatter(period[gc_srcid] / 3600, dist[gc_srcid] / pos_all[gcname[i]][3], marker=label[i],
                    color=color_list[i], label=label_gc[i], s=70)
    ax1.scatter(period_LW, 1e31 * L_LW, marker='x', color='grey', label='LW', s=60)
    counts_lx, bins_lx = np.histogram(L, bins_lx)
    lx_drc = ax3.hist(bins_lx[:-1], bins=bins_lx, weights=counts_lx / np.sum(counts_lx), histtype='step',
                      orientation='horizontal', lw=1.5, color='k',
                      linestyle='--')
    counts_dist, bins_rc = np.histogram(dist_rc, bins_rc)
    dist_drc = ax4.hist(bins_rc[:-1], bins=bins_rc, weights=counts_dist / np.sum(counts_dist), histtype='step',
                        orientation='horizontal',
                        lw=1.5, color='k', linestyle='--')

    ax1.plot([P_min, P_min], [0, 1e33], '--', color='grey')
    ax1.plot([P_gap[0], P_gap[0]], [0, 1e33], '-', lw=2., color='orange')
    ax1.plot([P_gap[1], P_gap[1]], [0, 1e33], '-', lw=2., color='orange')
    ax1.set_ylabel(r'Luminosity ($\rm erg~s^{-1}$)', hawk.font1)
    (period_sim, L_sim, Mdot_sim) = CV_model.read_Knigge_model()
    # ax5 = ax1.twinx()
    # ax5.set_ylabel(r'$\rm Log\dot{M_2}$ ($\rm M_\odot~ \rm year^{-1}$)', hawk.font1)
    # ax5.tick_params(labelsize=18)
    ax1.scatter(period_sim, L_sim, s=10, marker='.', color='red', label='Model')

    # ax1.set_xlabel('Period (h)',hawk.font1)
    ax1.loglog()
    ax2.loglog()
    ax3.set_yscale('log')
    ax4.set_yscale('log')
    ax1.set_xlim(1, 28)
    ax1.tick_params(labelsize=18)
    ax2.plot([P_min, P_min], [0, 150], '--', color='grey')
    ax2.plot([P_gap[0], P_gap[0]], [0, 150], '-', lw=2., color='orange')
    ax2.plot([P_gap[1], P_gap[1]], [0, 150], '-', lw=2., color='orange')
    ax2.plot([1, 26], [4, 4], '-', lw=1., color='c')
    ax2.plot([1, 26], [10, 10], '-', lw=1., color='c')
    ax2.fill_between([1, 26], [4, 4], [10, 10], facecolor='yellow', alpha=0.2)
    ax2.set_xlabel('Period (h)', hawk.font1)
    ax2.set_ylabel(r'R/$r_{c}$', hawk.font1)
    ax2.set_xlim(1, 28)
    ax2.set_ylim(1e-1, 150)
    ax2.tick_params(labelsize=18)
    ax4.tick_params(labelsize=18)
    ax4.xaxis.set_major_formatter(ticker.PercentFormatter(xmax=1, decimals=0))
    ax1.legend(loc='upper right', ncol=1, handletextpad=0.001, columnspacing=0.001, facecolor=None, edgecolor='grey')
    # ax1.legend(loc='right',ncol=1,facecolor=None,edgecolor='grey')
    plt.subplots_adjust(wspace=0, hspace=0.0)

    path_out = '/Users/baotong/Desktop/aas/GCall/figure/'
    if save:
        plt.savefig(path_out + 'GC_CV_profile.pdf', bbox_inches='tight', pad_inches=0.05)
    if show:
        plt.show()
    plt.figure(2)
    return period


def plot_P_L(save=0, show=1):
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

    fig, ax1 = plt.subplots(ncols=1, nrows=1, figsize=(10, 10))
    # ax1 = fig.add_subplot(211)
    # ax2 = fig.add_subplot(212)
    for i in range(len(gcname) - 1):
        gc_srcid = np.where(type == gcname[i])[0]
        if len(gc_srcid) == 0:
            continue
        else:
            if i == 0:
                s = 80
            else:
                s = 120
            ax1.scatter(period[gc_srcid] / 3600, L[gc_srcid], marker=label[i], color=color_list[i], label=gcname[i],
                        s=s)

    mCVpars = np.array(list(localCVs.mCVs.values()))
    period_mCV = mCVpars[:, 0];
    L_mCV = mCVpars[:, 1]
    non_mCVpars = np.array(list(localCVs.non_mCVs.values()))
    period_non_mCV = non_mCVpars[:, 0];
    L_non_mCV = non_mCVpars[:, 1]
    ax1.scatter(period_mCV, 10 ** L_mCV, marker='<', color='c', label='local mCVs',
                s=60)
    ax1.scatter(period_non_mCV, 10 ** L_non_mCV, marker='>', color='c', label='local non-mCVs',
                s=60)
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
    ax1.legend(loc='lower center', ncol=len(gcname), handletextpad=0.1, columnspacing=0.1, handlelength=1.5,
               edgecolor='grey')
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


def plot_P_L_type(save=0, show=1):
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


def plot_P_L_threereg(save=0, show=1):
    result_GC = pd.read_excel('/Users/baotong/Desktop/period_terzan5/candidate_allGC.xlsx', 'all')
    idex_GC = np.where(result_GC['judge'] == 'CV')[0]
    period_GC = np.array(result_GC['period_all'])[idex_GC]
    type_GC = np.array(result_GC['GC'])[idex_GC]
    L_GC = np.array(result_GC['L'])[idex_GC]
    print(len(period_GC))

    (ra, dec, seq, period_LW, L_LW, Lmin, Lmax, type) = plot_pXS.load_LW('result_LW')
    # period_LW = period[np.where((period > 3900) & (period < 40000))]
    L_LW = L_LW * 1.11423 * 1e31
    path_table = '/Users/baotong/Desktop/period/table/'
    result_NSC = pd.read_excel(path_table + 'final_all_del.csv', 'result_NSC_IG')
    line = result_NSC['line']

    ID_NSC = result_NSC['seq']
    period_NSC = result_NSC['P']
    L_NSC = result_NSC['L']
    L_NSC = L_NSC * 1.3 * 1e31
    ID_NSC = ID_NSC[np.where(line == 3)[0]]
    period_NSC = period_NSC[np.where(line == 3)[0]]
    L_NSC = L_NSC[np.where(line == 3)[0]]
    ##-----P_L-----##
    # P_gap = [7740.0 / 3600., 11048.0 / 3600.]
    # P_min = [4902.0 / 3600., 4986.0/3600]
    # y1 = [0, 2e33]
    # # ax2.text(7900 / 3600., 5e36, 'period gap')
    # plt.text(P_gap[0] + 0.25, y1[1], r'$\rm P_{gap}$',fontsize=18)
    # plt.text(P_min[1] - 0.15, y1[1], r'$\rm P_{min}$',fontsize=18)
    # plt.ylim(ymin=8e30,ymax=1e33)
    # plt.xlim(xmin=0.95,xmax=20)
    #
    # plt.fill_between(P_gap, y1[1], facecolor='yellow', alpha=0.2)
    # plt.fill_between(P_min, y1[1], facecolor='grey', alpha=0.2)
    #
    # plt.scatter(period_GC/3600.,L_GC,marker='o',label='Globular cluster')
    # plt.scatter(period_LW/3600., L_LW, marker='x',label='LW')
    # plt.scatter(period_NSC/3600., L_NSC, marker='^',label='NSC')
    # plt.loglog()
    # plt.show()

    ##--------hist--------##
    path_fits = '/Users/baotong/Desktop/period_LW/'
    RK = fits.open(path_fits + 'RK14.fit')
    orb = RK[1].data['Orb_Per']
    type1 = RK[1].data['Type1'];
    type2 = RK[1].data['Type2'];
    type3 = RK[1].data['Type3']
    orb = orb * 24;
    spin = RK[1].data['_3___Per']
    orb_DN = orb[np.where(type1 == 'DN')]
    orb_Polar = orb[np.where((type2 == 'AM') | (type3 == 'AM'))]
    orb_IP = orb[np.union1d(np.where((type2 == 'IP') | (type2 == 'DQ')), np.where((type3 == 'IP') | (type3 == 'DQ')))]
    spin_IP = spin[np.union1d(np.where((type2 == 'IP') | (type2 == 'DQ')), np.where((type3 == 'IP') | (type3 == 'DQ')))]
    spin_IP /= 3600.
    orb_all = np.concatenate((orb_DN, orb_IP, orb_Polar))
    orb_all = orb_all[np.where(orb_all > 0)]

    bins = np.logspace(np.log10(0.5), np.log10(30), 41)
    bins_2 = np.logspace(np.log10(0.8), np.log10(20), 21)
    bins_spin = np.logspace(np.log10(3 / 36.), np.log10(2), 31)
    bins_p = np.logspace(np.log10(0.5), np.log10(12), 16)
    bins_p = np.concatenate((bins_p, [15, 20., 25., 30.]))
    bins_NSC = np.concatenate(([0.08, 0.2, 0.4], bins_p))

    fig, (ax1, ax2) = plt.subplots(ncols=1, nrows=2, sharex='all', gridspec_kw={'height_ratios': [2, 1]},
                                   figsize=(9, 6))
    ax1.hist(orb_all, bins=bins, histtype='step', lw=2, color='red', linestyle='--')
    ax1.hist(period_NSC / 3600, bins=bins_NSC, histtype='step', lw=4, color='blue', linestyle='-')
    ax1.hist(period_GC / 3600, bins=bins_p, histtype='step', lw=2, linestyle='-', facecolor='grey',
             hatch='/', edgecolor='k', fill=False)
    ax1.hist(period_LW / 3600, bins=bins_p, histtype='step', lw=3, color='c', linestyle='-')

    ax1.legend(['CVs in Solar Neighborhood', 'CVs in NSC', 'CVs in GCs', 'CVs in LW'])
    P_min = 1.373333333
    P_gap = [7740.0 / 3600., 11448.0 / 3600.]
    ax1.set_ylim(8e-1, 360)
    ax1.set_xlim(0.06, 60)
    ax1.plot([P_min, P_min], [0, 220], '--', color='grey')
    ax1.plot([P_gap[0], P_gap[0]], [0, 220], '-', lw=2., color='orange')
    ax1.plot([P_gap[1], P_gap[1]], [0, 220], '-', lw=2., color='orange')
    ax1.text(P_gap[0] + 0.3, 250, 'gap', fontsize=14)
    ax1.text(P_min - 0.2, 250, 'minimum', fontsize=14)
    # ax1.set_xlabel('Period (hours)',font1)
    ax1.set_ylabel('Number of sources', plot_pXS.font1)
    ax1.tick_params(labelsize=16)
    ax1.set_xscale('log')
    ax1.set_yscale('log')

    # ax2=fig.add_subplot(212)
    ax2.hist(orb_all, bins=bins, histtype='step', lw=2, color='red', cumulative=1, density=1, linestyle='--')
    ax2.hist(period_NSC / 3600, bins=bins_NSC, histtype='step', lw=4, color='blue', cumulative=1, density=1,
             linestyle='-')
    ax2.hist(period_GC / 3600, bins=bins_p, histtype='step', lw=2, linestyle='-', cumulative=1, density=1,
             facecolor='grey',
             hatch='/', edgecolor='k', fill=False)
    ax2.hist(period_LW / 3600, bins=bins_p, histtype='step', lw=3, color='c', cumulative=1, density=1, linestyle='-')
    ax2.tick_params(labelsize=16)
    ax2.set_ylabel('CDF', plot_pXS.font1)
    ax2.set_xlabel('Period (hour)', plot_pXS.font1)
    plt.show()


def plot_spec():
    label = ['x', '^', 'v', 'o', 'D', '*', 's']
    color_list = ['r', 'g', 'b', 'k', 'orange', 'purple', 'magenta']
    result_all = pd.read_excel('/Users/baotong/Desktop/period_terzan5/candidate_allGC.xlsx', 'allbutTuc')
    ra = np.array(result_all['ra'])
    dec = np.array(result_all['dec'])
    seq = np.array(result_all['seq'])
    period = np.array(result_all['period_all'])
    type = np.array(result_all['GC'])
    dist = np.array(result_all['proj_dist'])
    counts = np.array(result_all['counts'])
    exptime = np.array(result_all['expT'])
    L = np.array(result_all['L'])
    kT = np.array(result_all['kT'])
    plt.figure(1)
    plt.scatter(kT, L)
    plt.loglog()
    plt.show()


def plot_src_info(save=0, show=1):
    label = ['x', '^', 'v', 'o', 'D', '*', 's', '.']
    color_list = ['r', 'g', 'b', 'k', 'orange', 'purple', 'magenta', 'cyan']
    for i in range(len(gcname)):
        path = '/Users/baotong/Desktop/period_' + gcname[i] + '/'
        src_info = np.loadtxt(path + 'src_info.txt')
        counts_all = src_info[:, 3];
        exptime_all = src_info[:, 4];
        VI_all = src_info[:, 5]
        print('bright_src=', len(np.where(counts_all > 100)[0]))
        bins_counts = np.logspace(0.5, 4, 20)
        plt.hist(counts_all, bins=bins_counts, histtype='step', lw=2, linestyle='-', color=color_list[i],
                 label=gcname[i])
        plt.semilogx()
        # plt.show()
        # plt.scatter(exptime_all,VI_all)
        # plt.semilogy()
        # plt.show()
    plt.legend()
    plt.show()


def plot_P_L_4reg_scatter(save=0, show=1):
    result_GC = pd.read_excel('/Users/baotong/Desktop/period_terzan5/candidate_allGC.xlsx', 'all')
    idex_GC = np.where(result_GC['judge'] == 'CV')[0]
    period_GC = np.array(result_GC['period_all'])[idex_GC]
    type_GC = np.array(result_GC['GC'])[idex_GC]
    L_GC = np.array(result_GC['L'])[idex_GC]
    print(len(period_GC))

    (ra, dec, seq, period_LW, L_LW, Lmin, Lmax, type) = plot_pXS.load_LW('result_LW')
    # period_LW = period[np.where((period > 3900) & (period < 40000))]
    L_LW = L_LW * 1.11423 * 1e31
    idex = np.where(type == 'CV')[0]
    seq = seq[idex]
    period_LW = period_LW[idex]
    L_LW = L_LW[idex]

    path_table = '/Users/baotong/Desktop/period/table/'
    result_NSC = pd.read_excel(path_table + 'final_all_del.csv', 'result_NSC_IG')
    line = result_NSC['line']
    ID_NSC = result_NSC['seq']
    period_NSC = result_NSC['P']
    L_NSC = result_NSC['L']
    L_NSC = L_NSC * 1.3 * 1e31
    ID_NSC = ID_NSC[np.where(line == 3)[0]]
    period_NSC = period_NSC[np.where(line == 3)[0]]
    L_NSC = L_NSC[np.where(line == 3)[0]]

    mCVpars = np.array(list(localCVs.mCVs.values()))
    period_mCV = mCVpars[:, 0];
    L_mCV = mCVpars[:, 1]
    non_mCVpars = np.array(list(localCVs.non_mCVs.values()))
    period_non_mCV = non_mCVpars[:, 0];
    L_non_mCV = non_mCVpars[:, 1]
    # -----P_L-----##
    P_gap = [7740.0 / 3600., 11048.0 / 3600.]
    P_min = [4902.0 / 3600., 4986.0 / 3600]
    y1 = [0, 2e33]
    # ax2.text(7900 / 3600., 5e36, 'period gap')
    fig, ax1 = plt.subplots(ncols=1, nrows=1, figsize=(8, 10))

    ax1.text(P_gap[0] + 0.25, y1[1], r'$\rm P_{gap}$', fontsize=18)
    ax1.text(P_min[1] - 0.15, y1[1], r'$\rm P_{min}$', fontsize=18)
    # plt.ylim(ymin=8e30,ymax=1e34)
    ax1.set_xlim(xmin=0.95, xmax=20)

    ax1.fill_between(P_gap, y1[1], facecolor='yellow', alpha=0.2)
    ax1.fill_between(P_min, y1[1], facecolor='grey', alpha=0.2)

    ax1.scatter(period_GC / 3600., L_GC, s=60, marker='o', facecolor='none', edgecolors='blue',
                label='Globular cluster')
    ax1.scatter(period_LW / 3600., L_LW, s=60, marker='x', color='k', edgecolors='k', label='LW')
    ax1.scatter(period_NSC / 3600., L_NSC, s=60, marker='^', color='c', edgecolors='c', label='NSC')
    ax1.scatter(period_mCV, 2.25 * 10 ** L_mCV, s=60, marker='s', color='magenta', edgecolors='magenta',
                label='lcoal mCVs')
    ax1.scatter(period_non_mCV, 1.76 * 10 ** L_non_mCV, s=60, marker='v', color='grey', edgecolors='grey',
                label='local non mCVs')
    (period_sim, L_sim, Mdot_sim) = CV_model.read_Knigge_model()
    ax2 = ax1.twinx()
    ax2.set_ylabel(r'$\rm Log\dot{M_2}$ ($\rm M_\odot~ \rm year^{-1}$)', hawk.font1)
    ax1.scatter(period_sim, 5 * L_sim, s=30, marker='.', color='red')
    ax2.tick_params(labelsize=18)

    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.legend()
    path_out = '/Users/baotong/Desktop/aas/GCall/figure/'
    if save:
        plt.savefig(path_out + 'reg4_CV_P_L.pdf', bbox_inches='tight', pad_inches=0.05)
    if show:
        plt.show()
    return None

def plot_surface_density(save=0,show=1):

    dist_rc = [];dist_all=[];dist_rc_core=[];dist_all_core=[]

    for i in range(len(gcname)):
        gc_srcid = np.where((judge=='CV')&(type == gcname[i]))[0]
        if i<3:  ##去掉omg cen
            dist_rc = np.concatenate((dist_rc, dist[gc_srcid] / pos_all[gcname[i]][2]))
            srcinfo = np.loadtxt(f'/Users/baotong/Desktop/period_{gcname[i]}/src_info.txt')
            srcinfo=srcinfo[np.where(srcinfo[:,3]>100)]
            if gcname[i]=='Tuc':dist_single=srcinfo[:,-1]/pos_all[gcname[i]][2]
            else:dist_single=srcinfo[:,6]/pos_all[gcname[i]][2]
            dist_all=np.concatenate((dist_all,dist_single))
        if i>=4:
            dist_rc_core = np.concatenate((dist_rc_core, dist[gc_srcid] / pos_all[gcname[i]][2]))
            srcinfo = np.loadtxt(f'/Users/baotong/Desktop/period_{gcname[i]}/src_info.txt')
            srcinfo=srcinfo[np.where(srcinfo[:,3]>100)]
            if gcname[i]=='Tuc':dist_single=srcinfo[:,-1]/pos_all[gcname[i]][2]
            else:dist_single=srcinfo[:,6]/pos_all[gcname[i]][2]
            dist_all_core=np.concatenate((dist_all_core,dist_single))

    bins_rc=np.logspace(np.log10(0.01),np.log10(7),6)
    bins2=np.logspace(np.log10(0.01),np.log10(10),11)
    x1=[(bins_rc[i]+bins_rc[i+1])/2 for i in range(len(bins_rc)-1)]
    x1err = [(bins_rc[i + 1] - bins_rc[i]) / 2 for i in range(len(bins_rc) - 1)]
    area1 = [(bins_rc[i + 1] ** 2 - bins_rc[i] ** 2) * np.pi for i in range(len(bins_rc) - 1)]
    x2 = [(bins2[i] + bins2[i + 1]) / 2 for i in range(len(bins2) - 1)]
    area2 = [(bins2[i + 1] ** 2 - bins2[i] ** 2) * np.pi for i in range(len(bins2) - 1)]
    x2err = [(bins2[i + 1] - bins2[i]) / 2 for i in range(len(bins2) - 1)]

    hist1=plt.hist(dist_rc,bins_rc)
    hist2=plt.hist(dist_all,bins2)
    hist3=plt.hist(dist_rc_core,bins_rc)
    hist4=plt.hist(dist_all_core,bins2)
    plt.close()
    y1=hist1[0];y2=hist2[0];y3=hist3[0];y4=hist4[0]
    print(dist_rc,dist_rc_core)
    print(np.sum(y1),np.sum(y3))
    y1_err = np.array(poisson_conf_interval(y1, interval='frequentist-confidence'))
    y1_err[0] = y1 - y1_err[0];y1_err[1] = y1_err[1] - y1
    y2_err = np.array(poisson_conf_interval(y2, interval='frequentist-confidence'))
    y2_err[0] = y2 - y2_err[0];
    y2_err[1] = y2_err[1] - y2
    y3_err = np.array(poisson_conf_interval(y3, interval='frequentist-confidence'))
    y3_err[0] = y3 - y3_err[0];
    y3_err[1] = y3_err[1] - y3
    y4_err = np.array(poisson_conf_interval(y4, interval='frequentist-confidence'))
    y4_err[0] = y4 - y4_err[0];
    y4_err[1] = y4_err[1] - y4

    y1=y1/area1;y1_err=y1_err/area1;y3=y3/area1;y3_err=y3_err/area1
    y2 = y2 / area2;y2_err = y2_err / area2;    y4= y4/ area2;y4_err = y4_err / area2
    # plt.hist(dist_rc,bins=bins,histtype='step')
    plt.figure(1,(10,8))
    plt.errorbar(x1,y1*20,xerr=x1err,yerr=y1_err*20,fmt='ro', capsize=3, elinewidth=3, ecolor='r', color='r',
                 markersize=4,label='Periodic CVs (normal)')
    # plt.errorbar(x2, y2, xerr=x2err, yerr=y2_err, fmt='ko', capsize=1, elinewidth=1, ecolor='k', color='k',
                 # markersize=4, label='Total (normal)')
    x1=np.array(x1);x1err=np.array(x1err)
    mask=np.where(y3>0)[0]
    print('mask=',mask)
    y3_err=np.array([y3_err[0][mask],y3_err[1][mask]])
    plt.errorbar(x1[mask],y3[mask]*20,xerr=x1err[mask],yerr=y3_err*20,fmt='ro', capsize=3, elinewidth=3, ecolor='g', color='g',
                 markersize=4,label='Periodic CVs (core-collapse)')
    # plt.errorbar(x2, y4, xerr=x2err, yerr=y4_err, fmt='ko', capsize=1, elinewidth=1, ecolor='grey', color='grey',
    #              markersize=4, label='Total (core-collapse)')

    (popt1, perr1) = spectrafit(x1, y1, y1_err[0])
    (popt2, perr2) = spectrafit(x1[mask], y3[mask], y3_err[0])

    plt.plot(x1, 10 * np.exp(f(np.log(x1), popt1[0], popt1[1])), '-.', color='k')
    plt.plot(x1, 20 *np.exp(f(np.log(x1), popt2[0], popt2[1])), '-.', color='r')

    # 定义插值函数
    f2 = interp1d(x2, y2-2, kind='linear')  # cubic 表示使用三次样条插值
    f4 = interp1d(x2, y4-2, kind='linear')  # cubic 表示使用三次样条插值
    # 生成平滑的曲线
    bins2_smooth = np.logspace(np.log10(0.01),np.log10(10),20)
    x2_smooth =np.array([(bins2_smooth[i] + bins2_smooth[i + 1]) / 2 for i in range(len(bins2_smooth) - 1)])
    x2_smooth=x2_smooth[np.where((x2_smooth<np.max(x2))&(x2_smooth>np.min(x2)))[0]]
    y2_smooth = f2(x2_smooth)
    y4_smooth = f4(x2_smooth)
    plt.plot(x2_smooth, y2_smooth, 'b.-', label='Smooth Curve')
    plt.plot(x2_smooth, y4_smooth, 'g.-',label='Smooth Curve')
    print('popt1=', popt1)
    print('perr1=', perr1)
    print('popt2=', popt2)
    print('perr2=', perr2)
    print('y1=',y1)
    print('y3=',y3)
    plt.legend()
    plt.semilogx()
    plt.semilogy()
    plt.show()

def read_GLres_src(save=0,show=1):
    fig, f1_axes = plt.subplots(ncols=2, nrows=3, sharey='row', sharex='col',
                                gridspec_kw={'height_ratios': [1, 1,1], 'width_ratios': [1, 1]},
                                figsize=(12, 15))
    f1_axes=f1_axes.flatten()
    for k in range(1,len(gcname)):
        ax = f1_axes[k-1]
        label_gc = ['47 Tuc', 'Terzan 5', 'M 28',r'$\omega~cen$', 'NGC 6397', 'NGC 6752', 'NGC 6266']
        srcid = seq[np.where(type == gcname[k])[0]]
        path = '/Users/baotong/Desktop/period_' + gcname[k] + '/'
        for j in range(len(srcid)):
            prob_GL = []
            for i in range(100):
                if os.path.isfile(path+'sim_GL/'+'result_3h_{0}_f{1}.txt'.format(str(srcid[j]),i)):
                    res=np.loadtxt(path+'sim_GL/'+'result_3h_{0}_f{1}.txt'.format(str(srcid[j]),i))
                elif os.path.isfile(path+'sim_GL/'+'result_10h_{0}_f{1}.txt'.format(str(srcid[j]),i)):
                    res=np.loadtxt(path+'sim_GL/'+'result_10h_{0}_f{1}.txt'.format(str(srcid[j]),i))
                else:
                    print('None')
                    continue
                prob_GL.append(res[2])
            prob_GL=np.array(prob_GL)
            bins = np.linspace(0, 1, 21)
            ax.hist(prob_GL,bins=20,histtype='step',label=f'#{j+1}',lw=2)
            fD=len(prob_GL[prob_GL>0.95])
            # print(fD)
            # print(np.where(prob_GL>0.9)[0])
        # plt.semilogx()
        ax.plot([0.95, 0.95], [0, 100], '--',color='grey')
        ax.text(0.5,60,f'{label_gc[k]}',hawk.font1)
        # ax.legend()
        if k==1:ax.legend(loc='center', ncol=3, handletextpad=0.2, columnspacing=0.1, facecolor=None,
                   edgecolor='grey',handlelength=0.5)
        elif k==2:ax.legend(loc='upper right', ncol=2, handletextpad=0.2, columnspacing=0.1, facecolor=None,
                   edgecolor='grey',handlelength=0.5)
        else:ax.legend(loc='upper right', ncol=1, handletextpad=0.2, columnspacing=0.1, facecolor=None,
                   edgecolor='grey',handlelength=0.5)
        # plt.xlabel(r'$P_{\rm GL}$',hawk.font1)
        # plt.ylabel('Number of trials',hawk.font1)
        # plt.tick_params(labelsize=16)
        ax.set_yscale('log')
        if k==5 or k==6:ax.set_xlabel(r'$P_{\rm GL}$',hawk.font1)
        if k==1 or k==3 or k==5:ax.set_ylabel('Number of trials',hawk.font1)
        ax.tick_params(labelsize=16)
    plt.subplots_adjust(wspace=0, hspace=0.0)
    if save:plt.savefig('/Users/baotong/Desktop/aas/GCall/figure/'+'hist_sim.pdf',bbox_inches='tight', pad_inches=0.05)
    if show:plt.show()
    else:plt.close()
if __name__ == '__main__':
    # plot_profile()
    # plot_P_L_profile(save=1, show=1)
    # plot_P_L(save=0,show=1)
    # plot_Phist(save=0,show=1)
    # plot_P_L_4reg_scatter(save=0,show=1)
    # plot_profile(save=0,show=1)
    # plot_Phist(save=1,show=1)
    # plot_src_info(show=1)
    # plot_P_L_threereg()
    plot_surface_density(save=0,show=1)
    # read_GLres_src(save=1,show=1)