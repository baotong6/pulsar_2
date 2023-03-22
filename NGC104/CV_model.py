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
from astropy import constants as c
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
color_list = ['c', 'g', 'b', 'k', 'orange', 'purple', 'magenta']


def WD_M_R(M):
    # 输入白矮星质量，输出白矮星半径
    R = 7.8 * 1e8 * ((1.44 / M) ** (2.0 / 3.0) - (M / 1.44) ** (2.0 / 3.0)) ** 0.5
    return R


def read_Knigge_model():
    path = '/Users/baotong/Desktop/period_Tuc/'
    file = fits.open(path + 'Knigge2011_revised.fit')
    Mdot = file[1].data['logMLR2']
    index = np.where(Mdot < 0)[0]
    period = file[1].data['Per'][index]
    Mdot = file[1].data['logMLR2'][index]
    print(Mdot)
    M1 = file[1].data['M1'][index]
    R1 = WD_M_R(M1)
    epsilon = 0.5
    eta = 10
    L = 3.8e30 * (epsilon / 0.5) * (M1 / 1) * (10 ** Mdot / (eta * 1e-12))

    return (period, L, Mdot)

def load_LW(label):
    res = pd.read_excel('/Users/baotong/Desktop/period_LW/' + 'final_all_del_add.xlsx', label)
    ra = np.array(res['ra'])
    dec = np.array(res['dec'])
    seq = np.array(res['seq'])
    period = np.array(res['P'])
    L = np.array(res['L_ast'])
    Lmin =L+ np.array(res['L_low'])
    Lmax =L+ np.array(res['L_high'])
    type = np.array(res['type'])
    return (ra,dec,seq,period,L,Lmin,Lmax,type)

def plot_P_L(save=0,show=1):
    result_all = pd.read_excel('/Users/baotong/Desktop/period_terzan5/candidate_allGC.xlsx', 'all')
    idex = np.where(result_all['judge'] == 'CV')[0]
    period = np.array(result_all['period_all'])[idex]
    type = np.array(result_all['GC'])[idex]
    # ra = np.array(result_all['ra'])[idex]
    # dec = np.array(result_all['dec'])[idex]
    # seq = np.array(result_all['seq'])[idex]
    # dist = np.array(result_all['proj_dist'])[idex]
    # counts = np.array(result_all['counts'])[idex]
    # exptime = np.array(result_all['expT'])[idex]
    L = np.array(result_all['L'])[idex]
    P_min = 7. / 6.
    P_gap = [7740.0 / 3600., 11448.0 / 3600.]
    print(len(L))

    (ra, dec, seq, period1, L1, Lmin1, Lmax1, type1) = load_LW('result_LW')
    period_LW_CV=period1[np.where(type1=='CV')]
    L_LW_CV=L1[np.where(type1=='CV')]
    L_LW_CV=L_LW_CV*1.11423

    fig, ax1 = plt.subplots(ncols=1, nrows=1,figsize=(8, 10))
    for i in range(len(gcname) - 1):
        gc_srcid = np.where(type == gcname[i])[0]
        if len(gc_srcid)==0:continue
        else:
            # if i == 0:
            #     s = 50
            # else:
            #     s = 80
            ax1.scatter(period[gc_srcid] / 3600, L[gc_srcid], marker=label[i], color=color_list[i], label=gcname[i],
                    s=80)
    ax1.scatter(period_LW_CV / 3600., L_LW_CV * 1e31, marker='x', s=80, color='grey', linewidths=2,label='LW')
    ax1.plot([P_min, P_min], [0, 1e33], '--', color='grey')
    ax1.plot([P_gap[0], P_gap[0]], [0, 1e33], '-', lw=2., color='orange')
    ax1.plot([P_gap[1], P_gap[1]], [0, 1e33], '-', lw=2., color='orange')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_ylabel(r'Luminosity ($\rm erg~s^{-1}$)', hawk.font1)
    # ax1.set_xlabel('Period (h)',hawk.font1)
    ax1.set_xlim(0.9, 14)
    ax1.set_ylim(1e31, 5e34)
    ax1.tick_params(labelsize=18)
    ax1.set_xlabel('Period (h)', hawk.font1)
    ax1.tick_params(labelsize=18)
    ax1.legend(loc='upper right',bbox_to_anchor=(1,0.7), ncol=4, handletextpad=0.1, columnspacing=0.1, handlelength=1.5, edgecolor='grey')
    plt.subplots_adjust(wspace=0, hspace=0.05)
    (period_sim,L_sim,Mdot_sim)=read_Knigge_model()
    # ax2 = ax1.twinx()
    # ax2.set_ylabel(r'$\rm Log\dot{M_2}$ ($\rm M_\odot~ \rm year^{-1}$)', hawk.font1)
    # ax2.scatter(period_sim,Mdot_sim,s=30,marker='.',color='red')
    # ax2.tick_params(labelsize=18)

    path_out = '/Users/baotong/Desktop/aas/GCall/figure/'
    if save:
        plt.savefig(path_out + 'GC_CV_PL_sim2.pdf', bbox_inches='tight', pad_inches=0.05)
    if show:
            plt.show()
    return period
plot_P_L(save=1,show=1)