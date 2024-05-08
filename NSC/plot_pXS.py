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
from matplotlib import cm
from astropy import units as u
from scipy import optimize as op
import os
import NGC104.plot_pXS as plot_pXS
import GC.localCVs as localCVs
import NGC104.CV_model as CV_model
import GC.SDSSCV as SDSS
import GC.MESA_sim as MESA
import hawkeye as hawk
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
    label = result_NSC['label']
    line=result_NSC['line']
    ID_NSC = result_NSC['seq']
    period_NSC = result_NSC['P']
    L_NSC = result_NSC['L']
    L_NSC = L_NSC * 1.3 * 1e31
    ID_NSC = ID_NSC[np.where(label == 1)[0]]
    period_NSC = period_NSC[np.where(label == 1)[0]]
    L_NSC = L_NSC[np.where(label == 1)[0]]
    period_NSC_CV = period_NSC[np.where(line == 3)[0]]
    L_NSC_CV = L_NSC[np.where(line == 3)[0]]
    print(len(period_NSC),len(period_NSC_CV))
    ##-----P_L-----##
    P_gap = [7740.0 / 3600., 11048.0 / 3600.]
    P_min = [4902.0 / 3600., 4986.0/3600]
    y1 = [0, 2e33]
    # ax2.text(7900 / 3600., 5e36, 'period gap')
    plt.text(P_gap[0] + 0.25, y1[1], r'$\rm P_{gap}$',fontsize=18)
    plt.text(P_min[1] - 0.15, y1[1], r'$\rm P_{min}$',fontsize=18)
    plt.ylim(ymin=8e30,ymax=2e33)
    plt.xlim(xmin=0.1,xmax=20)
    
    plt.fill_between(P_gap, y1[1], facecolor='yellow', alpha=0.2)
    plt.fill_between(P_min, y1[1], facecolor='grey', alpha=0.2)
    
    plt.scatter(period_GC/3600.,L_GC,marker='s',color='blue',s=80,label='Globular cluster')
    plt.scatter(period_LW/3600., L_LW, marker='^',color='orange',s=80,label='LW')
    plt.scatter(period_NSC/3600., L_NSC, marker='x',color='red',s=80,label='NSC')
    plt.scatter(period_NSC_CV/3600., L_NSC_CV, marker='o',facecolor='none',edgecolors='black',s=110,label='NSC_CV')
    plt.loglog()
    plt.legend()
    plt.show()

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
    bins_p = np.logspace(np.log10(0.5), np.log10(12), 18)
    bins_GC=np.logspace(np.log10(0.8), np.log10(15), 13)
    bins_p = np.concatenate((bins_p, [15, 20., 25., 30.]))
    bins_NSC = np.concatenate(([0.08, 0.2, 0.4], bins_p))

    fig, (ax1, ax2) = plt.subplots(ncols=1, nrows=2, sharex='all', gridspec_kw={'height_ratios': [2, 1]},
                                   figsize=(9, 6))
    ax1.hist(orb_all, bins=bins, histtype='step', lw=2, color='red', linestyle='--',label='CVs in Solar Neighborhood')
    ax1.hist(period_GC / 3600, bins=bins_GC, histtype='step', lw=4, color='green',label='CVs in GCs')
    ax1.hist(period_LW / 3600, bins=bins_p, histtype='step', lw=3, color='c', linestyle='-',label='CVs in Galactic Bulge')
    ax1.hist(period_NSC / 3600, bins=bins_NSC, histtype='step', lw=4, color='blue', linestyle='-',label='Periodic X-ray sources in NSC')
    ax1.hist(period_NSC_CV / 3600, bins=bins_NSC, histtype='step', lw=2, linestyle='-', facecolor='grey',
             hatch='/', edgecolor='k', fill=True,label='CVs in NSC')
    ax1.legend()
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
    ax2.hist(period_GC / 3600, bins=bins_GC, histtype='step', lw=4, color='green', linestyle='-', cumulative=1, density=1)
    ax2.hist(period_NSC_CV / 3600, bins=bins_NSC, histtype='step', lw=2, linestyle='-', facecolor='grey',
             hatch='/', edgecolor='k', fill=True,cumulative=1, density=1)
    ax2.hist(period_LW / 3600, bins=bins_p, histtype='step', lw=3, color='c', cumulative=1, density=1, linestyle='-')
    ax2.tick_params(labelsize=16)
    ax2.set_ylabel('CDF', plot_pXS.font1)
    ax2.set_xlabel('Period (hour)', plot_pXS.font1)
    plt.subplots_adjust(wspace=0, hspace=0.0)
    plt.show()

def ratio_error_with_limits(mean1, low_limit1, up_limit1, mean2, low_limit2, up_limit2):
    # 计算比值
    ratio = mean1 / mean2

    rel_error1 = max((mean1 - low_limit1) / mean1, (up_limit1 - mean1) / mean1)
    rel_error2 = max((mean2 - low_limit2) / mean2, (up_limit2 - mean2) / mean2)
    ratio_error = ratio * ((rel_error1)**2 + (rel_error2)**2)**0.5
    
    return ratio_error
def plot_eqw():
    line_ew64=[0.125229,0.152686]
    line_ew67=[0.189578,0.222416]
    line_ew70=[0.124765,0.161709]
    noline_ew64=[0, 0.114829]
    noline_ew67=[0.0921115,0.266668]
    noline_ew70=[0.0886128, 0.26263]
    line_int64=[1.45922e-07,1.69832E-07,1.93737e-07]
    line_int67=[2.21071e-07,2.49046E-07, 2.77412e-07]
    line_int70=[1.14169e-07,1.40764E-07,1.67041e-07]
    noline_int64=[0,6.89436E-09,1.64936e-08]
    noline_int67=[5.08082e-09,1.67537E-08,2.83066e-08]
    noline_int70=[2.77702e-09,1.45960E-08,2.63095e-08]

    SS_ew67=[241,78.3];IP_ew67=[107,16];Polar_ew67=[221,135]
    DNe_ew67=[438,84.6];ABs_ew67=[286,58.5];GRXE_ew67=[490,15]

    SS_I70on67=[0.45,0.11]
    IP_I70on67=[0.71,0.04]
    Polar_I70on67=[0.44,0.14]
    DNe_I70on67=[0.27,0.06]
    ABs_I70on67=[0.08,0.04]
    GRXE_I70on67=[0.2,0.08]
    errline_70on67=ratio_error_with_limits(line_int70[1],line_int70[0],line_int70[2],
                                           line_int67[1],line_int67[0],line_int67[2])
    errnoline_70on67=ratio_error_with_limits(noline_int70[1],noline_int70[0],noline_int70[2],
                                             noline_int67[1],noline_int67[0],noline_int67[2])
    line_I70on67=[line_int70[1]/line_int67[1],errline_70on67]
    noline_I70on67=[noline_int70[1]/noline_int67[1],errnoline_70on67]

    plt.errorbar(x=0.5*(line_ew67[0]+line_ew67[1])*1000,y=line_I70on67[0],
                 xerr=1000*(line_ew67[1]-line_ew67[0])*0.5,yerr=line_I70on67[1],fmt='o', capsize=5)
    plt.text(0.5*(line_ew67[0]+line_ew67[1])*1000,line_I70on67[0]*1.2,'CVs in NSC',fontsize=18)
    plt.errorbar(x=0.5*(noline_ew67[0]+noline_ew67[1])*1000,y=noline_I70on67[0],
                 xerr=1000*(noline_ew67[1]-noline_ew67[0])*0.5,yerr=noline_I70on67[1],fmt='o', capsize=5)
    plt.text(0.5*(noline_ew67[0]+noline_ew67[1])*1000,
             noline_I70on67[0]*1.2,'Unknown sources in NSC',fontsize=18)

    plt.errorbar(x=SS_ew67[0],y=SS_I70on67[0],xerr=SS_ew67[1],yerr=SS_I70on67[1],fmt='o', capsize=5)
    plt.text(x=SS_ew67[0]*1.2,y=SS_I70on67[0]*1.2,s='SSs',fontsize=18)
    plt.errorbar(x=IP_ew67[0],y=IP_I70on67[0],xerr=IP_ew67[1],yerr=IP_I70on67[1],fmt='o', capsize=5)
    plt.text(x=IP_ew67[0],y=IP_I70on67[0]*1.2,s='IPs',fontsize=18)
    plt.errorbar(x=Polar_ew67[0],y=Polar_I70on67[0],xerr=Polar_ew67[1],yerr=Polar_I70on67[1],fmt='o', capsize=5)
    plt.text(x=Polar_ew67[0],y=Polar_I70on67[0]*1.2,s='Polars',fontsize=18)
    plt.errorbar(x=DNe_ew67[0],y=DNe_I70on67[0],xerr=DNe_ew67[1],yerr=DNe_I70on67[1],fmt='o', capsize=5)
    plt.text(x=DNe_ew67[0],y=DNe_I70on67[0]*1.2,s='DNes',fontsize=18)
    plt.errorbar(x=ABs_ew67[0],y=ABs_I70on67[0],xerr=ABs_ew67[1],yerr=ABs_I70on67[1],fmt='o', capsize=5)
    plt.text(x=ABs_ew67[0],y=ABs_I70on67[0]*1.2,s='ABs',fontsize=18)
    plt.errorbar(x=GRXE_ew67[0],y=GRXE_I70on67[0],xerr=GRXE_ew67[1],yerr=GRXE_I70on67[1],fmt='o', capsize=5)
    plt.text(x=GRXE_ew67[0],y=GRXE_I70on67[0]*1.2,s='GRXE',fontsize=18)
    plt.xlabel('EW 6.7 (keV)',hawk.font1)
    plt.ylabel('I7.0/I6.7',hawk.font1)
    plt.tick_params(labelsize=18)
    plt.show()


if __name__=='__main__':
    # plot_P_L_threereg()
    plot_eqw()
