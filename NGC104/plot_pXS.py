#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import pandas as pd
from astropy.stats import poisson_conf_interval
from astropy.coordinates import SkyCoord
from astropy import units as u
import scipy
import hawkeye as hawk
from timing_comb import load_data
font1 = {'family': 'Normal',
         'weight': 'normal',
         'size': 18, }
plt.rc('legend',fontsize=14 )

label_all = ['47Tuc', 'terzan5', 'M28', 'omg_cen', 'NGC6397', 'NGC6752']
# pos_all = [[6.0236250, -72.0812833, 3.17 * 60, 3.17 / 8.8 * 60],  # 47Tuc
#            [267.0202083, -24.7790556, 0.72 * 60, 0.72 / 3.4 * 60],  # terzan5
#            [276.1363750, -24.8702972, 1.97 * 60, 1.97 / 8.2 * 60],  # M28
#            [201.69700, -47.47947, 5 * 60, 5 / 2.1 * 60],  # omega_cen
#            [265.17539, -53.67433, 2.9 * 60, 2.9 / 58 * 60],  # NGC 6397
#            [287.71713, -59.98455, 1.91, 1.91 / 11.24 * 60]]  # NGC 6752
pos_all={'47Tuc':[6.0236250, -72.0812833, 3.17 * 60, 3.17 / 8.8 * 60],
         'terzan5':[267.0202083, -24.7790556, 0.72 * 60, 0.72 / 3.4 * 60],
         'M28':[276.1363750, -24.8702972, 1.97 * 60, 1.97 / 8.2 * 60],
         'omega_cen':[201.69700, -47.47947, 5 * 60, 5 / 2.1 * 60],
         'NGC6397':[265.17539, -53.67433, 2.9 * 60, 2.9 / 58 * 60],
         'NGC6752':[287.71713, -59.98455, 1.91, 1.91 / 11.24 * 60]}

path='/Users/baotong/Desktop/period_Tuc/'
path_out = '/Users/baotong/Desktop/aas/pXS_Tuc/figure/'

def read_excel(label):
    res = pd.read_excel(path + 'result_0.5_8_all.xlsx', label)
    ra = np.array(res['RA'])
    dec = np.array(res['DEC'])
    seq = np.array(res['Seq'])
    period = np.array(res['P_out'])
    L = np.array(res['L'])
    Lmin = np.array(res['Lmin'])
    Lmax = np.array(res['Lmax'])
    type = np.array(res['type'])

    return (ra,dec,seq,period,L,Lmin,Lmax,type)

def load_LW(label):
    res = pd.read_excel('/Users/baotong/Desktop/period_LW/' + 'final_all_del_add.xlsx', label)
    ra = np.array(res['ra'])
    dec = np.array(res['dec'])
    seq = np.array(res['seq'])
    period = np.array(res['P'])
    L = np.array(res['L_ast'])
    Lmin =L+ np.array(res['L_low'])
    Lmax =L+ np.array(res['L_high'])
    # type = np.array(res['type'])
    return (ra,dec,seq,period,L,Lmin,Lmax)
def plot_RK_CV():
    bins=np.logspace(np.log10(0.5), np.log10(50), 71)
    bins_2=np.logspace(np.log10(0.1), np.log10(20), 31)
    bins_spin=np.logspace(np.log10(3/36.), np.log10(2), 21)

    path_fits = '/Users/baotong/Desktop/period_LW/'
    RK = fits.open(path_fits + 'RK14.fit')
    orb = RK[1].data['Orb_Per']
    type1 = RK[1].data['Type1'];type2 = RK[1].data['Type2'];type3 = RK[1].data['Type3']
    M1 = RK[1].data['M1'];M2 = RK[1].data['M2'];M1_M2 = RK[1].data['M1_M2']
    orb = orb * 24;spin = RK[1].data['_3___Per']

    orb_DN=orb[np.where(type1=='DN')]
    orb_Polar=orb[np.where((type2=='AM')|(type3=='AM'))]
    orb_IP=orb[np.union1d(np.where((type2=='IP')|(type2=='DQ')),np.where((type3=='IP')|(type3=='DQ')))]
    spin_IP=spin[np.union1d(np.where((type2=='IP')|(type2=='DQ')),np.where((type3=='IP')|(type3=='DQ')))]
    spin_IP/=3600.
    fig=plt.figure(1,figsize=(9,6))
    ax1=fig.add_subplot(111)
    ax1.plot()
    ax1.hist(orb_Polar,bins=bins,histtype='step',lw=1.5,color='blue',linestyle='--')
    ax1.hist(orb_DN,bins=bins,histtype = 'step',lw=1,color='red',linestyle='-')
    ax1.hist(orb_IP,bins=bins_2,histtype = 'step',lw=1.5,color='green',linestyle='dashdot')
    ax1.hist(spin_IP, bins = bins_spin, histtype = 'step',lw=1.5, color = 'purple',linestyle='dotted')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlim(7e-2,60)
    #print(len(np.where(spin_IP > 0)[0]))
    #plt.legend(['NSC','LW','Polar','DN','IP','Spin of IP'])
    # plt.show()
    return (fig,ax1)

def plot_CV_all(save=0,show=1):
    bins=np.logspace(np.log10(0.5), np.log10(30), 41)
    bins_2 =np.logspace(np.log10(0.8), np.log10(100), 41)
    bins_spin=np.logspace(np.log10(3/36.), np.log10(2), 31)
    bins_p=np.logspace(np.log10(0.5), np.log10(30), 21)

    path_fits = '/Users/baotong/Desktop/period_LW/'
    RK = fits.open(path_fits + 'RK14.fit')
    orb = RK[1].data['Orb_Per']
    type1 = RK[1].data['Type1'];type2 = RK[1].data['Type2'];type3 = RK[1].data['Type3']
    M1 = RK[1].data['M1'];M2 = RK[1].data['M2'];M1_M2 = RK[1].data['M1_M2']
    orb = orb * 24;spin = RK[1].data['_3___Per']
    orb_DN=orb[np.where(type1=='DN')]
    orb_Polar=orb[np.where((type2=='AM')|(type3=='AM'))]
    orb_IP=orb[np.union1d(np.where((type2=='IP')|(type2=='DQ')),np.where((type3=='IP')|(type3=='DQ')))]
    spin_IP=spin[np.union1d(np.where((type2=='IP')|(type2=='DQ')),np.where((type3=='IP')|(type3=='DQ')))]
    spin_IP/=3600.
    orb_all=np.concatenate((orb_DN,orb_IP,orb_Polar))
    orb_all=orb_all[np.where(orb_all>0)]
    print(len(orb_all))
    (ra, dec, seq, period, L, Lmin, Lmax, type)=read_excel('47Tuc')
    period_Tuc=period[np.where(type=='CV')]
    period_Tuc/=3600.

    (ra, dec, seq, period, L, Lmin, Lmax)=load_LW('result_LW')
    period_LW=period[np.where((period>3600)&(period<40000))]
    period_LW/=3600.

    # fig, axes = plt.subplots(2, 1, figsize=(15, 10), sharex='all')
    fig, (ax1, ax2) = plt.subplots(ncols=1, nrows=2, sharex='all',gridspec_kw={'height_ratios': [2, 1]}, figsize=(15, 10))
    # ax1=fig.add_subplot(211)
    ax1.plot()
    ax1.hist(orb_all,bins=bins,histtype='step',lw=2,color='red',linestyle='--')
    ax1.hist(period_Tuc, bins=bins_p, histtype='step', lw=2, linestyle='-',facecolor='c',
         hatch='/', edgecolor='k',fill=True)
    ax1.hist(period_LW, bins=bins_2, histtype='step', lw=3, color='blue', linestyle='-')
    # ax1.hist(spin_IP, bins = bins_spin, histtype = 'step',lw=1.5, color = 'purple',linestyle='-')
    ax1.legend(['Field CVs', 'CVs in 47 Tuc', 'CVs in Galactic Bulge'])

    P_min = 7./6.
    P_gap = [7600.0 / 3600., 11448.0 / 3600.]
    ax1.set_ylim(8e-1, 360)
    ax1.set_xlim(0.5, 60)
    ax1.plot([P_min, P_min], [0, 220], '--',color='grey')
    ax1.plot([P_gap[0], P_gap[0]], [0, 220], '-',lw=2.,color='orange')
    ax1.plot([P_gap[1], P_gap[1]], [0, 220], '-',lw=2.,color='orange')
    ax1.text(P_gap[0]+0.3, 250, 'gap',fontsize=14)
    ax1.text(P_min-0.2, 250, 'minimum',fontsize=14)
    # ax1.set_xlabel('Period (hours)',font1)
    ax1.set_ylabel('Number of sources',font1)
    ax1.tick_params(labelsize=16)
    ax1.set_xscale('log')
    ax1.set_yscale('log')

    # ax2=fig.add_subplot(212)
    ax2.hist(orb_all,bins=bins,histtype='step',lw=2,color='red',cumulative=1,density=1,linestyle='--')
    ax2.hist(period_Tuc, bins=bins_p, histtype='step', lw=2, linestyle='-',cumulative=1,density=1,facecolor='c',
         hatch='/', edgecolor='k',fill=True)
    ax2.hist(period_LW, bins=bins_2, histtype='step', lw=3, color='blue',cumulative=1,density=1, linestyle='-')
    # ax2.hist(spin_IP, bins = bins_spin, histtype = 'step',lw=1.5, color = 'purple',cumulative=1,density=1,linestyle='-')

    ax2.plot([P_min, P_min], [0, 1], '--',color='grey')
    ax2.plot([P_gap[0], P_gap[0]], [0, 1], '-',lw=2.,color='orange')
    ax2.plot([P_gap[1], P_gap[1]], [0, 1], '-',lw=2.,color='orange')

    ax2.set_xscale('log')
    ax2.set_xlabel('Period (hours)',font1)
    ax2.set_ylabel('CDF',font1)
    ax2.tick_params(labelsize=16)
    # ax2.set_yscale('log')
    if save:
        plt.savefig(path_out + '47Tuc_NP.eps', bbox_inches='tight', pad_inches=0.05)
    if show:
        plt.show()

def plot_NP(save=1,show=1):
    (fig,ax1)=plot_RK_CV()
    (ra, dec, seq, period, L, Lmin, Lmax, type)=read_excel('47Tuc')
    period/=3600.
    period=period[np.where(type=='CV')]
    bins_p=np.logspace(np.log10(2.3), np.log10(27), 9)
    print(bins_p)
    ax1.hist(period,bins=bins_p,histtype = 'step',lw=3,color='black')
    ax1.set_xlabel('Period (hours)',font1)
    ax1.set_ylabel('Number of sources',font1)
    plt.tick_params(labelsize=16)
    ax1.legend(['Polar', 'DN', 'IP', 'Spin of IP','CV in 47Tuc'])
    P_min = 7./6.
    P_gap = [7740.0 / 3600., 11448.0 / 3600.]
    ax1.plot([P_min, P_min], [0, 200], '--',color='grey')
    ax1.plot([P_gap[0], P_gap[0]], [0, 200], '-',lw=2.,color='grey')
    ax1.plot([P_gap[1], P_gap[1]], [0, 200], '-',lw=2.,color='grey')

    ax1.text(P_gap[0]+0.1, 230, 'gap',fontsize=12)
    ax1.text(P_min - 0.41, 220, 'minimum',fontsize=12)

    if save:
        plt.savefig(path_out+'47Tuc_NP.eps',bbox_inches='tight', pad_inches=0.0)
    if show:
        plt.show()

    return None

def plot_dist_profile():
    (ra, dec, seq, period, L, Lmin, Lmax, type)=read_excel('47Tuc')
    print(type)
    [ra_center,dec_center,rhl,rc]=pos_all['47Tuc']
    c1 = SkyCoord(ra * u.deg, dec * u.deg, frame='fk5')
    c2 = SkyCoord(ra_center * u.deg, dec_center * u.deg, frame='fk5')
    dist = c1.separation(c2)
    dist = dist.arcmin

    period_CV=period[np.where(type=='CV')]
    dist_CV=dist[np.where(type=='CV')]
    L_CV=L[np.where(type=='CV')]

    period_LB=period[np.where(type=='LMXB')]
    dist_LB=dist[np.where(type=='LMXB')]
    L_LB=L[np.where(type=='LMXB')]
    print(dist_LB)

    period_AB=period[np.where(type=='AB')]
    dist_AB=dist[np.where(type=='AB')]
    L_AB=L[np.where(type=='AB')]

    fig=plt.figure(1,figsize=(12,8))
    ax2=fig.add_subplot(211)
    ax2.scatter(period_CV/3600.,L_CV,marker='v',s=80,color='w',linewidths=2,edgecolors='red')
    ax2.scatter(period_LB/3600.,L_LB,marker='*',s=80,color='w',linewidths=2,edgecolors='green')
    ax2.scatter(period_AB/3600.,L_AB,marker='o',s=80,color='w',linewidths=2,edgecolors='purple')
    ax2.set_yscale('log')
    ax2.set_xscale('log')
    ax2.set_ylabel(r'Luminosity ($\rm erg~s^{-1}$)',font1)
    ax2.tick_params(labelsize=16)
    ax2.legend(['CV','LMXB','AB'],loc='best')
    P_gap = [7740.0 / 3600., 11448.0 / 3600.]
    P_min = [4902.0 / 3600., 4986.0/3600]
    y1 = [0, 2e33]
    # ax2.text(7900 / 3600., 5e36, 'period gap')

    ax2.text(P_gap[0] + 0.25, y1[1], r'$\rm P_{gap}$',fontsize=18)
    ax2.text(P_min[1] - 0.15, y1[1], r'$\rm P_{min}$',fontsize=18)
    ax2.set_ylim(ymax=4e33)
    ax2.fill_between(P_gap, y1[1], facecolor='yellow', alpha=0.2)
    ax2.fill_between(P_min, y1[1], facecolor='grey', alpha=0.2)

    ax3=fig.add_subplot(212)
    ax3.scatter(period_CV/3600.,dist_CV,marker='v',s=80,color='w',linewidths=2,edgecolors='red')
    ax3.scatter(period_LB/3600.,dist_LB,marker='*',s=80,color='w',linewidths=2,edgecolors='green')
    ax3.scatter(period_AB/3600.,dist_AB,marker='o',s=80,color='w',linewidths=2,edgecolors='purple')
    ax3.plot([P_min[0],30],[rhl/60,rhl/60],'--')
    ax3.plot([P_min[0],30],[rc/60,rc/60],'--')
    ax3.text(P_min[0]-0.1,rhl/60,r'$\rm r_h$',fontsize=15)
    ax3.text(P_min[0]-0.1,rc/60,r'$\rm r_c$',fontsize=15)

    ax3.set_yscale('log')
    ax3.set_xscale('log')
    ax3.set_ylabel(r'R ($\rm arcmin$)',font1)
    ax3.set_xlabel('Period (hours)',font1)
    ax3.tick_params(labelsize=16)

    ax2.get_shared_x_axes().join(ax2, ax3)
    y2 = [0, 4]
    # ax2.text(7900 / 3600., 5e36, 'period gap')
    ax3.fill_between(P_min, y2[1], facecolor='grey', alpha=0.2)
    ax3.fill_between(P_gap, y2[1], facecolor='yellow', alpha=0.2)
    plt.savefig(path_out+'47Tuc_profile.pdf',bbox_inches='tight', pad_inches=0.0)
    plt.show()

def plot_erosita_lightcurve(dataname,ecf):
    path_Tuc = f'/Users/baotong/eSASS/data/raw_data/47_Tuc/txt/txt_merge_psf{ecf}_0.2_5/'
    path = path_Tuc
    dataname = '{0}.txt'.format(dataname);epoch_file = path + 'epoch_src_' + dataname
    src_evt=np.loadtxt(path+dataname);epoch_info=np.loadtxt(epoch_file)
    CR=hawk.plot_longT_V(src_evt=src_evt, bkg_file=None,epoch_info=epoch_info)
    CR/=ecf/100.
    (useid, epoch_info_use)=hawk.choose_obs(epoch_info,flux_info=CR,
                                            flux_filter=1e-1,expT_filter=1000,
                                            if_flux_high=False, if_expT_high=True,obsID=[700014])
    src_evt_use =hawk.filter_obs(src_evt, useid)
    time=hawk.filter_energy(src_evt_use[:,0],src_evt_use[:,1],[200,5000])
    bin_len=1000
    lc=hawk.get_hist(time,len_bin=bin_len)

def plot_twopfold():
    period = 31200.27
    net_p = 0.999
    figurepath = '/Users/baotong/Desktop/aas/pXS_Tuc/figure/'
    time1=1;epoch_info_use1=1;
    time2=2;epoch_info_use2=2;
    hawk.phase_fold(time=time,epoch_info=epoch_info_use,net_percent=net_p,p_test=period,outpath=figurepath,bin=100,shift=0.5,label='',save=1,show=1)

if __name__=="__main__":
    plot_CV_all(save=1,show=1)
    # plot_dist_profile()
