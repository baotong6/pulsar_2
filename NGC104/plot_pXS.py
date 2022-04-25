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
from stingray.lightcurve import Lightcurve
import hawkeye as hawk
import rocket as rocket
from scipy import optimize as op
from timing_comb import load_data,get_lc_frombkgimg
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
    print('CV number',len(orb_all))
    print('short period CV',len(orb_all[np.where(orb_all<7700/3600.)]))
    print('long period CV',len(orb_all[np.where(orb_all>11448.0 / 3600)]))
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
    P_gap = [7740.0 / 3600., 11448.0 / 3600.]
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

def plot_dist_profile(save=0,show=1):
    (ra, dec, seq, period, L, Lmin, Lmax, type)=read_excel('47Tuc')
    (ra_AB03, dec_AB03, seq_AB03, period_AB03, L_AB03, Lmin_AB03, Lmax_AB03, type_AB03)=read_excel('47Tuc_AB')
    [ra_center,dec_center,rhl,rc]=pos_all['47Tuc']
    c1 = SkyCoord(ra * u.deg, dec * u.deg, frame='fk5')
    c2 = SkyCoord(ra_center * u.deg, dec_center * u.deg, frame='fk5')
    dist = c1.separation(c2)
    dist = dist.arcmin
    print(seq)
    print(dist*60)
    c3= SkyCoord(ra_AB03 * u.deg, dec_AB03 * u.deg, frame='fk5')
    dist_AB03 = c3.separation(c2)
    dist_AB03 = dist_AB03.arcmin


    period_CV=period[np.where(type=='CV')]
    dist_CV=dist[np.where(type=='CV')]
    L_CV=L[np.where(type=='CV')]

    period_LB=period[np.where(type=='LMXB')]
    dist_LB=dist[np.where(type=='LMXB')]
    L_LB=L[np.where(type=='LMXB')]

    period_AB=period[np.where(type=='AB')]
    dist_AB=dist[np.where(type=='AB')]
    L_AB=L[np.where(type=='AB')]

    period_AB03=period_AB03[np.where(type_AB03=='xAB')]*3600
    dist_AB03=dist_AB03[np.where(type_AB03=='xAB')]
    L_AB03=L_AB03[np.where(type_AB03=='xAB')]
    print(period_AB03)

    fig=plt.figure(1,figsize=(12,8))
    ax2=fig.add_subplot(211)
    ax2.scatter(period_CV/3600.,L_CV,marker='v',s=80,color='w',linewidths=2,edgecolors='red')
    ax2.scatter(period_LB/3600.,L_LB,marker='*',s=80,color='w',linewidths=2,edgecolors='green')
    ax2.scatter(period_AB/3600.,L_AB,marker='o',s=80,color='w',linewidths=2,edgecolors='purple')
    ax2.scatter(14366.89/3600.,4.32e30,marker='^',s=80,color='w',linewidths=2,edgecolors='black')
    ax2.set_yscale('log')
    ax2.set_xscale('log')
    ax2.set_ylabel(r'Luminosity ($\rm erg~s^{-1}$)',font1)
    ax2.tick_params(labelsize=16)
    ax2.legend(['CV','LMXB','AB','Src-No.481'],loc='best')
    P_gap = [7740.0 / 3600., 11048.0 / 3600.]
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
    ax3.scatter(14366.89/3600.,270.86/60,marker='^',s=80,color='w',linewidths=2,edgecolors='black')
    # ax3.scatter(period_AB03/3600.,dist_AB03,marker='o',s=80,color='w',linewidths=2,edgecolors='purple')

    ax3.plot([P_min[0],30],[rhl/60,rhl/60],'--')
    ax3.plot([P_min[0],30],[rc/60,rc/60],'--')
    ax3.plot([P_min[0],30],[1.7,1.7],'-',color='grey')
    ax3.plot([P_min[0],30],[40,40],'-',color='grey')

    ax3.text(P_min[0]-0.1,rhl/60,r'$\rm r_h$',fontsize=15)
    ax3.text(P_min[0]-0.1,rc/60,r'$\rm r_c$',fontsize=15)
    ax3.text(P_min[0] - 0.1, 1.2, r'$\rm R_{in}$', fontsize=15)
    ax3.text(P_min[0]-0.1,30,r'$\rm R_{out}$',fontsize=15)

    ax3.set_yscale('log')
    ax3.set_xscale('log')
    ax3.set_ylabel(r'R ($\rm arcmin$)',font1)
    ax3.set_xlabel('Period (hours)',font1)
    ax3.tick_params(labelsize=16)

    ax2.get_shared_x_axes().join(ax2, ax3)
    y2 = [0, 4]
    # ax2.text(7900 / 3600., 5e36, 'period gap')
    ax3.fill_between(P_min, 40, facecolor='grey', alpha=0.2)
    ax3.fill_between(P_gap, 40, facecolor='yellow', alpha=0.2)
    ax3.fill_between([P_min[0],30],y1=1.7,y2=40,facecolor='blue',alpha=0.1)

    # plt.figure(2)
    # plt.hist(dist_CV*60,bins=20,histtype='step',lw=2,color='red',cumulative=1,density=1,linestyle='--')
    # plt.hist(dist_AB03*60, bins=30, histtype = 'step', lw = 3, color = 'blue', cumulative = 1, density = 1, linestyle = '-')
    # plt.semilogx()
    # plt.show()
    if save:
        plt.savefig(path_out + '47Tuc_profile.pdf', bbox_inches='tight', pad_inches=0.05)
    if show:
        plt.show()

    # plt.figure(2)
    # plt.scatter(dist_CV/rhl,L_CV)
    # plt.show()

# def plot_twopfold():
#     period = 31200.27
#     net_p = 0.999
#     figurepath = '/Users/baotong/Desktop/aas/pXS_Tuc/figure/'
#     time1=1;epoch_info_use1=1;
#     time2=2;epoch_info_use2=2;
#     hawk.phase_fold(time=time,epoch_info=epoch_info_use,net_percent=net_p,p_test=period,outpath=figurepath,bin=100,shift=0.5,label='',save=1,show=1)
def plot_src_lc_singleobs(ifsin=None,figurepath=None,save=0,show=0):
    dataname='481'
    bin_len = 1000
    period=14366.89
    (src_evt_use,epoch_info_use)=load_data(dataname=dataname,ecf=75)
    time = src_evt_use[:, 0]
    lc=hawk.get_hist(time,len_bin=bin_len,tstart=epoch_info_use[:,0][0],tstop=epoch_info_use[:,1][-1])
    lc_net = get_lc_frombkgimg(int(dataname), src_evt_use, epoch_info_use, ecf=75, bin_len=bin_len)
    ## delete gap points ##
    index_list=[];index_net_list=[]
    for i in range(len(epoch_info_use)-1):
        index=np.where((lc.time>epoch_info_use[i][1])&(lc.time<epoch_info_use[i+1][0]))
        index_net=np.where((lc_net.time>epoch_info_use[i][1])&(lc_net.time<epoch_info_use[i+1][0]))
        index=list(index[0]);index_list.extend(index)
        index_net = list(index_net[0]);index_net_list.extend(index_net)

    a=np.delete(lc.time,index_list)
    b=np.delete(lc.counts,index_list)
    c=np.delete(lc_net.time,index_net_list)
    d=np.delete(lc_net.counts,index_net_list)

    fig, (ax1, ax2) = plt.subplots(ncols=1, nrows=2, sharex='all', gridspec_kw={'height_ratios': [1, 1],'hspace':0.1},
                                   figsize=(12, 8))
    # ax1.set_title(r'$T_0={0}$'.format(a[0]),font1)
    x=a-a[0]
    y2=b;y3=d
    y2_err = np.array(poisson_conf_interval(y2, interval='frequentist-confidence'))
    y2_err[0] = y2 - y2_err[0]
    y2_err[1] = y2_err[1] - y2
    y3_err = np.array(poisson_conf_interval(y3, interval='frequentist-confidence'))
    y3_err[0] = y3- y3_err[0]
    y3_err[1] = y3_err[1] - y3
    # plt.figure(1,(9,6))
    if ifsin:
        (popt,perr)=hawk.curvefit_sin(x,y2,0.5*(y2_err[0]+y2_err[1]),period)
        print(popt)
        # y1 = sin_temp(x,popt[0],popt[1],popt[2],popt[3])
        x1=np.linspace(x.min(),x.max(),10000)
        y1 = hawk.sin_temp(x1,period,6.28,8,8)
        y2 = hawk.sin_temp(x1, period*2, 6.28, 8, 8)
        plt.plot(x1,y1)
        plt.plot(y2)
    # absolute_sigma = True, sigma = yerr,

    ### add bkg ###
    obsIDlist=epoch_info_use[:,2].astype('int')
    blank = np.zeros(len(obsIDlist)) + 1
    record=0
    for i in range(len(epoch_info_use)):
        if src_evt_use.ndim < 2 or len(src_evt_use) < 10:
            blank[i] = 0
            continue
        obsid=obsIDlist[i]
        path='/Users/baotong/eSASS/data/raw_data/47_Tuc/txt/txt_psf{0}_{1}/'.format(75,obsid)
        src_info=np.loadtxt(path+'src_info.txt')
        bkg_cts_est_list = src_info[:, 2]
        bkg_cts_est = bkg_cts_est_list[np.where(src_info[:, 0] == int(dataname))][0]
        bkg=bkg_cts_est * bin_len / (epoch_info_use[:,1][i] - epoch_info_use[:,0][i])
        if i==0:
            ax1.plot([epoch_info_use[i][1]-lc.time[0],epoch_info_use[i][1]-lc.time[0]],[0,30],'--',color='grey')
            ax2.plot([epoch_info_use[i][1] - lc.time[0], epoch_info_use[i][1] - lc.time[0]], [0, 20], '--', color='grey')
        elif i==len(epoch_info_use)-1:
            ax1.plot([epoch_info_use[i][0]-lc.time[0],epoch_info_use[i][0]-lc.time[0]],[0,30],'--',color='grey')
            ax2.plot([epoch_info_use[i][0]-lc.time[0],epoch_info_use[i][0]-lc.time[0]],[0,20],'--',color='grey')
        else:
            ax1.plot([epoch_info_use[i][0]-lc.time[0],epoch_info_use[i][0]-lc.time[0]],[0,30],'--',color='grey')
            ax1.plot([epoch_info_use[i][1]-lc.time[0],epoch_info_use[i][1]-lc.time[0]],[0,30],'--',color='grey')
            ax2.plot([epoch_info_use[i][0]-lc.time[0],epoch_info_use[i][0]-lc.time[0]],[0,20],'--',color='grey')
            ax2.plot([epoch_info_use[i][1]-lc.time[0],epoch_info_use[i][1]-lc.time[0]],[0,20],'--',color='grey')
        print(bkg)
        b_1sigma = poisson_conf_interval(bkg, interval='frequentist-confidence').T
        bkg_y_low = b_1sigma[0];bkg_y_high = b_1sigma[1]
        ax1.fill_between(epoch_info_use[i][0:2]-a[0], bkg_y_low, bkg_y_high, facecolor='green', alpha=0.5)

    ax1.errorbar(x, y2,yerr=y2_err, fmt='co', capsize=4, elinewidth=2, ecolor='red',color='green')
    ax1.set_ylabel('Counts/bin',font1)
    ax1.tick_params(labelsize=16)
    x2=c-c[0]
    ax2.errorbar(x2, y3, yerr=y3_err, fmt='co', capsize=4, elinewidth=2, ecolor='red', color='green')
    x1 = np.linspace(x2.min(), x2.max(), 10000)
    y1 = hawk.sin_temp(x1, period, -1.4, 8, 8)
    y2= hawk.sin_temp(x1, period*2,6.28, 8, 8)
    ax2.plot(x1, y1)
    ax2.plot(x1,y2)
    ax2.set_xlabel(r'Time-$T_0$ (second)',font1)
    ax2.set_ylabel('Net-Counts/bin',font1)
    ax2.tick_params(labelsize=16)
    if save:plt.savefig(figurepath+f'{dataname}_lc_4longobs.pdf',bbox_inches='tight', pad_inches=0.01)
    if show:plt.show()
    else:plt.close()

def plot_LSandpfold():
    dataname='481'
    bin_len = 100
    period=14366.89
    (src_evt_use,epoch_info_use)=load_data(dataname=dataname,ecf=75)
    lc_net = get_lc_frombkgimg(int(dataname), src_evt_use, epoch_info_use, ecf=75, bin_len=bin_len)
    time = src_evt_use[:, 0]
    hawk.phase_fold(time=time,epoch_info=epoch_info_use,net_percent=net_p,p_test=period,outpath=figurepath,bin=100,shift=0.83,
                    label=dataname,text='Seq.162 (W58)',save=1,show=1)
    (FP, out_period, max_NormLSP)=hawk.get_LS(lc.time,lc.counts,freq=freq,outpath=figurepath, outname=str(dataname),save=1,show=1)
    return None

def plot_CR_all():
    path='/Users/baotong/Desktop/period_Tuc/'
    catalog=fits.open(path+'xray_properties-592.fits')[1].data
    Seq=catalog['Seq']
    ra=catalog['RAdeg']
    dec=catalog['Dedeg']
    src_cts=np.loadtxt(path+'xray592_srccts_t.txt')
    inter_srcID=rocket.select_src_bypos(Seq,ra,dec,ra_c=6.0236250,dec_c=-72.0812833,inter_radius=2.6*60,outpath=None,outname=None,save=0)
    inter_srcID=inter_srcID-1
    bright_src_index=np.intersect1d(np.where(src_cts>100)[0],inter_srcID)
    print(len(inter_srcID))
    print(len(bright_src_index))
    fig, axes = plt.subplots(1, 1, figsize=(9, 6), sharex='all', sharey='all')
    bins = np.logspace(0.1,4,30)
    bins[14]=101.0
    DIS_CTS = axes.hist(src_cts[inter_srcID], bins=bins, linewidth=3, histtype='step', color='green')
    DIS_bright = axes.hist(src_cts[bright_src_index], bins, histtype='step', color='red', facecolor='r', hatch='/',
                              edgecolor='k', fill=True)
    axes.set_yscale('log')
    axes.set_xscale('log')
    axes.tick_params(labelsize=16)
    plt.show()

def plot_CR_GCLW(save=0,show=1):
    catalog=fits.open(path+'xray_properties-592.fits')[1].data
    Seq=catalog['Seq']
    ra=catalog['RAdeg']
    dec=catalog['Dedeg']
    src_cts_GC=np.loadtxt(path+'xray592_srccts_t.txt')

    cat2=pd.read_excel('/Users/baotong/Desktop/period_LW/catalog_LW.xlsx')
    net_cts=cat2['net_cts'];bkg_cts=cat2['bkg_cts']
    src_cts_LW=net_cts+bkg_cts
    src_cts_LW=np.array(src_cts_LW)
    bright_src_GC=src_cts_GC[np.where(src_cts_GC>100)]
    bright_src_LW=src_cts_LW[np.where(src_cts_LW>100)]
    bins=np.logspace(2,4,30)
    plt.figure(1,(9,6))
    plt.hist(bright_src_GC,bins=bins,histtype='step',lw=2,cumulative=1,density=1, linestyle='-',color='red')
    plt.hist(bright_src_LW,bins=bins,histtype='step',lw=2,cumulative=1,density=1, linestyle='-',color='green')
    plt.legend(['47 Tuc','LW'],loc='upper left')
    plt.semilogx()
    plt.xlabel('Counts',font1)
    plt.ylabel('CDF',font1)
    plt.tick_params(labelsize=16)
    if save:plt.savefig(path_out+'CDF_GC_LW.pdf',bbox_inches='tight', pad_inches=0.01)
    if show:plt.show()
    else:plt.close()
    return None
if __name__=="__main__":
    # plot_CV_all(save=1,show=1)
    # plot_dist_profile(save=1,show=1)
    plot_src_lc_singleobs(figurepath=path_out,save=0,show=1)
    # plot_CR_all()
    # plot_CR_GCLW(save=1,show=1)