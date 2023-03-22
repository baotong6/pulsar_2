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
    L_1_8=np.array(res['L_1_8'])
    Lmin = np.array(res['Lmin'])
    Lmax = np.array(res['Lmax'])
    type = np.array(res['type'])

    return (ra,dec,seq,period,L,Lmin,Lmax,type,L_1_8)

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
    print('Too long period CV',len(orb_all[np.where(orb_all>12.)]))
    print('CV number',len(orb_all))
    print('short period CV',len(orb_all[np.where(orb_all<7700/3600.)]))
    print('long period CV',len(orb_all[np.where(orb_all>11448.0 / 3600)]))
    (ra, dec, seq, period, L, Lmin, Lmax, type,L_1_8)=read_excel('47Tuc')
    period_Tuc=period[np.where(type=='CV')]
    period_Tuc/=3600.

    (ra, dec, seq, period, L, Lmin, Lmax,type)=load_LW('result_LW')
    period_LW=period[np.where((period>3900)&(period<40000))]
    period_LW/=3600.
    print(period_LW)
    # fig, axes = plt.subplots(2, 1, figsize=(15, 10), sharex='all')
    fig, (ax1, ax2) = plt.subplots(ncols=1, nrows=2, sharex='all',gridspec_kw={'height_ratios': [2, 1]}, figsize=(15, 10))
    # ax1=fig.add_subplot(211)
    ax1.plot()
    ax1.hist(orb_all,bins=bins,histtype='step',lw=2,color='red',linestyle='--')
    ax1.hist(period_Tuc, bins=bins_p, histtype='step', lw=2, linestyle='-',facecolor='c',
         hatch='/', edgecolor='k',fill=True)
    ax1.hist(period_LW, bins=bins_2, histtype='step', lw=3, color='blue', linestyle='-')
    # ax1.hist(spin_IP, bins = bins_spin, histtype = 'step',lw=1.5, color = 'purple',linestyle='-')
    ax1.legend(['CVs in Solar Neighborhood', 'CVs in 47 Tuc', 'CVs in Galactic Bulge'])

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
        plt.savefig(path_out + '47Tuc_NP.pdf', bbox_inches='tight', pad_inches=0.05)
    if show:
        plt.show()

def plot_NP(save=1,show=1):
    (fig,ax1)=plot_RK_CV()
    (ra, dec, seq, period, L, Lmin, Lmax, type,L_1_8)=read_excel('47Tuc')
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
    (ra, dec, seq, period, L, Lmin, Lmax, type,L_1_8)=read_excel('47Tuc')
    print(L)
    (ra_AB03, dec_AB03, seq_AB03, period_AB03, L_AB03, Lmin_AB03, Lmax_AB03, type_AB03,L_1_8_AB)=read_excel('47Tuc_AB')
    [ra_center,dec_center,rhl,rc]=pos_all['47Tuc']
    c1 = SkyCoord(ra * u.deg, dec * u.deg, frame='fk5')
    c2 = SkyCoord(ra_center * u.deg, dec_center * u.deg, frame='fk5')
    dist = c1.separation(c2)
    dist = dist.arcmin
    # print(seq)
    # print(dist*60)
    c3= SkyCoord(ra_AB03 * u.deg, dec_AB03 * u.deg, frame='fk5')
    dist_AB03 = c3.separation(c2)
    dist_AB03 = dist_AB03.arcmin

    # L=L_1_8
    ##换成1-8keV的光度

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
    # print(period_AB03)


    (ra, dec, seq, period, L, Lmin, Lmax,type)=load_LW('result_LW')
    # print(type)
    period_LW_CV=period[np.where(type=='CV')]
    L_LW_CV=L[np.where(type=='CV')]
    L_LW_CV=L_LW_CV*1.11423


    fig=plt.figure(1,figsize=(15,10))
    ax2=fig.add_subplot(211)
    ax2.scatter(period_CV/3600.,L_CV,marker='v',s=80,color='w',linewidths=2,edgecolors='red',label='CV')
    ax2.scatter(period_LB/3600.,L_LB,marker='*',s=80,color='w',linewidths=2,edgecolors='green',label='LMXB')
    ax2.scatter(period_AB/3600.,L_AB,marker='o',s=80,color='w',linewidths=2,edgecolors='purple',label='AB')
    # ax2.scatter(14366.89/3600.,4.32236E+30,marker='^',s=80,color='w',linewidths=2,edgecolors='c',label='Src-No.481')
    ax2.scatter(period_LW_CV/3600.,L_LW_CV*1e31,marker='v',s=80,color='w',linewidths=2,edgecolors='black',label='CV in Galactic Bulge')

    ax2.set_yscale('log')
    ax2.set_xscale('log')
    ax2.set_ylabel(r'Luminosity ($\rm erg~s^{-1}$)',font1)
    ax2.tick_params(labelsize=16)
    ax2.legend()
    # ax2.legend(['CV','LMXB','AB','Src-No.481','CV in Galactic Bulge'],loc='best')
    P_gap = [7740.0 / 3600., 11048.0 / 3600.]
    P_min = [4902.0 / 3600., 4986.0/3600]
    y1 = [0, 2e33]
    # ax2.text(7900 / 3600., 5e36, 'period gap')

    ax2.text(P_gap[0] + 0.25, y1[1], r'$\rm P_{gap}$',fontsize=18)
    ax2.text(P_min[1] - 0.15, y1[1], r'$\rm P_{min}$',fontsize=18)
    ax2.set_ylim(ymax=4e33)
    ax2.set_xlim(xmin=0.95,xmax=35)
    ax2.fill_between(P_gap, y1[1], facecolor='yellow', alpha=0.2)
    ax2.fill_between(P_min, y1[1], facecolor='grey', alpha=0.2)
    ax2.set_xticks([1,10,20,30])
    ax2.set_xticklabels(['1','10','20','30'])

    ax3=fig.add_subplot(212)
    ax3.scatter(period_CV/3600.,dist_CV,marker='v',s=80,color='w',linewidths=2,edgecolors='red',label='CV')
    ax3.scatter(period_LB/3600.,dist_LB,marker='*',s=80,color='w',linewidths=2,edgecolors='green',label='LMXB')
    ax3.scatter(period_AB/3600.,dist_AB,marker='o',s=80,color='w',linewidths=2,edgecolors='purple',label='AB')
    # ax3.scatter(14366.89/3600.,270.86/60,marker='^',s=80,color='w',linewidths=2,edgecolors='c',label='Src-No.481')
    # ax3.scatter(period_AB03/3600.,dist_AB03,marker='o',s=80,color='w',linewidths=2,edgecolors='purple')

    ax3.plot([P_min[0],30],[rhl/60,rhl/60],'--')
    ax3.plot([P_min[0],30],[rc/60,rc/60],'--')
    ax3.plot([P_min[0],30],[1.7,1.7],'-',color='grey')
    ax3.plot([P_min[0],30],[40,40],'-',color='grey')

    ax3.text(P_min[0]-0.1,rhl/60,r'$\rm r_h$',fontsize=15)
    ax3.text(P_min[0]-0.1,rc/60,r'$\rm r_c$',fontsize=15)
    ax3.text(P_min[0] - 0.1, 1.2, r'$\rm R_{in}$', fontsize=15)
    ax3.text(P_min[0]-0.15,30,r'$\rm R_{out}$',fontsize=15)

    ax3.set_yscale('log')
    ax3.set_xscale('log')
    ax3.set_ylabel(r'R ($\rm arcmin$)',font1)
    ax3.set_xlabel('Period (hours)',font1)
    ax3.tick_params(labelsize=16)

    # ax2.get_shared_x_axes().join(ax2, ax3)
    y2 = [0, 4]
    # ax2.text(7900 / 3600., 5e36, 'period gap')
    ax3.fill_between(P_min, 40, facecolor='grey', alpha=0.2)
    ax3.fill_between(P_gap, 40, facecolor='yellow', alpha=0.2)
    ax3.fill_between([P_min[0],30],y1=1.7,y2=40,facecolor='blue',alpha=0.1)
    ax3.set_xlim(xmin=0.95,xmax=35)
    ax3.set_xticks([1,10,20,30])
    ax3.set_xticklabels(['1','10','20','30'])
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
        # plt.plot(x1,y2)
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
    y1 = hawk.sin_temp(x1, period, 6.28, 7, 7)
    y2= hawk.sin_temp(x1, period*2,6.28, 8, 8)
    ax2.plot(x1, y1)
    # ax2.plot(x1,y2)
    ax2.set_xlabel(r'Time (second)',font1)
    ax2.set_ylabel('Net-Counts/bin',font1)
    ax2.tick_params(labelsize=16)
    if save:plt.savefig(figurepath+f'{dataname}_lc_4longobs.pdf',bbox_inches='tight', pad_inches=0.01)
    if show:plt.show()
    else:plt.close()

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

def plot_CV_Temp(save=0,show=1):
    (ra, dec, seq, period, L, Lmin, Lmax, type, L_1_8) = read_excel('47Tuc')
    (ra_AB03, dec_AB03, seq_AB03, period_AB03, L_AB03, Lmin_AB03, Lmax_AB03, type_AB03, L_1_8_AB) = read_excel(
        '47Tuc_AB')
    [ra_center, dec_center, rhl, rc] = pos_all['47Tuc']
    c1 = SkyCoord(ra * u.deg, dec * u.deg, frame='fk5')
    c2 = SkyCoord(ra_center * u.deg, dec_center * u.deg, frame='fk5')
    dist = c1.separation(c2)
    dist = dist.arcmin
    print(seq)
    print(dist * 60)
    c3 = SkyCoord(ra_AB03 * u.deg, dec_AB03 * u.deg, frame='fk5')
    dist_AB03 = c3.separation(c2)
    dist_AB03 = dist_AB03.arcmin

    L = L_1_8
    ##换成1-8keV的光度
    period_CV = period[np.where(type == 'CV')]
    dist_CV = dist[np.where(type == 'CV')]
    L_CV = L[np.where(type == 'CV')]

    period_LB = period[np.where(type == 'LMXB')]
    dist_LB = dist[np.where(type == 'LMXB')]
    L_LB = L[np.where(type == 'LMXB')]

    period_AB = period[np.where(type == 'AB')]
    dist_AB = dist[np.where(type == 'AB')]
    L_AB = L[np.where(type == 'AB')]

    period_AB03 = period_AB03[np.where(type_AB03 == 'xAB')] * 3600
    dist_AB03 = dist_AB03[np.where(type_AB03 == 'xAB')]
    L_AB03 = L_AB03[np.where(type_AB03 == 'xAB')]
    # print(period_AB03)

    (ra, dec, seq, period, L, Lmin, Lmax, type) = load_LW('result_LW')
    print(type)
    period_LW_CV = period[np.where(type == 'CV')]
    L_LW_CV = L[np.where(type == 'CV')]

    CV_T_LW=np.array([40,40,11,40,7,42,40,21,
                      40,40,40,40,40,40,40,40,40,
                      4,40,40])
    low_T_LW=np.array([40,40,6,40,4,10,40,5.8,
                       40,40,40,40,40,40,40,40,40,
                       1,40,40])
    high_T_LW=np.array([50,50,64,50,19,85,50,77,
                        50,50,50,50,50,50,50,50,50,
                        32,50,50,])

    CV_T_Tuc=np.array([7,9,40,4.6,13,20,1.6,19,4.4,2.1,11])
    low_T_Tuc=np.array([6,6,40,3.8,9,12,1.2,10,4.0,1.8,9])
    high_T_Tuc=np.array([9,14,50,5.8,23,45,2.1,54,5.4,2.4,14])

    ## plot fig1 ###
    '''
    fig, (ax1, ax2) = plt.subplots(ncols=1, nrows=2, sharex='all',gridspec_kw={'height_ratios': [1, 1]}, figsize=(6, 6))
    for i in range(len(period_LW_CV)):
        if CV_T_LW[i]==40:
            ax1.errorbar(period_LW_CV[i]/3600, CV_T_LW[i],yerr=np.row_stack((CV_T_LW[i]-low_T_LW[i],high_T_LW[i]-CV_T_LW[i])),
                         fmt='co', uplims=False, lolims=True,
                         capsize=4, elinewidth=2, ecolor='k',color='k')
        else:
            ax1.errorbar(period_LW_CV[i]/3600, CV_T_LW[i],yerr=np.row_stack((CV_T_LW[i]-low_T_LW[i],high_T_LW[i]-CV_T_LW[i])), fmt='co',
                         capsize=4, elinewidth=2, ecolor='k',color='k')
    ax1.errorbar(period_CV/3600, CV_T_Tuc, yerr=np.row_stack((CV_T_Tuc-low_T_Tuc, high_T_Tuc-CV_T_Tuc)), fmt='co', capsize=4, elinewidth=2,
                 ecolor='r', color='r')
    ax1.set_xscale('log')
    ax1.set_yscale('log')

    # ax1.set_xlim(1,15)
    index1=np.where(period_LW_CV<2*3600)
    index2=np.where((period_LW_CV<10*3600)&(period_LW_CV>2*3600))
    index3=np.where(period_LW_CV>10*3600)

    ind1=np.where(period_CV<2*3600)
    ind2=np.where((period_CV<10*3600)&(period_CV>2*3600))
    ind3=np.where(period_CV>10*3600)
    ax2.errorbar(1.5, CV_T_LW[index1].mean(),xerr=np.row_stack([0.5,0.5]),
                 yerr=np.row_stack(((CV_T_LW[index1]-low_T_LW[index1]).mean(),(high_T_LW[index1]-CV_T_LW[index1]).mean())),
                 fmt='co', capsize=4, elinewidth=2, ecolor='k',color='k',label='CVs in Galactic Bulge')
    ax2.errorbar(6, CV_T_LW[index2].mean(),xerr=np.row_stack([4,4]),
                 yerr=np.row_stack(((CV_T_LW[index2]-low_T_LW[index2]).mean(),(high_T_LW[index2]-CV_T_LW[index2]).mean())),
                 fmt='co', capsize=4, elinewidth=2, ecolor='k',color='k')
    ax2.errorbar(12, CV_T_LW[index3].mean(),xerr=np.row_stack([2,2]),
                 yerr=np.row_stack(((CV_T_LW[index3]-low_T_LW[index3]).mean(),(high_T_LW[index3]-CV_T_LW[index3]).mean())),
                 fmt='co', capsize=4, elinewidth=2, ecolor='k',color='k')
    ax2.errorbar(1, CV_T_Tuc[ind1].mean(),xerr=np.row_stack([0.5,0.5]),yerr=np.row_stack(((CV_T_Tuc[ind1]-low_T_Tuc[ind1]).mean(),(high_T_Tuc[ind1]-CV_T_Tuc[ind1]).mean())),
                 fmt='co', capsize=4, elinewidth=2, ecolor='r',color='r',label='CVs in 47 Tuc')
    ax2.errorbar(6, CV_T_Tuc[ind2].mean(),xerr=np.row_stack([4,4]),yerr=np.row_stack(((CV_T_Tuc[ind2]-low_T_Tuc[ind2]).mean(),(high_T_Tuc[ind2]-CV_T_Tuc[ind2]).mean())),
                 fmt='co', capsize=4, elinewidth=2, ecolor='r',color='r')
    ax2.errorbar(12, CV_T_Tuc[ind3].mean(),xerr=np.row_stack([2,2]),yerr=np.row_stack(((CV_T_Tuc[ind3]-low_T_Tuc[ind3]).mean(),(high_T_Tuc[ind3]-CV_T_Tuc[ind3]).mean())),
                 fmt='co', capsize=4, elinewidth=2, ecolor='r',color='r')
    ax2.legend()
    # ax2.set_xlim(1,15)
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax1.set_ylabel(r'$\rm T_b$ (keV)',font1)
    # ax1.set_xlabel('Period (hours)',font1)
    ax1.tick_params(labelsize=16)
    ax2.set_ylabel(r'$\rm T_b$ (keV)',font1)
    ax2.set_xlabel('Period (hours)',font1)
    ax2.tick_params(labelsize=16)
    '''

    ## plot fig2 ###
    (fig,ax1)=plt.subplots(ncols=1, nrows=1, figsize=(9,6))
    for i in range(len(period_LW_CV)):
        if CV_T_LW[i]==40:
            ax1.errorbar(period_LW_CV[i]/3600, CV_T_LW[i],yerr=np.row_stack((CV_T_LW[i]-low_T_LW[i],high_T_LW[i]-CV_T_LW[i])),
                         fmt='ks', uplims=False, lolims=True,
                         capsize=2, elinewidth=0.5, capthick=0.5,ecolor='k',color='k',markersize=3,markerfacecolor='none',markeredgecolor='grey')
        else:
            ax1.errorbar(period_LW_CV[i]/3600, CV_T_LW[i],yerr=np.row_stack((CV_T_LW[i]-low_T_LW[i],high_T_LW[i]-CV_T_LW[i])), fmt='ks',
                         capsize=2, elinewidth=0.5, capthick=0.5,ecolor='k',color='k',markersize=3,markerfacecolor='none',markeredgecolor='grey')
    ax1.errorbar(period_CV/3600, CV_T_Tuc, yerr=np.row_stack((CV_T_Tuc-low_T_Tuc, high_T_Tuc-CV_T_Tuc)), fmt='rs', capsize=2, elinewidth=0.5,
                 ecolor='r', color='r',capthick=0.5,markersize=3,markerfacecolor='none',markeredgecolor='r')

    index1 = np.where(period_LW_CV < 2 * 3600)
    index2 = np.where((period_LW_CV < 10 * 3600) & (period_LW_CV > 2 * 3600))
    index3 = np.where(period_LW_CV > 10 * 3600)

    ind1 = np.where(period_CV < 2 * 3600)
    ind2 = np.where((period_CV < 10 * 3600) & (period_CV > 2 * 3600))
    ind3 = np.where(period_CV > 10 * 3600)
    ax1.errorbar(1.5, CV_T_LW[index1].mean(), xerr=np.row_stack([0.5, 0.5]),
                 yerr=np.row_stack(
                     ((CV_T_LW[index1] - low_T_LW[index1]).mean(), (high_T_LW[index1] - CV_T_LW[index1]).mean())),
                 fmt='ks', capsize=4, elinewidth=2, ecolor='k', color='k', markersize=12,label='CVs in Galactic Bulge')
    ax1.errorbar(6, CV_T_LW[index2].mean(), xerr=np.row_stack([4, 4]),
                 yerr=np.row_stack(
                     ((CV_T_LW[index2] - low_T_LW[index2]).mean(), (high_T_LW[index2] - CV_T_LW[index2]).mean())),
                 fmt='ks', capsize=4, elinewidth=2, ecolor='k', color='k',markersize=12)
    ax1.errorbar(12, CV_T_LW[index3].mean(), xerr=np.row_stack([2, 2]),
                 yerr=np.row_stack(
                     ((CV_T_LW[index3] - low_T_LW[index3]).mean(), (high_T_LW[index3] - CV_T_LW[index3]).mean())),
                 fmt='ks', capsize=4, elinewidth=2, ecolor='k', color='k',markersize=12)
    ax1.errorbar(1, CV_T_Tuc[ind1].mean(), xerr=np.row_stack([0.5, 0.5]), yerr=np.row_stack(
        ((CV_T_Tuc[ind1] - low_T_Tuc[ind1]).mean(), (high_T_Tuc[ind1] - CV_T_Tuc[ind1]).mean())),
                 fmt='rs', capsize=4, elinewidth=2, ecolor='r', color='r',markersize=12, label='CVs in 47 Tuc')
    ax1.errorbar(6, CV_T_Tuc[ind2].mean(), xerr=np.row_stack([4, 4]), yerr=np.row_stack(
        ((CV_T_Tuc[ind2] - low_T_Tuc[ind2]).mean(), (high_T_Tuc[ind2] - CV_T_Tuc[ind2]).mean())),
                 fmt='rs', capsize=4, elinewidth=2, ecolor='r', color='r',markersize=12)
    ax1.errorbar(20, CV_T_Tuc[ind3].mean(), xerr=np.row_stack([10, 10]), yerr=np.row_stack(
        ((CV_T_Tuc[ind3] - low_T_Tuc[ind3]).mean(), (high_T_Tuc[ind3] - CV_T_Tuc[ind3]).mean())),
                 fmt='rs', capsize=4, elinewidth=2, ecolor='r', color='r',markersize=12)
    ax1.legend()
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_ylabel(r'$\rm T_b$ (keV)',font1)
    ax1.set_xlabel('Period (hours)',font1)
    ax1.tick_params(labelsize=16)

    if save:
        plt.savefig(path_out + 'CV_Temp.pdf', bbox_inches='tight', pad_inches=0.05)
    if show:
        plt.show()

def adapt_bin(dist,cxb):
    # criteria net>10 and SN>3
    bin_step=0.01
    bin_lf=0;bin_rt=bin_lf+bin_step
    bin_rt_list=[]
    while bin_rt<np.max(dist):
        temp_cts=len(np.where((dist<bin_rt) & (dist>bin_lf))[0])
        temp_area=(bin_rt**2-bin_lf**2)*np.pi
        net_cts=temp_cts-cxb*temp_area
        if net_cts<0 or temp_cts==0:
            bin_rt+=bin_step
            continue
        else:
            SN=net_cts/np.sqrt(temp_cts)
            if net_cts>5 or SN >3:
                bin_rt_list.append(bin_rt)
                bin_lf=bin_rt;bin_rt+=bin_step
                continue
            else:
                bin_rt+=bin_step
                continue

    return bin_rt_list

def f(x,a,b):
    logS=a*x+b  #S~V^a
    return logS
def spectrafit(x,y,error):
    # popt, pcov = op.curve_fit(f, np.log(x), np.log(y),absolute_sigma=True,sigma=np.log(error))
    popt, pcov = op.curve_fit(f, np.log(x), np.log(y))
    perr = np.sqrt(np.diag(pcov))
    logydata=f(np.log(x),popt[0],popt[1])
    ydata=np.exp(logydata)
    return (popt,perr)

def plot_surface_density(save=0,show=1):
    cat=fits.open('/Users/baotong/Desktop/period_Tuc/cheng2019_Tuc.fit')
    ra=cat[1].data['RAJ2000'];dec=cat[1].data['DEJ2000'];F058=cat[1].data['F0_5-8'];L058=cat[1].data['L0_5-8']
    ## filtering periodic sources ##
    filter_index=np.array([217,414,185,366,423,232,263,331,273,317,162,252,290,198,312,229])-1
    ra=np.delete(ra,filter_index);dec=np.delete(dec,filter_index);F058=np.delete(F058,filter_index);L058=np.delete(L058,filter_index)
    ##--------------------------##
    [ra_center,dec_center,rhl,rc]=pos_all['47Tuc']
    c1 = SkyCoord(ra * u.deg, dec * u.deg, frame='fk5')
    c2 = SkyCoord(ra_center * u.deg, dec_center * u.deg, frame='fk5')
    dist_all = c1.separation(c2)
    dist_all=dist_all.arcmin

    # dist=np.array([8.47522594,49.89957386,23.63156505,35.67687066,58.05747122,3.39603853,3.94875538,15.18944144,
    #                10.12621419,10.80148101,21.69163174,3.32388461,20.65438152,53.67657374,16.24158616,
    #                14.81336603,4.48256703,270.8604792])
    dist=np.array([8.47522594,49.89957386,23.63156505,35.67687066,58.05747122,3.39603853,3.94875538,15.18944144,
                   10.12621419,10.80148101,21.69163174,3.32388461,53.67657374,
                   14.81336603,4.48256703,270.8604792])
    # dist=np.array([49.89957386,23.63156505,35.67687066,58.05747122,15.18944144,
    #                10.12621419,10.80148101,3.32388461,53.67657374,
    #                14.81336603,4.48256703])
    dist=dist/60
    cheng_profile=pd.read_excel(path+'profile.xlsx')
    x_f=cheng_profile.iloc[0:16,0];x_f_err=cheng_profile.iloc[0:16,1]
    y_f=cheng_profile.iloc[0:16,2];y_f_err=cheng_profile.iloc[0:16,3]
    x_f_cxb=cheng_profile.iloc[0:16,4];y_f_cxb=cheng_profile.iloc[0:51,5]
    # y1_net=y1
    # for i in range(len(x1)):
    #     index=np.argmin(x1_cxb-x1[i])
    #     y1_net[i]-=y1cxb[index]
    # x2=cheng_profile.iloc[0:13,15];x2err=cheng_profile.iloc[0:13,16]
    # y2=cheng_profile.iloc[0:13,17];y2err=cheng_profile.iloc[0:13,18]
    # x2_cxb=x2;y2cxb=np.zeros(len(x2))+0.25
    # y2_net=y2-y2cxb

    index_faint=np.where(F058<1e-6);index_bright=np.where(F058>1e-6)
    # index_faint = np.where(L058 < 5e30);index_bright = np.where(L058 > 5e30)
    print('bright group:',len(index_bright[0]))
    print('faint group:',len(index_faint[0]))
    bin_rt_list=adapt_bin(dist_all[index_bright],0.25)
    bin_rt_list=np.concatenate(([0],bin_rt_list))
    bin_rt_list=bin_rt_list[:-1];bin_rt_list[-1]=3.17  ## convert the hist range to FoV
    # print('bin_rt_list',bin_rt_list)
    x1=[(bin_rt_list[i]+bin_rt_list[i+1])/2 for i in range(len(bin_rt_list)-1)]
    x1err = [(bin_rt_list[i + 1] - bin_rt_list[i]) / 2 for i in range(len(bin_rt_list) - 1)]
    area1 = [(bin_rt_list[i + 1] ** 2 - bin_rt_list[i] ** 2) * np.pi for i in range(len(bin_rt_list) - 1)]
    histfaint=plt.hist(dist_all[index_faint],bin_rt_list)
    plt.close()
    histbright=plt.hist(dist_all[index_bright],bin_rt_list)
    plt.close()
    y1=histfaint[0];y2=histbright[0]
    y1_err = np.array(poisson_conf_interval(y1, interval='frequentist-confidence'))
    y1_err[0] = y1 - y1_err[0];y1_err[1] = y1_err[1] - y1
    y2_err = np.array(poisson_conf_interval(y2, interval='frequentist-confidence'))
    y2_err[0] = y2 - y2_err[0];y2_err[1] = y2_err[1] - y2
    print('y2=',y2)
    y1=y1/area1
    y1_net=y1
    for i in range(len(x1)-1):
        index=np.argmin(x_f_cxb-x1[i+1])
        y1_net[i]-=y_f_cxb[index]
    print(y1_net)
    y2_net=y2/area1-0.25
    y1_err=y1_err/area1;y2_err=y2_err/area1
    print('bins_group:',bin_rt_list)
    bins=[0,0.1,0.2,0.36,1.0]
    hist1=plt.hist(dist,bins)
    plt.close()
    num1=hist1[0]
    y=num1
    x=[(bins[i]+bins[i+1])/2 for i in range(len(bins)-1)]
    area=[(bins[i+1]**2-bins[i]**2)*np.pi for i in range(len(bins)-1)]
    y_err = np.array(poisson_conf_interval(y, interval='frequentist-confidence'))
    y_err[0] = y - y_err[0]
    y_err[1] = y_err[1] - y
    print(y_err)
    y=y/area
    y_err=y_err/area
    xerr=[(bins[i+1]-bins[i])/2 for i in range(len(bins)-1)]
    (popt, perr) = spectrafit(x1[0:3], y1_net[0:3], y1_err[0][0:3])
    print('Faint group')
    print('popt=',popt)
    print('perr=',perr)
    x1=np.concatenate(([1e-2],x1))
    plt.figure(1, (8, 8))
    plt.plot(x1[0:4], np.exp(f(np.log(x1[0:4]), popt[0], popt[1])),'-.',color='r')
    # plt.fill_between(x1[0:7],np.exp(f(np.log(x1[0:7]), popt[0]-perr[0], popt[1]-perr[1])),
    #                  np.exp(f(np.log(x1[0:7]), popt[0]+perr[0], popt[1]+perr[1])),facecolor = 'r', alpha = 0.4)
    x1 = x1[1:]
    (popt, perr) = spectrafit(x1[0:3], y2_net[0:3], y2_err[0][0:3])
    print('Bright group')
    print('popt=',popt)
    print('perr=',perr)

    x1=np.concatenate(([1e-2],x1))
    plt.plot(x1[0:4], np.exp(f(np.log(x1[0:4]), popt[0], popt[1])),'-.',color='g')
    # plt.fill_between(x1[0:7],np.exp(f(np.log(x1[0:7]), popt[0]-perr[0], popt[1]-perr[1])),
    #                  np.exp(f(np.log(x1[0:7]), popt[0]+perr[0], popt[1]+perr[1])),facecolor = 'g', alpha = 0.4)
    # x1 = x1[1:]
    (popt, perr) = spectrafit(x[0:3], y[0:3], y_err[0][0:3])
    x = np.concatenate(([1e-2], x))
    plt.plot(x[0:4], np.exp(f(np.log(x[0:4]), popt[0], popt[1])),'-.',linewidth=2,color='k')

    print('Periodic group')
    print('popt=',popt)
    print('perr=',perr)
    # plt.fill_between(x[0:4],np.exp(f(np.log(x[0:4]), popt[0]-perr[0], popt[1]-perr[1])),
    #                  np.exp(f(np.log(x[0:4]), popt[0]+perr[0], popt[1]+perr[1])),facecolor = 'k', alpha = 0.4)
    x1 = x1[1:];x=x[1:]
    plt.errorbar(x1,y1_net,xerr=x1err,yerr=y1_err,fmt='ro', capsize=1, elinewidth=1, ecolor='r', color='r',markersize=4,label='Faint group')
    plt.errorbar(x1,y2_net,xerr=x1err,yerr=y2_err,fmt='go', capsize=1, elinewidth=1, ecolor='g', color='g',markersize=4,label='Bright group')
    plt.errorbar(x,y,xerr=xerr,yerr=y_err,fmt='ks', capsize=2, elinewidth=2, ecolor='k', color='k',markersize=8,label=r'$\rm Periodic~CVs$')
    plt.plot([3.17,3.17],[0,1000],'--')
    plt.plot([3.17/8.8,3.17/8.8],[0,1000],'--')
    plt.text(3.17/8.8,1400,r'$r_c$',font1,horizontalalignment='center',verticalalignment='center')
    plt.text(3.17,1400,r'$r_h$',font1,horizontalalignment='center',verticalalignment='center')
    plt.legend()
    plt.semilogy()
    plt.semilogx()
    plt.xlim(1e-2,7)
    plt.ylim(1e-2,5e3)
    plt.xlabel('R (arcmin)',font1)
    plt.ylabel(r'$\rm Number~of~source~per~arcmin^2$',font1)
    plt.tick_params(labelsize=16)
    if save:
        plt.savefig(path_out + 'profile_3group.pdf', bbox_inches='tight', pad_inches=0.05)
    if show:
        plt.show()
    '''
    plt.errorbar(x1,y1_net,xerr=x1err,yerr=y1err,fmt='ro', capsize=2, elinewidth=2, ecolor='r', color='r',markersize=4,label='Faint group')
    plt.errorbar(x2,y2_net,xerr=x2err,yerr=y2err,fmt='go', capsize=2, elinewidth=2, ecolor='g', color='g',markersize=4,label='Bright group')
    plt.errorbar(x,y,xerr=xerr,yerr=y_err,fmt='ks', capsize=2, elinewidth=2, ecolor='k', color='k',markersize=4,label='Periodic binaries')
    plt.legend()
    plt.semilogy()
    plt.xlabel('R (arcmin)',font1)
    plt.ylabel('Number of source per arcmin2',font1)
    plt.tick_params(labelsize=16)
    plt.show()

    bins1=np.array([0.01,0.08,0.2,0.3,0.4,0.5,0.7,1.0,2.0,3.0,4.0,5.0,6.0,7.0])
    bins1=np.linspace(0.01,7.0,12)
    x1=[(bins1[i]+bins1[i+1])/2 for i in range(len(bins1)-1)]
    x1err = [(bins1[i + 1] - bins1[i]) / 2 for i in range(len(bins1) - 1)]
    area1 = [(bins1[i + 1] ** 2 - bins1[i] ** 2) * np.pi for i in range(len(bins1) - 1)]
    histfaint=plt.hist(dist_all[index_faint],bins1)
    plt.close()
    histbright=plt.hist(dist_all[index_bright],bins1)
    plt.close()
    y1=histfaint[0];y2=histbright[0]
    y1_err = np.array(poisson_conf_interval(y1, interval='frequentist-confidence'))
    y1_err[0] = y1 - y1_err[0];y1_err[1] = y1_err[1] - y1
    y2_err = np.array(poisson_conf_interval(y2, interval='frequentist-confidence'))
    y2_err[0] = y2 - y2_err[0];y2_err[1] = y2_err[1] - y2
    y1=y1/area1;y2=y2/area1
    y1_err=y1_err/area1;y2_err=y2_err/area1
    print('y2=',y2)

    plt.errorbar(x1,y1,xerr=x1err,yerr=y1_err,fmt='ro', capsize=2, elinewidth=2, ecolor='r', color='r',markersize=4,label='Faint group')
    plt.errorbar(x1,y2,xerr=x1err,yerr=y2_err,fmt='go', capsize=2, elinewidth=2, ecolor='g', color='g',markersize=4,label='Bright group')
    plt.errorbar(x,y,xerr=xerr,yerr=y_err,fmt='ks', capsize=2, elinewidth=2, ecolor='k', color='k',markersize=4,label='Periodic binaries')
    plt.legend()
    plt.semilogy()
    plt.xlabel('R (arcmin)',font1)
    plt.ylabel('Number of source per arcmin2',font1)
    plt.tick_params(labelsize=16)
    plt.show()
    '''

if __name__=="__main__":
    plot_CV_all(save=1,show=1)
    # plot_dist_profile(save=1,show=1)
    # plot_src_lc_singleobs(figurepath=path_out,save=1,show=1)
    # plot_CR_all()
    # plot_CR_GCLW(save=1,show=1)
    # plot_CV_Temp(save=1,show=1)
    # plot_surface_density(save=1,show=1)