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

font1 = {'family': 'Normal',
         'weight': 'normal',
         'size': 18, }
plt.rc('legend',fontsize=15 )

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
def plot_RK_CV():
    bins=np.logspace(np.log10(0.5), np.log10(100), 51)
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
    ax1.set_xscale('log');ax1.set_yscale('log')
    ax1.set_xlim(7e-2,60)
    #print(len(np.where(spin_IP > 0)[0]))
    #plt.legend(['NSC','LW','Polar','DN','IP','Spin of IP'])
    ax1.legend(['Polar', 'DN', 'IP', 'Spin of IP'])
    P_min = 7./6.
    P_gap = [7740.0 / 3600., 11448.0 / 3600.]
    ax1.plot([P_min, P_min], [0, 200], '--',color='grey')
    ax1.plot([P_gap[0], P_gap[0]], [0, 200], '-',lw=2.,color='grey')
    ax1.plot([P_gap[1], P_gap[1]], [0, 200], '-',lw=2.,color='grey')

    ax1.text(P_gap[0]+0.1, 230, 'gap',fontsize=12)
    ax1.text(P_min - 0.41, 220, 'minimum',fontsize=12)

    # plt.show()
    return (fig,ax1)

def read_excel(label):
    res = pd.read_excel(path + 'result_0.5_8_all.xlsx', label)
    ra = np.array(res['RA'])
    dec = np.array(res['DEC'])
    seq = np.array(res['seq'])
    period = np.array(res['P_out'])
    L = np.array(res['L'])
    Lmin = np.array(res['Lmin'])
    Lmax = np.array(res['Lmax'])
    type = np.array(res['type'])

    return (ra,dec,seq,period,L,Lmin,Lmax,type)

def plot_NP(save=1,show=1):
    (fig,ax1)=plot_RK_CV()
    (ra, dec, seq, period, L, Lmin, Lmax, type)=read_excel('47Tuc')
    period/=3600.
    print(len(period))
    period=period[np.where(type=='CV')]
    bins_p=np.logspace(np.log10(0.2), np.log10(20), 21)
    ax1.hist(period,bins=bins_p,histtype = 'step',lw=3,color='black')
    ax1.set_xlabel('Period (hours)',font1)
    ax1.set_ylabel('Number of sources',font1)
    plt.tick_params(labelsize=16)
    ax1.legend(['Polar', 'DN', 'IP', 'Spin of IP','47 Tuc'])
    if save:
        plt.savefig(path_out+'47Tuc_NP.eps',bbox_inches='tight', pad_inches=0.0)
    if show:
        plt.show()

    return None

def plot_dist_profile():
    (ra, dec, seq, period, L, Lmin, Lmax, type)=read_excel('47Tuc')
    [ra_center,dec_center,rhl,rc]=pos_all['47Tuc']
    c1 = SkyCoord(ra * u.deg, dec * u.deg, frame='fk5')
    c2 = SkyCoord(ra_center * u.deg, dec_center * u.deg, frame='fk5')
    dist = c1.separation(c2)
    dist = dist.arcmin

    period=period[np.where(type=='CV')]
    dist=dist[np.where(type=='CV')]
    L=L[np.where(type=='CV')]

    fig=plt.figure(1,figsize=(9,6))
    ax2=fig.add_subplot(211)
    ax2.scatter(period/3600.,L,marker='v',s=80,color='w',linewidths=2,edgecolors='red')
    ax2.set_yscale('log')
    ax2.set_ylabel(r'Luminosity ($\rm erg~s^{-1}$)',font1)
    ax2.tick_params(labelsize=16)
    P_gap = [7740.0 / 3600., 11448.0 / 3600.]
    x = P_gap
    y1 = [0, 5e32]
    # ax2.text(7900 / 3600., 5e36, 'period gap')
    ax2.fill_between(x, y1[1], facecolor='yellow', alpha=0.2)

    ax3=fig.add_subplot(212)
    ax3.scatter(period/3600.,dist,marker='v',s=80,color='w',linewidths=2,edgecolors='red')
    ax3.plot([0,15],[rhl/60,rhl/60],'--')
    ax3.plot([0,15],[rc/60,rc/60],'--')
    ax3.set_yscale('log')
    ax3.set_ylabel(r'R ($\rm arcmin$)',font1)
    ax3.set_xlabel('Period (hours)',font1)
    ax3.tick_params(labelsize=16)
    ax2.get_shared_x_axes().join(ax2, ax3)
    y2 = [0, 4]
    # ax2.text(7900 / 3600., 5e36, 'period gap')
    ax3.fill_between(x, y2[1], facecolor='yellow', alpha=0.2)
    plt.savefig(path_out+'47Tuc_profile.pdf',bbox_inches='tight', pad_inches=0.0)
    plt.show()
if __name__=="__main__":
    # plot_NP()
    plot_dist_profile()
