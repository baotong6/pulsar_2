#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import string
#import correct as correct
from datetime import datetime
from scipy.interpolate import lagrange
from scipy import interpolate
from scipy.optimize import curve_fit
import pandas as pd
from astropy.stats import poisson_conf_interval
import scipy
import hawkeye as hawk
font1 = {'family': 'Normal',
         'weight': 'normal',
         'size': 18, }
plt.rc('legend',fontsize=15 )
path='/Users/baotong/Desktop/period_Tuc/'
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
    #print(len(np.where(spin_IP > 0)[0]))
    #plt.legend(['NSC','LW','Polar','DN','IP','Spin of IP'])
    ax1.legend(['Polar', 'DN', 'IP', 'Spin of IP'])
    # plt.show()
    return ax1

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
    path_out='/Users/baotong/Desktop/aas/pXS_Tuc/figure/'
    ax1=plot_RK_CV()
    (ra, dec, seq, period, L, Lmin, Lmax, type)=read_excel('47Tuc')
    period/=3600.
    print(len(period))
    period=period[np.where(type=='CV')]
    print(period*3600)
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
if __name__=="__main__":
    plot_NP()