#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import string
from scipy.interpolate import lagrange
from scipy import interpolate
from scipy.optimize import curve_fit
import pandas as pd
# import read_csv as data
# from scipy.interpolate import Spline

path_fits='/Users/baotong/Desktop/period_LW/'
path_fig='/Users/baotong/Desktop/aas/NSC_pCV/figure/'
RK=fits.open(path_fits+'RK14.fit')
orb=RK[1].data['Orb_Per']
type1=RK[1].data['Type1']
type2=RK[1].data['Type2']
type3=RK[1].data['Type3']
M1=RK[1].data['M1']
M2=RK[1].data['M2']
M1_M2=RK[1].data['M1_M2']
orb=orb*24
spin=RK[1].data['_3___Per']

Lum_LW=np.array([31.26,32.36,31.82,32.35,31.75,31.97
                        ,31.83,31.25,31.54,31.33,31.98
                        ,31.88,31.93,32.11,31.93
                        ,32.05,31.62,32.92,31.73,32.78,31.97,31.94,30.88])
orb_LW=data.P_LW/3600.
orb_NSC=data.P_NSC_IG/3600.

def plot_P_N():
    orb_DN=orb[np.where(type1=='DN')]
    orb_Polar=orb[np.where((type2=='AM')|(type3=='AM'))]
    orb_IP=orb[np.union1d(np.where((type2=='IP')|(type2=='DQ')),np.where((type3=='IP')|(type3=='DQ')))]
    #print(len(orb_IP))
    spin_IP=spin[np.union1d(np.where((type2=='IP')|(type2=='DQ')),np.where((type3=='IP')|(type3=='DQ')))]
    spin_IP/=3600.
    #print(len(spin_IP))
    #print(orb_LW)

    bins=np.logspace(np.log10(0.5), np.log10(100), 51)
    bins_2=np.logspace(np.log10(0.1), np.log10(20), 31)
    bins_spin=np.logspace(np.log10(3/36.), np.log10(2), 21)
    orb_LW32=[]
    orb_LW33=[]
    for i in range(len(Lum_LW)):
        if Lum_LW[i]>32:
            orb_LW33.append(orb_LW[i])
        else:
            orb_LW32.append(orb_LW[i])
    LW_plt=plt.hist(orb_LW,bins=bins_2,histtype = 'step')
    LW_plt32=plt.hist(orb_LW32,bins=bins_2,histtype = 'step')
    LW_plt33= plt.hist(orb_LW32, bins = bins_2, histtype = 'step')
    NSC_plt=plt.hist(orb_NSC, bins = bins_2, histtype = 'step')
    plt.close()

    plt.figure(1,(10,6))
    plt.xlim(0.07,100)
    P_min = 4920.0 / 3600.
    P_gap = [7740.0 / 3600., 11448.0 / 3600.]

    #plt.step(NSC_plt[1], np.append(NSC_plt[0], 0), linestyle='-', lw=4, where='post', color='black')

    #plt.step(LW_plt[1],np.append(LW_plt[0],0),linestyle='--',lw=3,where='post',color='brown')
    plt.loglog()
    spin_IP=spin_IP[np.where(spin_IP > 0)]
    # plt.hist(orb,bins=bins,histtype = 'step',color='black')
    plt.hist(orb_Polar,bins=bins,histtype='step',lw=2,color='blue')
    plt.hist(orb_DN,bins=bins,histtype = 'step',lw=2,color='red')
    plt.hist(orb_IP,bins=bins,histtype = 'step',lw=1,color='green')
    plt.hist(spin_IP, bins = bins_spin, histtype = 'step',lw=1, color = 'purple')

    print(len(spin_IP))
    print(len(orb_DN))
    #print(len(np.where(spin_IP > 0)[0]))

    #plt.legend(['NSC','LW','Polar','DN','IP','Spin of IP'])
    plt.legend(['Polar', 'DN', 'IP', 'Spin of IP'])
    #print(len(orb_Polar),len(orb_IP),len(orb_DN))
    plt.tick_params(labelsize = 15)
    plt.plot([P_min, P_min], [0, 200], '--',color='grey')
    plt.plot([P_gap[0], P_gap[0]], [0, 200], '-',lw=2.,color='grey')
    plt.plot([P_gap[1], P_gap[1]], [0, 200], '-',lw=2.,color='grey')
    plt.text(P_gap[0]+0.1, 220, 'gap',fontsize=12)
    plt.text(P_min - 0.41, 220, 'minimum',fontsize=12)

    plt.xlabel('Period (hours)',fontsize=20)
    plt.ylabel('Number of sources',fontsize=20)

    plt.savefig(path_fig+'N_P.eps')
    N=len(orb_IP)+len(orb_Polar)+len(orb_DN)
    #print(N)
    plt.show()
# plot_P_N()

def plot_P_RK():
    bins=np.logspace(np.log10(0.5), np.log10(100), 51)
    bins_2=np.logspace(np.log10(0.2), np.log10(14), 41)
    orb_DN = orb[np.where(type1 == 'DN')]
    orb_Polar = orb[np.where(type2 == 'AM')]
    orb_IP = orb[np.union1d(np.where((type2 == 'IP') | (type2 == 'DQ')), np.where((type2 == 'IP') | (type2 == 'DQ')))]
    orb_CV=np.concatenate((orb_DN,orb_Polar,orb_IP))

    orb_AMCV=orb[np.where(type1 == 'AC')]
    P_min = 7./6.
    P_gap = [7740.0 / 3600., 11448.0 / 3600.]

    print(len(orb_DN))
    print(len(np.where((orb_DN>P_gap[0])&(orb_DN<P_gap[1]))[0]))
    print(len(np.where(orb_Polar<P_gap[0])[0]))
    print(len(np.where(orb_IP<P_gap[0])[0]))
    print(len(orb_Polar),len(orb_IP))
    plt.figure(1,(9,6))
    plt.loglog()
    plt.xlabel('Period (hours)',fontsize=20)
    plt.ylabel('Number of sources',fontsize=20)

    # plt.hist(orb,bins=bins,histtype = 'step',color='black')
    # plt.hist(orb_Polar, bins=bins, histtype='step', lw=2, color='blue')
    # plt.hist(orb_DN, bins=bins, histtype='step', lw=2, color='red')
    # plt.hist(orb_IP, bins=bins, histtype='step', lw=2, color='green')
    plt.hist(orb_CV,bins=bins,histtype='step',lw=2,color='red')
    plt.hist(orb_AMCV,bins=bins_2,histtype='step',lw=2,color='black')
    plt.legend(['CV','AM CVn'])
    plt.plot([P_min, P_min], [0, 200], '--',color='grey')
    plt.plot([P_gap[0], P_gap[0]], [0, 200], '-',lw=2.,color='grey')
    plt.plot([P_gap[1], P_gap[1]], [0, 200], '-',lw=2.,color='grey')

    plt.text(P_gap[0]+0.1, 230, 'gap',fontsize=12)
    plt.text(P_min - 0.41, 220, 'minimum',fontsize=12)
    plt.savefig('/Users/baotong/Desktop/aas/作业格式/figure/'+'NP.eps')

    plt.show()

#plot_P_RK()

def plot_P_M():

    orb_DN_M1=orb[np.intersect1d(np.where(type1=='DN'),np.where(M1!=0))]
    # M1TOM2_DN=M1_M2[np.intersect1d(np.where(type1=='DN'),np.where(M1_M2!=0))]
    M1_DN=M1[np.intersect1d(np.where(type1=='DN'),np.where(M1!=0))]

    orb_IP_spin_M1=orb[np.intersect1d(np.union1d(np.where((type2=='IP')|(type2=='DQ')),np.where((type2=='IP')|(type2=='DQ'))),np.where(M1!=0),np.where(spin!=0))]
    spin_IP_M1=spin[np.intersect1d(np.union1d(np.where(type2=='IP'),np.where(type3=='IP')),np.where(M1!=0),np.where(spin!=0))]
    M1_spin_IP=M1[np.intersect1d(np.union1d(np.where(type2=='IP'),np.where(type3=='IP')),np.where(M1!=0),np.where(spin!=0))]
    plt.scatter(spin_IP_M1/3600.,orb_IP_spin_M1)
    # plt.scatter(M1_DN, orb_DN_M1 ** (-2))
    plt.show()

def plot_GL():
    a=np.random.random(10)*200
    x=np.linspace(0,0.9,10)
    plt.bar(x,a,width=0.1,align='edge',color='white',edgecolor='black')

    x+=0.05

    xnew = np.linspace(x.min(), x.max(), 300)  # 300 represents number of points to make between T.min and T.max
    func = interpolate.interp1d(x, a, kind = 'cubic')
    ynew = func(xnew)
    plt.plot(xnew, ynew)
    plt.xlim(0,1)
    plt.savefig(path_fits+'GL.eps')

    plt.show()

def get_dist_DN():
    orb_DN = orb[np.where(type1 == 'DN')]
    bins_DN = np.array([4944., 7740., 30000.])/3600.
    a =plt.hist(orb_DN,bins_DN,histtype = 'step')[0]

    plt.show()
    print(a)
#plot_P_N()
# get_dist_DN()
plot_P_RK()