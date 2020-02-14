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
import read_csv as data
# from scipy.interpolate import Spline
import spline

path_fits='/Users/baotong/Desktop/period_LW/'
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

def plot_P_N():
    orb_DN=orb[np.where(type1=='DN')]
    orb_IP=orb[np.union1d(np.where(type2=='IP'),np.where(type3=='IP'))]
    orb_LW=data.P_LW/3600.
    bins=np.logspace(np.log10(0.5), np.log10(100), 51)
    bins_2=np.logspace(np.log10(0.5), np.log10(15), 21)
    LW_plt=plt.hist(orb_LW,bins=bins_2,histtype = 'step')
    plt.close()

    plt.figure(1,(10,7))
    P_min = 4944.0 / 3600.
    P_gap = [7740.0 / 3600., 11448.0 / 3600.]

    plt.step(LW_plt[1],np.append(LW_plt[0],0),linestyle='-',lw=3,where='post')
    plt.loglog()

    plt.hist(orb,bins=bins,histtype = 'step',color='black')
    plt.hist(orb_DN,bins=bins,histtype = 'step',color='red')
    plt.hist(orb_IP,bins=bins,histtype = 'step',color='green')

    plt.legend(['LW','CV','DN','IP'])

    plt.plot([P_min, P_min], [0, 100], '--')
    plt.plot([P_gap[0], P_gap[0]], [0, 200], '-',lw=2.,color='yellow')
    plt.plot([P_gap[1], P_gap[1]], [0, 200], '-',lw=2.,color='yellow')
    plt.text(P_gap[0]-0.15, 220, 'period gap')
    plt.text(P_min - 0.41, 100, 'period minum')

    plt.xlabel('period (hours)')
    plt.ylabel('number of sources')
    plt.savefig(path_fits+'N_P.eps')
    plt.show()

def plot_P_M():

    orb_DN_M1=orb[np.intersect1d(np.where(type1=='DN'),np.where(M1!=0))]
    # M1TOM2_DN=M1_M2[np.intersect1d(np.where(type1=='DN'),np.where(M1_M2!=0))]
    M1_DN=M1[np.intersect1d(np.where(type1=='DN'),np.where(M1!=0))]

    orb_IP_spin_M1=orb[np.intersect1d(np.union1d(np.where(type2=='IP'),np.where(type3=='IP')),np.where(M1!=0),np.where(spin!=0))]
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

plot_GL()
