#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import pylab as pl
import string
from datetime import datetime
from scipy.interpolate import lagrange
from scipy import interpolate
from scipy.fftpack import fft,ifft
import scipy.signal as ss
import random
import pandas as pd
import read_csv as data
ID_LW=data.ID_LW
P_LW=data.P_LW
label_LW=data.label_LW

ID_LW_old = []
P_LW_old = []
ID_LW_new = []
P_LW_new = []

def plot_Pdet_Hong():
    path_LW = '/Users/baotong/Desktop/period_LW/simulation/Pdet_Hong/'
    # for i in range(len(label_LW)):
    #     if label_LW[i] == 0:
    #         ID_LW_new.append(ID_LW[i])
    #         P_LW_new.append(P_LW[i])
    #     else:
    #         ID_LW_old.append(ID_LW[i])
    #         P_LW_old.append(P_LW[i])
    # ID_LW_old =np.array(ID_LW_old)
    # P_LW_old = np.array(P_LW_old)
    # Pdet_old=np.array([0.801,0.096,0.999,0.139,0.675,0.661,0.991,0.833,0.998,0.257])

    ID_LW_old=np.array([153002,20,194,191,196,16,114,206,152,42])
    P_LW_old=np.array([10342.3,5130.57309,7448.98,8546.2781,6335.85,4728.902,12002.7,4886.79,6597.55,5261.93])
    Pdet_old=np.array([0.998,0.999,0.991,0.833,0.675,0.801,0.257,0.096,0.661,0.139])

    print(ID_LW_old)
    sim_N=100
    threshold=0.99
    detect_rate=[]

    for i in range(len(ID_LW_old)):
        detect=0.0
        period_real=P_LW_old[i]
        period_get=[]
        for j in range(1,sim_N+1):
            temp_info = np.loadtxt(path_LW + '{0}/result_sim_{1}.txt'.format(str(ID_LW_old[i]),str(j)))
            if temp_info[2] > threshold and 0.99 * period_real < temp_info[4] < 1.01 * period_real:
                detect += 1
                period_get.append(temp_info[4])

        detect_rate.append(detect/sim_N)
    #plt.figure(1,(9,7))
    x=np.linspace(1,10,10)
    xtick=[1,2,3,4,5,6,7,8,9,10]
    plt.xticks(xtick)
    plt.tick_params(labelsize=15)
    plt.plot(x,detect_rate,marker='o')
    plt.plot(x,Pdet_old,marker='o')
    plt.xlabel('Source ID',fontsize=15)
    plt.ylabel('Detection rate',fontsize=15)
    plt.legend(['GL method','LS method'])
    # for k in range(len(x)):
    #     plt.plot([x[k],x[k]],[0,Pdet_old[k]],'--')
    plt.savefig('/Users/baotong/Desktop/aas/mod_MN_pCV/figure/sim_LW/Pdet_com.eps',bbox_inches='tight')
    plt.show()
    print(Pdet_old)
    print(detect_rate)
plot_Pdet_Hong()

def plot_P_P():
    path_out = '/Users/baotong/Desktop/aas/mod_MN_pCV/figure/LW/'
    P_LW_LS = np.array([10342.30,5131.14,7448.40,8535.67,6341.54,4728.53,12075.81,4890.23,6597.41,5262.05])
    P_LW_GL = np.array([10342.30,5130.57,7448.98,8546.28,6335.85,4728.90,12002.70,4886.79,6597.55, 5261.93])
    ratio=P_LW_LS/P_LW_GL
    #plt.figure(1, (7,5))
    #plt.figure(1,(9,7))
    plt.ylim(0.99,1.01)
    ytick=[0.990,0.995,1.000,1.005,1.010]
    plt.yticks(ytick)
    plt.tick_params(labelsize=15)
    plt.xlabel(r'$P_{GL}$',fontsize=15)
    plt.ylabel(r'$\frac{P_{LS}}{P_{GL}}$',rotation='horizontal',fontsize=20,horizontalalignment='right')
    #plt.ylabel(r'$\frac{P_{LS}}{P_{GL}}$',fontsize=15)
    plt.scatter(P_LW_GL,ratio)
    plt.plot([4000,12500],[1,1],'--',color='orange')
    plt.savefig(path_out+'P_comp.eps',bbox_inches='tight')
    plt.show()
#plot_P_P()

def Pdet_LW():
    path_LW = '/Users/baotong/Desktop/period_LW/simulation/Pdet_LW/'
    sim_N=100
    threshold=0.9
    detect_rate=[]

    for i in range(len(ID_LW)):
        detect=0.0
        period_real=P_LW[i]
        period_get=[]
        for j in range(1,sim_N+1):
            temp_info = np.loadtxt(path_LW + '{0}/result_sim_{1}.txt'.format(str(ID_LW[i]),str(j)))
            if temp_info[2] > threshold and 0.001 * period_real < temp_info[4] < 1.001 * period_real:
                detect += 1
                period_get.append(temp_info[4])

        detect_rate.append(detect/sim_N)

    print(detect_rate)

#Pdet_LW()


