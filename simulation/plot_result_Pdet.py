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
            if temp_info[2] > threshold and 0.01 * period_real < temp_info[4] < 1.01 * period_real:
                detect += 1
                period_get.append(temp_info[4])

        detect_rate.append(detect/sim_N)

    x=np.linspace(1,10,10)
    plt.plot(x,detect_rate,marker='o')
    plt.plot(x,Pdet_old,marker='o')
    plt.xlabel('Source ID')
    plt.ylabel('Detection Rate')
    plt.legend(['GL method','LS method'])
    plt.savefig('/Users/baotong/Desktop/period_LW/simulation/fig_Pdet/Pdet_com.eps')
    plt.show()
    print(Pdet_old)
    print(detect_rate)
plot_Pdet_Hong()

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


