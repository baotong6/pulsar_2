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
# import read_data as data
import scipy.stats as stats

def trans(t,p_test,shift=0):
    ti=t
    #print ti
    #p_test = 1.0/5.500550055005501e-06
    #p_test=4945.45
    v = 1.0 /p_test
    p=1.0/v
    # pdot=-vdot/(v*v)
    # vddot = 2.0 * pdot * pdot / (p * p * p)
    turns = v * ti #+ vdot * ti * ti / 2.0 + vddot * ti * ti * ti / 6.0
    turns += shift
    #初始相位
    for i in range(len(turns)):
        turns[i]=turns[i] - int(turns[i])
    return turns

def get_T_in_mbins(epoch_file,w,m,fi=0.):
    #m==bin
    T=2*np.pi/w
    T_in_perbin = np.zeros(m)
    # 每个bin的总积分时间
    tbin = T/m
    # 每个bin的时间长度
    epoch_info = np.loadtxt(epoch_file)
    t_start = np.array([epoch_info[0]])
    t_end = np.array([epoch_info[1]])
    N_bin_t_start=t_start/tbin+m*fi/(2*np.pi)
    N_bin_t_end=t_end/tbin+m*fi/(2*np.pi)
    intN_bin_t_start=np.floor(N_bin_t_start)+1
    intN_bin_t_end=np.floor(N_bin_t_end)
    intN_bin_t_start=intN_bin_t_start.astype(int)
    intN_bin_t_end=intN_bin_t_end.astype(int)
    for i in range(len(N_bin_t_start)):
        if intN_bin_t_end[i]>=intN_bin_t_start[i]:
            T_in_perbin+=int((intN_bin_t_end[i]-intN_bin_t_start[i])/m)*tbin
            #print(intN_bin_t_start[i]-1)
            T_in_perbin[np.mod(intN_bin_t_start[i],m)-1]+=(intN_bin_t_start[i]-N_bin_t_start[i])*tbin
            T_in_perbin[np.mod(intN_bin_t_end[i],m)]+=(N_bin_t_end[i]-intN_bin_t_end[i])*tbin
            rest=np.mod(intN_bin_t_end[i]-intN_bin_t_start[i],m)
            for k in range(rest):
                T_in_perbin[int(np.mod((intN_bin_t_start[i] + k), m))] += tbin
            #print(rest)
        else:
            T_in_perbin[np.mod(intN_bin_t_start[i],m)-1]+=(N_bin_t_end[i]-N_bin_t_start[i])*tbin
    return T_in_perbin

def compute_bin(Tlist, m, w, fi=0):
    n = np.zeros(m, 'int')
    j = np.floor(m * np.mod(w * Tlist + fi, 2 * np.pi) / (2 * np.pi))
    j.astype(int)
    for u in range(0, m):
        n[u] = np.size(np.extract(j == u, j))
    return n

# def get_S_inEP(time,freq,bin):
#     turns=trans(time,1./freq)
#     for i in len(turns):

path='/Volumes/pulsar/xmm_CV/result/'

name=['HT_CAS','OY_CAR','QZ_VIR','RU_PEG','SS_AUR',
      'V893_SCO','VW_HYI','YZ_CNC','V405_PEG','AO_PSC',
      'DW_CNC','FO_AQR','HT_CAM','J1509_6649','J1649_3307',
      'J1719_4100','J1817_2508','J1830_1232','XY_ARI']
period_all=[6363.1008,5453.65,414.,32365.44,15793.91,
        6563.0304,6417.0144,7499.52,15348.7008,12927.6864,
        5166.1152,17457.984,5159.1168,21202.56,13020.48,
        14420.16,5512.32,19344.96,21833.0208]
i=16
dataname=name[i]+'_pn'+'.txt'
epoch_file=path+'epoch_'+dataname
period=period_all[i]
time=np.loadtxt(path+dataname)[:,0]
p_test=period
P_all=np.linspace(0.99*p_test,1.01*p_test,5000)
bin=100
S=[]
for i in range(len(P_all)):
    w=2*np.pi/P_all[i]
    cts_in_bin=compute_bin(time,bin,w)
    T_in_perbin=get_T_in_mbins(epoch_file,w,bin)
    cts_in_bin=np.array(cts_in_bin)
    T_in_perbin=np.array(T_in_perbin)
    R=len(time)/(sum(T_in_perbin))
    Rj=cts_in_bin/T_in_perbin
    #print(Rj)
    S.append(sum((Rj-R)**2/(R/T_in_perbin)))
plt.title('EP for {0}'.format(dataname[:-4]))
plt.xlabel('time')
plt.ylabel('S_value')
plt.plot([P_all[0],P_all[-1]],[60.4,60.4],'--')
plt.legend(['99%'])
plt.step(P_all,S,'red')
#path_fig='/Users/baotong/Desktop/li_pr/result_final/fig_ND/EP/'
#plt.savefig(path_fig+'EP_for_{0}.eps'.format(dataname[:-4]))
plt.show()
