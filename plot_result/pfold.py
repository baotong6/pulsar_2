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
import read_csv as data

def filter_energy(time,energy,band):
    T=time
    E=energy
    i=0
    if len(time)!=len(energy):
        print('error')
        return None
    else:
        while i <len(E):
            if E[i]<band[0] or E[i]>band[1]:
                E=np.delete(E,i)
                T=np.delete(T,i)
            else:
                i+=1

    return T

def get_T_in_mbins(epoch_file,w,m,fi):
    T=2*np.pi/w
    T_in_perbin = np.zeros(m)
    # 每个bin的总积分时间
    tbin = T/m
    # 每个bin的时间长度
    epoch_info = np.loadtxt(epoch_file)
    t_start = epoch_info[:, 0]
    t_end = epoch_info[:, 1]
    ID = epoch_info[:, 2]
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

def phase_fold(data_file,epoch_file,p_test,bin,net_percent,shift,label):
    time=np.loadtxt(data_file)[:,0]
    time_bkg=np.loadtxt('513_bkg.txt')[:,0]
    energy=np.loadtxt(data_file)[:,1]
    time=filter_energy(time,energy,[1000,8000])
    T_in_perbin = get_T_in_mbins(epoch_file, 2 * np.pi / p_test, bin, 0.0)
    def trans(t,p_test,shift):
        ti =t
        v = 1.0 /p_test
        turns = v * ti
        turns += shift
        # 初始相位
        for i in range(len(turns)):
            turns[i] = turns[i] - int(turns[i])
        return turns

    turns=trans(time,p_test,shift)
    turns_b=trans(time_bkg,p_test,shift)
    loc=np.zeros(bin)
    loc_b= np.zeros(bin)
    for index in turns:
        loc[int(index*bin)] += 1
    for index in turns_b:
        loc_b[int(index * bin)] += 1

    x = np.array([(i / bin + 0.5 / bin) for i in range(bin)])
    AM=1-min(loc)/max(loc)
    A0=AM/(2-AM)
    print('A0={0}'.format(A0))

    src_bkg = 1 - net_percent
    bkg_y = len(time) * src_bkg
    bkg_y_low=bkg_y-bkg_y**0.5
    bkg_y_high = bkg_y+bkg_y**0.5

    bkg_y /= bin
    bkg_y_low/=bin
    bkg_y_high/=bin

    fig=plt.figure()
    ax1 = fig.add_subplot(111)

    #plt.plot([0, 2], [bkg_y, bkg_y], '--')
    bkg_x = [0, 2]
    plt.fill_between(bkg_x, bkg_y_low, bkg_y_high,facecolor = 'yellow', alpha = 0.2)
    # y3=np.concatenate((loc_b,loc_b))
    # y3*=2.1590383026593E-05/5.0169307007741E-05


    x2=np.concatenate((x,x+1))
    y2=np.concatenate((loc,loc))
    correct_gap = T_in_perbin / (sum(T_in_perbin) / len(T_in_perbin))
    y2 /= np.concatenate((correct_gap, correct_gap))

    #plt.figure(1,(8,8))
    plt.title("#{0} P={1},cts={2}".format(label,str(p_test),str(len(time))), fontsize = 18)
    plt.xlabel('phase')
    plt.ylabel('counts/bin')
    plt.ylim(0,(np.max(y2)+np.max(y2)**0.5)*1.05)
    plt.step(x2,y2,color='red')
    #plt.step(x2,y3,'--',color='green')
    #plt.errorbar(x2 - 0.5 / bin, y3, yerr = y3 ** 0.5, fmt = '.', capsize = 1, elinewidth = 0.5, ecolor = 'green')
    plt.errorbar(x2 - 0.5 / bin, y2, yerr = y2 ** 0.5, fmt = '.', capsize = 1, elinewidth = 1, ecolor = 'red')

    ax2 = ax1.twinx()
    x22=x2
    y22=y2/np.mean(y2)
    yhigh=(np.max(y2)+np.max(y2)**0.5)*1.05/np.mean(y2)
    ax2.set_ylabel('Normalized flux')
    ax2.plot([0,2],[1.0,1.0],'--',color='green')
    ax2.set_ylim([0,yhigh])

    path_out='/Users/baotong/Desktop/aas/V63/figure/LW/'

    plt.savefig(path_out+'pfold_lc_{0}_second.eps'.format(dataname[0:-4]))
# phase_fold('513_bkg.txt','LW_epoch.txt',5334.75593,bin = 30, net_percent = 0.9, shift = 0.2, label = 12)

path_NSC='/Users/baotong/Desktop/period/txt_all_obs_IG/'
path_LW='/Users/baotong/Desktop/period_LW/txt_all_obs/'
path_ND='/Users/baotong/Desktop/period_gc/txt_all_obs/'
#epoch_file=path+'txt_90/'+'SgrA_I_epoch.txt'
#epoch_file=path+'txt_2_8k/'+'SgrA_I_epoch.txt'
path=path_LW
period=5130.57309*2
dataname='20.txt'
net_p=0.892
epoch_file = path + 'epoch_src_' + dataname
phase_fold(path + dataname, epoch_file, period, bin = 20, net_percent = net_p, shift = 0.2, label = '6')
# for i in range(len(data.ID_LW)):
# index=23
# for i in range(index,index+1):
#     period=data.P_LW[i]
#     net_p=data.net_percent_LW[i]
#     dataname=str(data.ID_LW[i])+'.txt'
#     if dataname[-7:-4]=='001' or dataname[-7:-4]=='002':
#         dataname=dataname[0:-7]+dataname[-4:]
#     epoch_file = path + 'epoch_src_' + dataname
#     phase_fold(path + dataname, epoch_file, period, bin = 20, net_percent = net_p, shift = 0.5, label = i + 1)
#
plt.show()