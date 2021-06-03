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
# import mpmath
font1 = {'family': 'Normal',
         'weight': 'normal',
         'size': 18, }
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

    return (T,E)

def plot_longT_V(data_file,bkg_file,epoch_file):
    epoch_info = np.loadtxt(epoch_file)
    t_start = epoch_info[:, 0]
    t_end = epoch_info[:, 1]
    obsID = epoch_info[:, 2]
    expT = epoch_info[:, 3]
    time = np.loadtxt(data_file)
    bkg_time=np.loadtxt(bkg_file)
    cts=[];bkg_cts=[]
    for i in range(len(obsID)):
        cts.append(len(np.where(time[:,2]==obsID[i])[0]))
        bkg_cts.append((len(np.where(time[:,2]==obsID[i])[0])))
    cts=np.array(cts);bkg_cts=np.array(bkg_cts)
    CR=(cts-16./25*bkg_cts)/expT
    print(obsID[np.where(CR>0.0006)])
    plt.plot(t_start,CR,marker='+')
    plt.show()

def phase_fold(data_file,bkg_file,epoch_file,p_test,bin,shift=0):
    srcevt=np.loadtxt(data_file)
    bkgevt=np.loadtxt(bkg_file)
    epoch = np.loadtxt(epoch_file)
    useid = epoch[:, 2][8:12]


    time = [];
    bkg_time = [];
    energy = [];
    bkg_energy = []
    srcevt[:, -1] = srcevt[:, -1].astype('int');
    bkgevt[:, -1] = bkgevt[:, -1].astype('int')

    for k in range(len(useid)):
        time = np.concatenate((time, srcevt[:, 0][np.where(srcevt[:, -1] == useid[k])]))
        energy = np.concatenate((energy, srcevt[:, 1][np.where(srcevt[:, -1] == useid[k])]))
        bkg_time = np.concatenate((bkg_time, bkgevt[:, 0][np.where(bkgevt[:, -1] == useid[k])]))
        bkg_energy = np.concatenate((bkg_energy, bkgevt[:, 1][np.where(bkgevt[:, -1] == useid[k])]))
    (time,energy)=filter_energy(time,energy,[1000,7000])
    (bkg_time,bkg_energy)=filter_energy(bkg_time,bkg_energy,[1000,7000])

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
    turns_b=trans(bkg_time,p_test,shift)
    loc=np.zeros(bin);loc_b=np.zeros(bin)
    loc_b= np.zeros(bin)
    for index in turns:
        loc[int(index*bin)] += 1
    for index in turns_b:
        loc_b[int(index * bin)] += 1

    print(loc)
    print(loc_b)
    # loc=loc-loc_b
    # bkg_y = len(time) * src_bkg
    # bkg_y_low=bkg_y-bkg_y**0.5
    # bkg_y_high = bkg_y+bkg_y**0.5

    # bkg_y /= bin
    # b_1sigma = poisson_conf_interval(bkg_y, interval='frequentist-confidence').T
    # bkg_y_low=b_1sigma[0]
    # bkg_y_high=b_1sigma[1]
    # bkg_x = [0, 2]

    fig=plt.figure(1,(10,7.5))
    ax1 = fig.add_subplot(111)

    # plt.fill_between(bkg_x, bkg_y_low, bkg_y_high,facecolor = 'blue', alpha = 0.5)
    # y3=np.concatenate((loc_b,loc_b))
    # y3*=2.1590383026593E-05/5.0169307007741E-05

    x = np.linspace(0, 1, bin + 1)[:-1]+0.5/bin
    x2=np.concatenate((x,x+1))
    y2=np.concatenate((loc,loc))

    y2_err=np.array(poisson_conf_interval(y2,interval='frequentist-confidence'))
    y2_err[0]=y2-y2_err[0]
    y2_err[1]=y2_err[1]-y2

    #label='F1'
    # plt.title("#{0} P={1:.2f},C={2}".format(label,p_test,str(len(time))), fontsize = 18)
    plt.xlabel('phase',font1)
    plt.ylabel('counts/bin',font1)
    plt.tick_params(labelsize = 18)
    plt.ylim(0,(np.max(y2)+np.max(y2)**0.5)*1.05)

    plt.step(x2, y2, where='mid', label='mid',color='red')
    #plt.errorbar(x2 - 0.5 / bin, y3, yerr = y3 ** 0.5, fmt = '.', capsize = 1, elinewidth = 0.5, ecolor = 'green')
    plt.errorbar(x2 , y2, yerr = y2_err, fmt = '.', capsize = 1, elinewidth = 1, ecolor = 'red')

    ax2 = ax1.twinx()
    x22=x2
    y22=y2/np.mean(y2)
    yhigh=(np.max(y2)+np.max(y2)**0.5)*1.05/np.mean(y2)
    ax2.set_ylabel('Normalized flux',font1)
    ax2.plot([0,2],[1.0,1.0],'--',color='green')
    ax2.set_ylim([0,yhigh])
    ax2.tick_params(labelsize=18)

    path_out='/Users/baotong/Desktop/CDFS/report1/'
    plt.show()
    #plt.close()
# phase_fold('513_bkg.txt','LW_epoch.txt',5334.75593,bin = 30, net_percent = 0.9, shift = 0.2, label = 12)

path_xmmCDFS='/Users/baotong/Desktop/CDFS/xmm_txt/'
path=path_xmmCDFS
period=994.78
dataname='XID89_pn_all_obs.txt'
bkgname='bkg_XID89_pn_all_obs.txt'
epochname='epoch_XID89_pn_all_obs.txt'
net_p=0.9

phase_fold(path + dataname, path+bkgname, path+epochname,period, bin = 10, shift = 0.6)
# plot_longT_V(path+dataname,path+bkgname,path+epochname)