#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import string
from datetime import datetime
from scipy.interpolate import lagrange
from scipy import interpolate
from scipy import optimize
from scipy.optimize import curve_fit
import pandas as pd
import linecache
from astropy.timeseries import LombScargle
from astropy.stats import poisson_conf_interval
from tkinter import _flatten

font1 = {'family': 'Normal',
         'weight': 'normal',
         'size': 18, }
def get_line_context(file_path, line_number):
    return linecache.getline(file_path, line_number).strip()


def get_VI(cts_rate_m, cts_rate_low, cts_rate_high):
    cts_rate = cts_rate_m
    low = cts_rate_low
    high = cts_rate_high
    i = 0
    if len(cts_rate) != len(low) != len(high):
        print('error')
        return None
    else:
        while i < len(low):
            if low[i] == 0:
                low = np.delete(low, i)
                high = np.delete(high, i)
                cts_rate = np.delete(cts_rate, i)
            else:
                i += 1
    if len(cts_rate) == 0:
        return 0
    VI = np.max(cts_rate) / np.min(cts_rate)
    return VI
def pfold(time,P,flux):
    def trans(t, p_test, shift=0.5):
        ti = t
        v = 1.0 / p_test
        turns = v * ti
        turns += shift
        # 初始相位
        for i in range(len(turns)):
            turns[i] = turns[i] - int(turns[i])
        return turns
    turns=trans(time,P)
    plt.figure(1,(8,6))
    plt.xlabel('phase')
    plt.ylabel('photon flux')
    #print(useobs[np.where((turns<0.5)&(turns>0.2))])
    plt.errorbar(turns,flux[0],yerr=[flux[0]-flux[1],flux[2]-flux[0]],
                 fmt = 'o', capsize = 3, elinewidth = 1, color='red',ecolor = 'red')
    plt.errorbar(turns+1., flux[0], yerr = [flux[0] - flux[1], flux[2] - flux[0]],
                 fmt = 'o', capsize = 3, elinewidth = 1, color = 'red', ecolor = 'red')
    #plt.savefig('/Volumes/pulsar/WR/1671/pfold_2_4.eps')
    plt.show()
def plot_lc_ap(mode):
    #path_out='/Users/baotong/Desktop/aas/MN_pCV/figure/'+mode+'/'
    if mode=='CDFS':
        src_ID=['711']
        path='/Users/baotong/Desktop/CDFS/aprates/'
        filename='CDFS_epoch.txt'
        obs_ID=np.loadtxt(path+filename)[:,2]
        obs_time=(np.loadtxt(path+filename)[:,0]+np.loadtxt(path+filename)[:,1])/2
        time=obs_time/86400+2449352.5-2400000.5
        os.chdir(path)

    for k in range(len(src_ID)):
        ID = int(src_ID[k])
        ID_note = int(src_ID[k])
        use_obs = []
        use_obs_time = []
        cts_rate_low = []
        cts_rate_high = []
        cts_rate = []
        for i in range(len(obs_ID)):
            outname = str(ID) + '_' + str(int(obs_ID[i])) + '_out.par'
            if os.path.exists(outname):
                use_obs.append(obs_ID[i])
                use_obs_time.append(time[i])
                a = get_line_context(outname, 15)[18:-3]
                a_low = get_line_context(outname, 16)[25:-3]
                a_high = get_line_context(outname, 17)[25:-3]
                if a == 'INDEF':
                    a = 0
                if a_low == 'INDEF':
                    a_low = 0
                if a_high == 0:
                    continue
                cts_rate.append(float(a))
                cts_rate_low.append(float(a_low))
                cts_rate_high.append(float(a_high))

            else:
                print('continue')
                continue
        VI = get_VI(cts_rate, cts_rate_low, cts_rate_high)
        cts_rate=np.array(cts_rate)*1e7
        cts_rate_low=np.array(cts_rate_low)*1e7
        cts_rate_high=np.array(cts_rate_high)*1e7

        path_aprates='/Users/baotong/Desktop/CDFS/aprates/'
        out_info = np.column_stack((cts_rate, cts_rate_low, cts_rate_high, use_obs))
        np.savetxt(path_aprates + str(ID) + '.txt', out_info, fmt='%10.3f %10.3f %10.3f %10d')

plot_lc_ap('CDFS')

def plot_lc(source_id):
    period =76175.45748114
    bin_len=5000
    path='/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8/'
    evt_file=np.loadtxt(path+'{0}.txt'.format(source_id))
    epoch_file=np.loadtxt(path+'epoch_src_{0}.txt'.format(source_id))
    aprates_file=np.loadtxt('/Users/baotong/Desktop/CDFS/aprates/{0}.txt'.format(source_id))
    EP1_file = np.loadtxt('/Users/baotong/Desktop/CDFS/aprates/CDFS_epoch_ep1.txt')
    EP2_file = np.loadtxt('/Users/baotong/Desktop/CDFS/aprates/CDFS_epoch_ep2.txt')
    EP3_file = np.loadtxt('/Users/baotong/Desktop/CDFS/aprates/CDFS_epoch_ep3.txt')
    EP4_file = np.loadtxt('/Users/baotong/Desktop/CDFS/aprates/CDFS_epoch_ep4.txt')
    index=[len(EP1_file),len(EP1_file)+len(EP2_file),len(EP1_file)+len(EP2_file)+len(EP3_file)]
    VI=get_VI(aprates_file[:,0],aprates_file[:,1],aprates_file[:,2])

    x1=(EP1_file[:,0]+EP1_file[:,1])/2;y1=aprates_file[0:index[0]][:,0];y1_low=aprates_file[0:index[0]][:,1];y1_high=aprates_file[0:index[0]][:,2]
    x2=(EP2_file[:,0]+EP2_file[:,1])/2;y2=aprates_file[index[0]:index[1]][:,0];y2_low=aprates_file[index[0]:index[1]][:,1];y2_high=aprates_file[index[0]:index[1]][:,2]
    x3=(EP3_file[:,0]+EP3_file[:,1])/2;y3=aprates_file[index[1]:index[2]][:,0];y3_low=aprates_file[index[1]:index[2]][:,1];y3_high=aprates_file[index[1]:index[2]][:,2]
    x4=(EP4_file[:,0]+EP4_file[:,1])/2;y4=aprates_file[index[2]:][:,0];y4_low=aprates_file[index[2]:][:,1];y4_high=aprates_file[index[2]:][:,2]

    x1 = x1/ 86400 + 2449352.5 - 2400000.5;x2 = x2 / 86400 + 2449352.5 - 2400000.5;x3=x3 / 86400 + 2449352.5 - 2400000.5;x4 = x4 / 86400 + 2449352.5 - 2400000.5
    # gap1=np.mod(EP2_file[:,0][0]-EP1_file[:,1][-1],period);gap2=np.mod(EP3_file[:,0][0]-EP2_file[:,1][-1],period);gap3=np.mod(EP4_file[:,0][0]-EP3_file[:,1][-1],period)
    # x2-=(EP2_file[:,0][0]-EP1_file[:,1][-1]-gap1)/86400
    # x3-=(EP3_file[:,0][0]-EP2_file[:,1][-1]-gap2)/86400
    # x4-=(EP4_file[:,0][0]-EP3_file[:,1][-1]-gap3)/86400
    color_use=['red','blue','green','purple']
    pfold((epoch_file[:,0]+epoch_file[:,1])/2,period,[aprates_file[:,0],aprates_file[:,1],aprates_file[:,2]])
    fig= plt.figure(1)
    plt.semilogy()
    plt.title('{0}, VI={1}'.format(source_id,round(VI, 2)), fontsize=20)
    for i in range(4):
        x=[x1,x2,x3,x4][i];y=[y1,y2,y3,y4][i];y_low=[y1_low,y2_low,y3_low,y4_low][i];y_high=[y1_high,y2_high,y3_high,y4_high][i]
        for j in range(len(x)):
            if y_low[j] == 0:
                plt.annotate("", xy = (x[j], y_high[j] * 0.5),
                             xytext = (x[j], y_high[j]), color = color_use[i],
                             weight = "bold",
                             arrowprops = dict(arrowstyle = "->", connectionstyle = "arc3", color =color_use[i]))
            else:
                plt.errorbar(x[j], y[j],
                        yerr = [[y[j] - y_low[j]], [y_high[j] - y[j]]],
                        fmt = 'o', capsize = 3, elinewidth = 1, color=color_use[i],ecolor = color_use[i])
    plt.xlabel('time')
    plt.ylabel('counts rate')
    plt.show()
    # plt.errorbar(x,y,xerr=xerr,yerr=yerr)
    # plt.show()
plot_lc('711')
