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
from scipy.optimize import curve_fit
import pandas as pd
import linecache
from astropy.timeseries import LombScargle

font1 = {'family': 'Normal',
         'weight': 'normal',
         'size': 18, }

obs_ID_I = np.loadtxt('SgrA_I_epoch.txt')[:, 2]
obs_time_I = (np.loadtxt('SgrA_I_epoch.txt')[:, 0] + np.loadtxt('SgrA_I_epoch.txt')[:, 1]) / 2
time_I = obs_time_I / 86400 + 2449352.5 - 2400000.5

obs_ID_G = np.loadtxt('SgrA_G_epoch.txt')[:, 2]
obs_time_G = (np.loadtxt('SgrA_G_epoch.txt')[:, 0] + np.loadtxt('SgrA_G_epoch.txt')[:, 1]) / 2
time_G = obs_time_G / 86400 + 2449352.5 - 2400000.5

path_out='/Users/baotong/Desktop/aas/V63/figure/NSC_IG/'
def get_line_context(file_path, line_number):
    return linecache.getline(file_path, line_number).strip()

def get_plot_ap(ID):
    src_ID=ID
    label = int(src_ID)
    path_NSC_I='/Volumes/pulsar/SgrA/merge_data/spectra/spectra_p/aprates/'
    path_NSC_G='/Volumes/pulsar/SgrAGRT/merge_data/spectra/aprates/'
    path_aprates='/Users/baotong/Desktop/period/aprates/'
    def get_cts_rate(src_ID,obs_ID,path,time):
        os.chdir(path)
        ID = int(src_ID)
        label = int(src_ID)
        if str(int(ID))[-3:] == '001' or str(int(ID))[-3:] == '002':
            ID = str(ID)[:-3]
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

                cts_rate.append(float(a)*1e7)
                cts_rate_low.append(float(a_low)*1e7)
                cts_rate_high.append(float(a_high)*1e7)

            else:
                continue
        return (cts_rate,cts_rate_low,cts_rate_high,use_obs_time,use_obs)
    (CR_I,CR_I_L,CR_I_H,usetime_I,useobs_I)=get_cts_rate(src_ID,obs_ID_I,path_NSC_I,time_I)
    (CR_G,CR_G_L,CR_G_H,usetime_G,useobs_G)=get_cts_rate(src_ID,obs_ID_G,path_NSC_G,time_G)

    usetime = np.concatenate((usetime_I, usetime_G))
    CR_IG = np.concatenate((CR_I, CR_G))
    CR_low_IG=np.concatenate((CR_I_L,CR_G_L))
    CR_high_IG = np.concatenate((CR_I_H, CR_G_H))
    useobs=np.concatenate((useobs_I,useobs_G))
    out_info=np.column_stack((CR_IG,CR_low_IG,CR_high_IG,useobs))
    np.savetxt(path_aprates+str(ID)+'.txt',out_info,fmt='%10.3f %10.3f %10.3f %10d')

    def pfold(time,P,flux,):
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
    #pfold(usetime,7.5188,[CR_IG,CR_low_IG,CR_high_IG])

    def get_LS(time,flux):
        x=time
        y=flux
        freq=1./np.linspace(0.1, 100, 20000)
        print(time)

        LS = LombScargle(x, y, dy = 1, normalization = 'standard', fit_mean = True,
                         center_data = True).power(freq, method = 'cython')
        FP_99 = LombScargle(x, y, dy = 1, normalization = 'standard',
                            fit_mean = True, center_data = True).false_alarm_level(
            0.01,
            minimum_frequency = freq[0], maximum_frequency = freq[-1])
        FP_90 = LombScargle(x, y, dy = 1, normalization = 'standard',
                            fit_mean = True, center_data = True).false_alarm_level(
            0.1,
            minimum_frequency = freq[0], maximum_frequency = freq[-1])
        FP_68 = LombScargle(x, y, dy = 1, normalization = 'standard',
                            fit_mean = True, center_data = True).false_alarm_level(
            0.32,
            minimum_frequency = freq[0], maximum_frequency = freq[-1])

        plt.semilogx()
        plt.step(freq, LS)
        print(max(LS))
        print()
        plt.plot([freq[0], freq[-1]], [FP_99, FP_99], '--')
        plt.plot([freq[0], freq[-1]], [FP_90, FP_90], '--')
        plt.plot([freq[0], freq[-1]], [FP_68, FP_68], '--')
        #plt.savefig('/Volumes/pulsar/WR/1671/LS.eps')
        plt.show()
    #get_LS(usetime,CR_IG)

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
        if len(cts_rate)==0:
            return 0
        VI = np.max(cts_rate) / np.min(cts_rate)
        return VI
    VI=get_VI(np.concatenate((CR_I,CR_G)),np.concatenate((CR_I_L,CR_G_L)),np.concatenate((CR_I_H,CR_G_H)))
    #VI_I = get_VI(CR_I,CR_I_L,CR_I_H)
    #VI_G = get_VI(CR_G,CR_G_L,CR_G_H)

    plt.figure(1,(10,7.5))
    plt.semilogy()
    plt.xlim(49850.,55050.)
    plt.ylim(ymin=1e-1,ymax=np.max(CR_I_H))
    plt.ylim(ymin=1e-1, ymax=5e2)
    #plt.ylim(ymax=max(np.max(CR_I_H),np.max(CR_G_H)))
    for i in range(len(usetime_I)):
        #print(usetime_I[i])
        if CR_I_L[i] == 0:
            #ax.plot([use_obs_time[i] - 15, use_obs_time[i] + 15], [cts_rate_high[i], cts_rate_high[i]], color = 'red')
            plt.annotate("", xy = (usetime_I[i], CR_I_H[i] * 0.5),
                         xytext = (usetime_I[i], CR_I_H[i]), color = "red",
                         weight = "bold",
                         arrowprops = dict(arrowstyle = "->", connectionstyle = "arc3", color = "red"))

        else:
            plt.errorbar(usetime_I[i], CR_I[i],
                    yerr = [[CR_I[i] - CR_I_L[i]], [CR_I_H[i] - CR_I[i]]],
                    fmt = 'o', capsize = 3, elinewidth = 1, color='red',ecolor = 'red')

    for i in range(len(usetime_G)):
        if CR_G_L[i] == 0:
            #ax.plot([use_obs_time[i] - 15, use_obs_time[i] + 15], [cts_rate_high[i], cts_rate_high[i]], color = 'red')
            plt.annotate("", xy = (usetime_G[i], CR_G_H[i] * 0.5),
                         xytext = (usetime_G[i], CR_G_H[i]), color = "green",
                         weight = "bold",
                         arrowprops = dict(arrowstyle = "->", connectionstyle = "arc3", color = "green"))

        else:
            plt.errorbar(usetime_G[i], CR_G[i],
                    yerr = [[CR_G[i] - CR_G_L[i]], [CR_G_H[i] - CR_G[i]]],
                    fmt = 'o', capsize = 1, elinewidth = 1, color='green',ecolor = 'green')
    plt.title('{0} #{1}, VI={2}'.format('NSC', ID, round(VI, 2)), fontsize = 20)
    plt.xlabel('MJD', font1)
    plt.tick_params(labelsize=18)
    # plt.ylabel('photon flux ( ' + r'$\rm ph \,s^{-1} \, cm^{-2}$' + ')', fontsize = 15)
    plt.ylabel('photon flux ( ' + r'$\rm 10^{-7} ph \,s^{-1} \, cm^{-2}$' + ' )', font1)
    #plt.ylim(1,100)
    #plt.savefig(path_out + str(label) + '_lc.eps')
    plt.show()
    #plt.close()

ID_NSC_IG=['214001','1748','2560','116','1083','1080',
           '1502','2380','1873','3242','2268','3370',
           '2478','1182','2961','2532','1266','1206',
           '1180','3067','1854','442001','1624','2508',
           '3357','2841','1219','2672','2422','1487',
           '1853','3120','1133','2730','1538001','1084',
           '147','1514001','2525','1529','3483','2157',
           '6','2344','1674','442002','1525','2238',
           '1538002','1677','3564','2199','2338','3107',
           '1634','214002','3401','1628','790','1514002',
           '1941','1769','3596','2574','973','2187','307']

ID_noG=['116','3596','214','3564']
ID_faint=['1180','1182','1628','1538','2961','1487','1514']
ID_new=['2338','3107','2238','3401','1857','307','2380']
ID_new_2=['3370','147','1769','2574','1529','1941','1525',
             '6','3483','790','1674','442','3242','350','1748',
             '2268','2478','1854','1873','1080','1083']
ID_WR=['1671','2148']

#ID_trans=['1620','1612','1626','1906','2176','2178','2276','3236']
ID_trans=['1857']
for i in ID_trans:
    get_plot_ap(i)
# get_plot_ap('1487')