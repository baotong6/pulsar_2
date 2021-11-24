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
import plot_result.read_csv as data
import linecache
from astropy.timeseries import LombScargle

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
    path_NSC_I='/Volumes/pulsar/SgrA/merge_data/spectra/spectra_p/aprates/'
    path_NSC_G='/Volumes/pulsar/SgrAGRT/merge_data/spectra/aprates/'

    def get_cts_rate(src_ID,obs_ID,path,time,band):
        os.chdir(path)
        ID = int(src_ID)
        ID_note = int(src_ID)
        if str(int(ID))[-3:] == '001' or str(int(ID))[-3:] == '002':
            ID = str(ID)[:-3]
        use_obs = []
        use_obs_time = []
        cts_rate_low = []
        cts_rate_high = []
        cts_rate = []
        for i in range(len(obs_ID)):
            outname = str(ID) + '_' + str(int(obs_ID[i])) + '_out_{0}.par'.format(band)
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
                continue
        return (cts_rate,cts_rate_low,cts_rate_high,use_obs_time,use_obs)

    def get_CR_band(band):
        (CR_I, CR_I_L, CR_I_H, usetime_I, useobs_I) = get_cts_rate(src_ID, obs_ID_I, path_NSC_I, time_I, band)
        (CR_G, CR_G_L, CR_G_H, usetime_G, useobs_G) = get_cts_rate(src_ID, obs_ID_G, path_NSC_G, time_G, band)
        usetime = np.concatenate((usetime_I, usetime_G))
        CR_IG = np.concatenate((CR_I, CR_G))
        CR_low_IG = np.concatenate((CR_I_L, CR_G_L))
        CR_high_IG = np.concatenate((CR_I_H, CR_G_H))
        useobs = np.concatenate((useobs_I, useobs_G))

        return [usetime,CR_IG,CR_low_IG,CR_high_IG]

    def get_LS(time,flux):
        x=time
        y=flux
        freq=1./np.linspace(50, 400, 350)

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
        plt.savefig('/Volumes/pulsar/WR/1671/LS_4_8.eps')
        plt.show()

    #get_LS(get_CR_band('4_8')[0],get_CR_band('4_8')[1])


    CR_H=get_CR_band('4_8')
    CR_S=get_CR_band('2_4')
    print(CR_H[0])
    print(CR_S[0])
    # HR=(CR_H[1]-CR_S[1])/(CR_H[1]+CR_S[1])
    # error=np.sqrt((4*CR_S[1]**2/(CR_H[1]+CR_S[1])**4)*(CR_H[2]-CR_H[1])**2+
    #               (4*CR_H[1]**2/(CR_H[1]+CR_S[1])**4)*(CR_S[2]-CR_S[1])**2)
    HR=CR_H[1]/CR_S[1]

    error=np.sqrt((CR_H[2]-CR_H[1])**2/CR_S[1]**2+(CR_S[2]-CR_S[1])**2*CR_H[1]**2/CR_S[1]**4)

    def pfold(time,P,y):
        def trans(t, p_test, shift=0.0):
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
        plt.xlabel('Phase')
        plt.ylabel('HR')
        plt.errorbar(turns, y,
                     yerr = error,
                     fmt = 'o', capsize = 3, elinewidth = 1, color = 'red', ecolor = 'red')

        plt.savefig('/Volumes/pulsar/WR/1671/pfold_HR_ratio.eps')
        plt.show()

    pfold(CR_H[0],189.,HR)


ID_WR=['1671']
for i in ID_WR:
    get_plot_ap(i)
