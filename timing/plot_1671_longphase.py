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

ratio_2=0.6444334590223477
ratio_3=0.9579565096645598

def get_line_context(file_path, line_number):
    return linecache.getline(file_path, line_number).strip()

src_ID='1671'

obs_ID_I = np.loadtxt('SgrA_I_epoch.txt')[:, 2]
obs_time_I = (np.loadtxt('SgrA_I_epoch.txt')[:, 0] + np.loadtxt('SgrA_I_epoch.txt')[:, 1]) / 2
time_I = obs_time_I / 86400 + 2449352.5 - 2400000.5

obs_ID_G = np.loadtxt('SgrA_G_epoch.txt')[:, 2]
obs_time_G = (np.loadtxt('SgrA_G_epoch.txt')[:, 0] + np.loadtxt('SgrA_G_epoch.txt')[:, 1]) / 2
time_G = obs_time_G / 86400 + 2449352.5 - 2400000.5

path='/Users/baotong/Desktop/period/txt_all_obs_IG/'
obsid=np.loadtxt(path+'epoch_src_{0}.txt'.format(src_ID))[:,2]
obs_time=(np.loadtxt(path+'epoch_src_{0}.txt'.format(src_ID))[:,0]+
          np.loadtxt(path+'epoch_src_{0}.txt'.format(src_ID))[:,1])/2
time_IG = obs_time / 86400 + 2449352.5 - 2400000.5


def get_obs_in_binP(obs,time,P,bin):
    def trans(t, p_test, shift = 0.0):
        ti = t
        v = 1.0 / p_test
        turns = v * ti
        turns += shift
        # 初始相位
        for i in range(len(turns)):
            turns[i] = turns[i] - int(turns[i])
        return turns

    turns = trans(time, P)
    obs_bin=[]
    for i in range(int(bin)):
        obs_bin.append(obs[np.where((turns < (i+1)*1./bin) & (turns > i*1./bin))])
    return obs_bin

def write_HR_txt(bin):
    ID = src_ID;
    obsbin = get_obs_in_binP(obsid, time_IG, 189., bin = bin)
    def get_cts_and_exp(obsid,mode,band):
        if mode == 'NSC':
            path = '/Volumes/pulsar/SgrA/merge_data/spectra/spectra_p'
            os.chdir(path)
            if band=='2_4':temp_1='2';temp_2=[137,273];ratio=ratio_2
            elif band=='4_8':temp_1='3';temp_2=[274,548];ratio=ratio_3
            else:print('error')
            pi_file_name = str(ID) + '_' + str(int(obsid)) + '.pi'
            b_pi_file_name = 'b_' + pi_file_name
            if os.path.exists(pi_file_name) and os.path.exists(b_pi_file_name):
                pi = fits.open(pi_file_name)
                b_pi = fits.open(b_pi_file_name)
                backscale = pi[1].header['BACKSCAL']
                b_backscale = b_pi[1].header['BACKSCAL']
                backscale *= 8192 ** 2
                b_backscale *= 8192 ** 2

                pi_data = pi[1].data
                b_pi_data = b_pi[1].data
                src_cts = ratio*sum(pi_data['COUNTS'][temp_2[0]:temp_2[1]])
                bkg_cts = ratio*sum(b_pi_data['COUNTS'][temp_2[0]:temp_2[1]])

        if mode == 'NSC_G':
            path = '/Volumes/pulsar/SgrAGRT/merge_data/spectra/'
            os.chdir(path)
            if band=='2_4':temp_1='2';temp_2=[137,273]
            elif band=='4_8':temp_1='3';temp_2=[274,548]
            else:print('error')
            pi_file_name = str(ID) + '_' + str(int(obsid)) + '.pi'
            b_pi_file_name = 'b_' + pi_file_name
            if os.path.exists(pi_file_name) and os.path.exists(b_pi_file_name):
                pi = fits.open(pi_file_name)
                b_pi = fits.open(b_pi_file_name)
                backscale = pi[1].header['BACKSCAL']
                b_backscale = b_pi[1].header['BACKSCAL']
                backscale *= 8192 ** 2
                b_backscale *= 8192 ** 2

                pi_data = pi[1].data
                b_pi_data = b_pi[1].data
                src_cts = sum(pi_data['COUNTS'][temp_2[0]:temp_2[1]])
                bkg_cts = sum(b_pi_data['COUNTS'][temp_2[0]:temp_2[1]])

        return (src_cts,bkg_cts,backscale,b_backscale)

    src_cts_soft_out=[]
    bkg_cts_soft_out=[]
    src_cts_hard_out=[]
    bkg_cts_hard_out=[]
    ratio_area_out=[]

    for i in range(len(obsbin)):
        src_cts_soft = 0;
        bkg_cts_soft = 0;
        src_cts_hard = 0;
        bkg_cts_hard = 0;
        ratio_area=0;

        for obs in obsbin[i]:
            if obs in obs_ID_I:
                aprates_temp_info_soft = get_cts_and_exp(obs, 'NSC', '2_4')
                aprates_temp_info_hard = get_cts_and_exp(obs, 'NSC', '4_8')
            elif obs in obs_ID_G:
                aprates_temp_info_soft = get_cts_and_exp(obs, 'NSC_G', '2_4')
                aprates_temp_info_hard = get_cts_and_exp(obs, 'NSC_G', '4_8')
            else:
                print('error')
            src_cts_soft += aprates_temp_info_soft[0]
            bkg_cts_soft += aprates_temp_info_soft[1]
            ratio_area = aprates_temp_info_soft[3]/aprates_temp_info_soft[2]

            src_cts_hard += aprates_temp_info_hard[0]
            bkg_cts_hard += aprates_temp_info_hard[1]

        src_cts_soft_out.append(src_cts_soft)
        bkg_cts_soft_out.append(bkg_cts_soft)
        src_cts_hard_out.append(src_cts_hard)
        bkg_cts_hard_out.append(bkg_cts_hard)
        ratio_area_out.append(ratio_area)

    all_cts_info = np.column_stack((src_cts_soft_out, src_cts_hard_out,
                                        bkg_cts_soft_out, bkg_cts_hard_out,
                                        ratio_area_out, ratio_area_out))
    path='/Volumes/pulsar/WR/1671/HR_LONGBIN/'
    np.savetxt(path  + 'HR_1671_long_bin{0}.txt'.format(bin), all_cts_info,
                   fmt = '%5d %5d %5d %5d %5.2f %5.2f')

#write_HR_txt(10)
def plot_HR_result(bin):
    path = '/Volumes/pulsar/WR/1671/HR_LONGBIN/'
    filename='HR_1671_long_bin20_result.txt'
    file=np.loadtxt(path+filename)
    HR=file[:,0]
    HR_low=file[:,1]
    HR_high=file[:,2]

    plt.xlabel('phase')
    plt.ylabel('HR')

    x=np.linspace(0,1-1./bin,bin)+0.5/bin
    plt.errorbar(x,HR,yerr=[HR-HR_low,HR_high-HR],fmt = 'o',
                 capsize = 3, elinewidth = 1, color = 'red', ecolor = 'red')
    plt.savefig(path+'HR_1671_long_bin20.eps')
    plt.show()

plot_HR_result(20)



