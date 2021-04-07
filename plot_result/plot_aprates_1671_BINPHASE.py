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

def read_region(regname):
    reg_file=[]
    with open(regname, 'r') as file_to_read:
        while True:
            lines = file_to_read.readline()  # 整行读取数据
            reg_file.append(lines)
            if not lines:
                break
                pass
    region = reg_file[-2][7:-2]
    reg_x, reg_y, reg_r = [float(i) for i in region.split(',')]
    return [reg_x, reg_y, reg_r]

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

def write_aprate_txt(bin):
    ID = src_ID;
    obsbin = get_obs_in_binP(obsid, time_IG, 189., bin = bin)
    def get_cts_and_exp(obsid,mode,band):
        if mode == 'NSC':
            path = '/Volumes/pulsar/SgrA/merge_data/spectra/spectra_p'
            os.chdir(path)
            if band=='2_4':temp_1='2';temp_2=[137,273]
            elif band=='4_8':temp_1='3';temp_2=[274,548]
            else:print('error')
            pi_file_name = str(ID) + '_' + str(int(obsid)) + '.pi'
            b_pi_file_name = 'b_' + pi_file_name
            expmap_name = 'reproj_evt2_sou_' + str(int(obsid)) + '_i{0}.fits'.format(temp_1)
            reg = read_region(str(int(ID)) + '.reg')
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
                print(src_cts, bkg_cts)
                expmap = fits.open(expmap_name)
                exptime_all = expmap[0].data
                exptime_all = exptime_all.T
                img_x = reg[0]
                img_y = reg[1]
                ##特别注意，这里的expmap，physical坐标与image坐标一致(忘了是为啥了)
                r = reg[2] ** 0.707
                ##为了简便，取个圆的内接正方形对exposure map平均吧
                exp_src = exptime_all[np.arange(int(img_x - r), int(img_x + r))]
                exp_src = exp_src[:, np.arange(int(img_y - r), int(img_y + r))]
                exp_s = np.mean(exp_src)
        if mode == 'NSC_G':
            path = '/Volumes/pulsar/SgrAGRT/merge_data/spectra/'
            os.chdir(path)
            if band=='2_4':temp_1='2';temp_2=[137,273]
            elif band=='4_8':temp_1='3';temp_2=[274,548]
            else:print('error')
            pi_file_name = str(ID) + '_' + str(int(obsid)) + '.pi'
            b_pi_file_name = 'b_' + pi_file_name
            expmap_name = 'reproj_evt2_sou_' + str(int(obsid)) + '_i{0}.fits'.format(temp_1)
            reg = read_region(str(int(ID)) + '.reg')
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
                print(src_cts, bkg_cts)
                expmap = fits.open(expmap_name)
                exptime_all = expmap[0].data
                exptime_all = exptime_all.T
                img_x = reg[0]-2896.
                img_y = reg[1]-2896.
                ##特别注意，这里的expmap，physical坐标与image坐标一致(忘了是为啥了)
                r = reg[2] ** 0.707
                ##为了简便，取个圆的内接正方形对exposure map平均吧
                exp_src = exptime_all[np.arange(int(img_x - r), int(img_x + r))]
                exp_src = exp_src[:, np.arange(int(img_y - r), int(img_y + r))]
                exp_s = np.mean(exp_src)

        return (src_cts,bkg_cts,exp_s,backscale,b_backscale)

    for i in range(len(obsbin)):
        src_cts_soft = 0;
        bkg_cts_soft = 0;
        exp_soft = 0;
        src_cts_hard = 0;
        bkg_cts_hard = 0;
        exp_hard = 0;
        backscale_soft=0;
        b_backscale_soft=0;
        backscale_hard=0;
        b_backscale_hard=0;
        for obs in obsbin[i]:
            if obs in obs_ID_I:
                aprates_temp_info_soft=get_cts_and_exp(obs,'NSC','2_4')
                aprates_temp_info_hard = get_cts_and_exp(obs, 'NSC', '4_8')
            elif obs in obs_ID_G:
                aprates_temp_info_soft=get_cts_and_exp(obs,'NSC_G','2_4')
                aprates_temp_info_hard = get_cts_and_exp(obs, 'NSC_G', '4_8')
            else:print('error')
            src_cts_soft+=aprates_temp_info_soft[0]
            bkg_cts_soft += aprates_temp_info_soft[1]
            exp_soft+=aprates_temp_info_soft[2]
            backscale_soft+=aprates_temp_info_soft[3]
            b_backscale_soft+=aprates_temp_info_soft[4]

            src_cts_hard+=aprates_temp_info_hard[0]
            bkg_cts_hard += aprates_temp_info_hard[1]
            exp_hard+=aprates_temp_info_hard[2]
            backscale_hard+=aprates_temp_info_hard[3]
            b_backscale_hard+=aprates_temp_info_hard[4]

        os.chdir('/Volumes/pulsar/WR/1671/')

        aprates_text = 'aprates n={0} m={1} A_s={2} A_b={3} alpha=0.9 beta=0.02 T_s=1 ' \
                       'E_s={4} eng_s=1 flux_s=1 T_b=1 E_b={5} eng_b=1 flux_b=1 clobber=yes ' \
                       'outfile={6} conf=0.68'.format(src_cts_soft, bkg_cts_soft, backscale_soft, b_backscale_soft, exp_soft, exp_soft,
                                                        str(int(ID)) + '_bin' + str(i) + '_out_2_4.par')
        with open('aprates/bin{0}/'.format(bin) + 'run_' + str(int(ID)) + '_bin' + str(i)  + '_2_4.e', 'w+') as f:
            f.writelines(aprates_text)

        aprates_text = 'aprates n={0} m={1} A_s={2} A_b={3} alpha=0.9 beta=0.02 T_s=1 ' \
                       'E_s={4} eng_s=1 flux_s=1 T_b=1 E_b={5} eng_b=1 flux_b=1 clobber=yes ' \
                       'outfile={6} conf=0.68'.format(src_cts_hard, bkg_cts_hard, backscale_hard, b_backscale_hard, exp_hard, exp_hard,
                                                        str(int(ID)) + '_bin' + str(i) + '_out_4_8.par')
        with open('aprates/bin{0}/'.format(bin) + 'run_' + str(int(ID)) + '_bin' + str(i)  + '_4_8.e', 'w+') as f:
            f.writelines(aprates_text)



#write_aprate_txt(bin=20)

def plot_aprates(bin):
    ID = src_ID;
    path='/Volumes/pulsar/WR/1671/aprates/bin{0}/'.format(bin)
    os.chdir(path)
    cts_rate_low_s = [];cts_rate_low_h = [];
    cts_rate_high_s = [];cts_rate_high_h = [];
    cts_rate_s = [];cts_rate_h = [];
    for i in range(bin):
        filename_s=str(ID)+'_bin{0}'.format(str(i))+'_out_2_4.par'
        filename_h = str(ID) + '_bin{0}'.format(str(i)) + '_out_4_8.par'
        if os.path.exists(filename_s):
            a = get_line_context(path+filename_s, 15)[18:-3]
            a_low = get_line_context(path+filename_s, 16)[25:-3]
            a_high = get_line_context(path+filename_s, 17)[25:-3]

            b = get_line_context(path+filename_h, 15)[18:-3]
            b_low = get_line_context(path+filename_h, 16)[25:-3]
            b_high = get_line_context(path+filename_h, 17)[25:-3]
        else:
            a=a_low=a_low=b=b_low=b_high=1e-5


        cts_rate_s.append(float(a))
        cts_rate_low_s.append(float(a_low))
        cts_rate_high_s.append(float(a_high))

        cts_rate_h.append(float(b))
        cts_rate_low_h.append(float(b_low))
        cts_rate_high_h.append(float(b_high))
    CR_S=np.array([np.array(cts_rate_low_s),np.array(cts_rate_s),np.array(cts_rate_high_s)])
    CR_H = np.array([np.array(cts_rate_low_h), np.array(cts_rate_h), np.array(cts_rate_high_h)])
    HR = CR_H[1] / CR_S[1]
    print(HR)
    error = np.sqrt((CR_H[2] - CR_H[1]) ** 2 / CR_S[1] ** 2 + (CR_S[2] - CR_S[1]) ** 2 * CR_H[1] ** 2 / CR_S[1] ** 4)

    #error=np.zeros(len(HR))+0.01
    # error=CR_H[1]*(CR_S[2] - CR_S[1]) /CR_S[1] ** 2+(CR_H[2] - CR_H[1]) /CR_S[1]
    # HR=(CR_H[1]-CR_S[1])/(CR_H[1]+CR_S[1])
    # error=np.sqrt((4*CR_S[1]**2/(CR_H[1]+CR_S[1])**4)*(CR_H[2]-CR_H[1])**2+
    #               (4*CR_H[1]**2/(CR_H[1]+CR_S[1])**4)*(CR_S[2]-CR_S[1])**2)
    x=np.linspace(0,1./bin*(bin-1),bin)+0.5/bin
    plt.xlabel('phase')
    plt.ylabel('H/S')
    plt.errorbar(x, HR,
                 yerr = error,
                 fmt = 'o', capsize = 3, elinewidth = 1, color = 'red', ecolor = 'red')

    plt.savefig('phase_HS_bin{0}.eps'.format(bin))
    plt.show()
plot_aprates(bin=20)