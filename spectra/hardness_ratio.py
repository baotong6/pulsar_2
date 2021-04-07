#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import string
path='/Users/baotong/Desktop/period/txt_all_obs_IG/'
src_ID='1671'
obsid=np.loadtxt(path+'epoch_src_{0}.txt'.format(src_ID))[:,2]
obsid_1671_05_1=np.array([15611., 15612. , 2954.  ,2943. , 3663. , 3392.  ,3393.,  3665.,  3549.,  4683.,
  4684.  ,5360.  ,6113. , 5950. , 5951. , 5952.,  5953.  ,5954.,  6642. , 6363.,
  6643. , 6644.  ,6645. , 7554. , 7555.,  7557.,  7558. , 9174. , 9173. ,13016.,
 13017. ,14941. ,14942., 13856. ,13857., 13854. ,14413. ,13855. ,14414. ,13847.,
 14427. ,13848. ,13849., 13846., 14438., 13845. ,14462. ,14463. ,13851. ,15568.,
 13843. ,15570. ,14468.])
obsid_1671_0_05=np.setdiff1d(obsid,obsid_1671_05_1)
obsinfo = np.loadtxt(path + 'epoch_src_{0}.txt'.format(src_ID))
obs_time = (obsinfo[:, 0] + obsinfo[:, 1]) / 2
time = obs_time / 86400 + 2449352.5 - 2400000.5

obs_ID_I = np.loadtxt('SgrA_I_epoch.txt')[:, 2]
obs_time_I = (np.loadtxt('SgrA_I_epoch.txt')[:, 0] + np.loadtxt('SgrA_I_epoch.txt')[:, 1]) / 2
time_I = obs_time_I / 86400 + 2449352.5 - 2400000.5

obs_ID_G = np.loadtxt('SgrA_G_epoch.txt')[:, 2]
obs_time_G = (np.loadtxt('SgrA_G_epoch.txt')[:, 0] + np.loadtxt('SgrA_G_epoch.txt')[:, 1]) / 2
time_G = obs_time_G / 86400 + 2449352.5 - 2400000.5

ratio_2=0.6444334590223477
ratio_3=0.9579565096645598
def get_HR_txt():
    src_cts_soft=[]
    bkg_cts_soft=[]
    src_cts_hard=[]
    bkg_cts_hard=[]
    ratio_area=[]

    path_spec='/Volumes/pulsar/WR/1671/HR/'
    for obs in obsid:
        pi_file_name = str(src_ID) + '_' + str(int(obs)) + '.pi'
        b_pi_file_name = 'b_' + pi_file_name
        os.chdir(path_spec)
        if os.path.exists(pi_file_name) and os.path.exists(b_pi_file_name):
            pi = fits.open(pi_file_name)
            b_pi = fits.open(b_pi_file_name)
            backscale = pi[1].header['BACKSCAL']
            b_backscale = b_pi[1].header['BACKSCAL']
            backscale *= 8192 ** 2
            b_backscale *= 8192 ** 2

            pi_data = pi[1].data
            b_pi_data = b_pi[1].data
            ratio_area.append(b_backscale/backscale)
            if obs in obs_ID_I:
                src_cts_soft.append(ratio_2*sum(pi_data['COUNTS'][137:273]))
                bkg_cts_soft.append(ratio_2*sum(b_pi_data['COUNTS'][137:273]))
                src_cts_hard.append(ratio_3*sum(pi_data['COUNTS'][274:548]))
                bkg_cts_hard.append(ratio_3*sum(b_pi_data['COUNTS'][274:548]))
            elif obs in obs_ID_G:
                src_cts_soft.append(sum(pi_data['COUNTS'][137:273]))
                bkg_cts_soft.append(sum(b_pi_data['COUNTS'][137:273]))
                src_cts_hard.append(sum(pi_data['COUNTS'][274:548]))
                bkg_cts_hard.append(sum(b_pi_data['COUNTS'][274:548]))
            else:
                print('error')

    all_cts_info=np.column_stack((src_cts_soft,src_cts_hard,bkg_cts_soft,bkg_cts_hard,ratio_area,ratio_area))
    np.savetxt(path_spec+'HR_1671_info_com.txt',all_cts_info,fmt='%5d %5d %5d %5d %5.2f %5.2f')
#get_HR_txt()

def plot_HR(src_ID='1671',P=189.):
    obsinfo = np.loadtxt(path + 'epoch_src_{0}.txt'.format(src_ID))
    obs_time = (obsinfo[:, 0] + obsinfo[:, 1]) / 2
    time = obs_time / 86400 + 2449352.5 - 2400000.5

    obs_ID_I = np.loadtxt('SgrA_I_epoch.txt')[:, 2]
    obs_time_I = (np.loadtxt('SgrA_I_epoch.txt')[:, 0] + np.loadtxt('SgrA_I_epoch.txt')[:, 1]) / 2
    time_I = obs_time_I / 86400 + 2449352.5 - 2400000.5

    obs_ID_G = np.loadtxt('SgrA_G_epoch.txt')[:, 2]
    obs_time_G = (np.loadtxt('SgrA_G_epoch.txt')[:, 0] + np.loadtxt('SgrA_G_epoch.txt')[:, 1]) / 2
    time_G = obs_time_G / 86400 + 2449352.5 - 2400000.5

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


    path_spec = '/Volumes/pulsar/WR/1671/HR/'
    filename='HR_1671_result_com.txt'
    file=np.loadtxt(path_spec+filename)
    HR=file[:,0]
    HR_low=file[:,1]
    HR_high=file[:,2]

    plt.scatter(time,HR)
    plt.show()

    plt.figure(1,(8,6))

    plt.xlabel('phase')
    plt.ylabel('HR')
    ratio=1.0
    for i in range(len(turns)):
        if time[i] in time_I:
            plt.errorbar(turns[i], HR[i], yerr = [[HR[i] - HR_low[i]], [HR_high[i] - HR[i]]],
                         fmt = 'o', capsize = 3, elinewidth = 1, color = 'red', ecolor = 'red')
            plt.errorbar(turns[i] + 1., HR[i], yerr =  [[HR[i] - HR_low[i]], [HR_high[i] - HR[i]]],
                         fmt = 'o', capsize = 3, elinewidth = 1, color = 'red', ecolor = 'red')
        else:
            plt.errorbar(turns[i], HR[i]*ratio, yerr =  [[HR[i]*ratio - HR_low[i]*ratio], [HR_high[i]*ratio - HR[i]*ratio]],
                         fmt = 'o', capsize = 3, elinewidth = 1, color = 'green', ecolor = 'green')
            plt.errorbar(turns[i] + 1., HR[i]*ratio, yerr =  [[HR[i]*ratio - HR_low[i]*ratio], [HR_high[i]*ratio - HR[i]*ratio]],
                         fmt = 'o', capsize = 3, elinewidth = 1, color = 'green', ecolor = 'green')

    plt.savefig('/Volumes/pulsar/WR/1671/HR_com.eps')
    plt.show()

plot_HR()