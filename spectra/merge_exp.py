#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import string
# import correct as correct
from datetime import datetime
from scipy.interpolate import lagrange
from scipy import interpolate
from scipy.fftpack import fft, ifft
import scipy.signal as ss
import random
import pandas as pd

obs_ID_I = np.loadtxt('SgrA_I_epoch.txt')[:, 2]
obs_time_I = (np.loadtxt('SgrA_I_epoch.txt')[:, 0] + np.loadtxt('SgrA_I_epoch.txt')[:, 1]) / 2
time_I = obs_time_I / 86400 + 2449352.5 - 2400000.5

obs_ID_G = np.loadtxt('SgrA_G_epoch.txt')[:, 2]
obs_time_G = (np.loadtxt('SgrA_G_epoch.txt')[:, 0] + np.loadtxt('SgrA_G_epoch.txt')[:, 1]) / 2
time_G = obs_time_G / 86400 + 2449352.5 - 2400000.5

path_I='/Volumes/pulsar/SgrA/merge_data/spectra/spectra_p/'
path_G='/Volumes/pulsar/SgrAGRT/merge_data/spectra/'
reg=[4201.0,4563.0,5.509019423583278]
def merge_exp():
    EXP_I=[]
    for obs in obs_ID_I:
        expmap=fits.open(path_I+'reproj_evt2_sou_{0}_i4.fits'.format(int(obs)))

        exptime_all = expmap[0].data
        exptime_all = exptime_all.T
        img_x = reg[0]
        img_y = reg[1]
        r = reg[2] ** 0.707
        ##为了简便，取个圆的内接正方形做平均吧
        exp_src = exptime_all[np.arange(int(img_x - r), int(img_x + r))]
        exp_src = exp_src[:, np.arange(int(img_y - r), int(img_y + r))]
        exp_s = np.mean(exp_src)
        EXP_I.append(exp_s)
        print('run')

    EXP_I=np.sum(np.array(EXP_I))

    EXP_G=[]
    for obs in obs_ID_G:
        expmap=fits.open(path_G+'reproj_evt2_sou_{0}_i4.fits'.format(int(obs)))

        exptime_all = expmap[0].data
        exptime_all = exptime_all.T
        img_x = reg[0] - 2896
        img_y = reg[1] - 2896
        r = reg[2] ** 0.707
        ##为了简便，取个圆的内接正方形做平均吧
        exp_src = exptime_all[np.arange(int(img_x - r), int(img_x + r))]
        exp_src = exp_src[:, np.arange(int(img_y - r), int(img_y + r))]
        exp_s = np.mean(exp_src)
        EXP_G.append(exp_s)
        print('run')

    EXP_G=np.sum(np.array(EXP_G))

    print('EXP_I={0}'.format(EXP_I))
    print('EXP_G={0}'.format(EXP_G))

def get_ratio_from2and3():
    def get_exp_mean(reg,expfile,mode):
        expmap = fits.open(expfile)
        exptime_all = expmap[0].data
        exptime_all = exptime_all.T
        if mode=='G':
            img_x = reg[0] - 2896
            img_y = reg[1] - 2896
        elif mode=='I':
            img_x = reg[0]
            img_y = reg[1]
        r = reg[2] ** 0.707
        ##为了简便，取个圆的内接正方形做平均吧
        exp_src = exptime_all[np.arange(int(img_x - r), int(img_x + r))]
        exp_src = exp_src[:, np.arange(int(img_y - r), int(img_y + r))]
        exp_s = np.mean(exp_src)
        return exp_s
    path='/Volumes/pulsar/WR/1671/HR/'
    os.chdir(path)
    exp_2I=get_exp_mean(reg,'SgrA_e2_exp.fits','I')
    exp_3I = get_exp_mean(reg, 'SgrA_e3_exp.fits', 'I')
    exp_2G = get_exp_mean(reg, 'SgrAGRT_2_exp.fits', 'G')
    exp_3G = get_exp_mean(reg, 'SgrAGRT_3_exp.fits', 'G')
    print('exp_2I={0}'.format(exp_2I))
    print('exp_3I={0}'.format(exp_3I))
    print('exp_2G={0}'.format(exp_2G))
    print('exp_3G={0}'.format(exp_3G))
    print('ratio_2={0}'.format(exp_2G/exp_2I))
    print('ratio_3={0}'.format(exp_3G/exp_3I))

get_ratio_from2and3()



