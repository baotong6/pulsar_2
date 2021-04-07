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
path='/Volumes/pulsar/WR/1671/txt/'
def write_HR_txt(bin,srctxt,bkgtxt,period=48016.90195):
    src_cts_soft=[]
    bkg_cts_soft=[]
    src_cts_hard=[]
    bkg_cts_hard=[]
    ratio_area=[]

    def trans(t, p_test, shift = 0.0):
        ti = t
        v = 1.0 / p_test
        turns = v * ti
        turns += shift
        # 初始相位
        for i in range(len(turns)):
            turns[i] = turns[i] - int(turns[i])
        turns=np.array(turns)
        return turns

    src_txt=np.loadtxt(path+srctxt)
    bkg_txt=np.loadtxt(path+bkgtxt)
    src_soft_txt=src_txt[np.where((src_txt[:,1]>2000)&(src_txt[:,1]<4000))]
    bkg_soft_txt=bkg_txt[np.where((bkg_txt[:,1]>2000)&(bkg_txt[:,1]<4000))]
    src_hard_txt = src_txt[np.where((src_txt[:, 1] > 4000) & (src_txt[:, 1] < 8000))]
    bkg_hard_txt = bkg_txt[np.where((bkg_txt[:, 1] > 4000) & (bkg_txt[:, 1] < 8000))]
    for i in range(bin):
        src_cts_soft.append(ratio_2*len(src_soft_txt[np.where((trans(src_soft_txt[:,0],period)>i*1./bin)
                 &(trans(src_soft_txt[:,0],period)<(i+1)*1./bin))]))
        src_cts_hard.append(ratio_3*len(src_hard_txt[np.where((trans(src_hard_txt[:, 0], period) > i * 1. / bin)
                                                      & (trans(src_hard_txt[:, 0], period) < (i + 1) * 1. / bin))]))
        bkg_cts_soft.append(ratio_2*len(bkg_soft_txt[np.where((trans(bkg_soft_txt[:, 0],  period) > i * 1. / bin)
                                                      & (trans(bkg_soft_txt[:, 0],  period) < (i + 1) * 1. / bin))]))
        bkg_cts_hard.append(ratio_3*len(bkg_hard_txt[np.where((trans(bkg_hard_txt[:, 0],  period) > i * 1. / bin)
                                                      & (trans(bkg_hard_txt[:, 0],  period) < (i + 1) * 1. / bin))]))

        ratio_area.append(10.79)

    print(src_cts_soft,src_cts_hard,bkg_cts_soft,bkg_cts_hard)
    all_cts_info = np.column_stack((src_cts_soft,src_cts_hard, bkg_cts_soft,bkg_cts_hard, ratio_area, ratio_area))
    np.savetxt(path +'HR_BIN/'+ 'HR_1671_I_bin{0}.txt'.format(bin), all_cts_info, fmt = '%5d %5d %5d %5d %5.2f %5.2f')


#write_HR_txt(20,'1671_src_I.txt','1671_bkg_I.txt',period=48016.90195)

def plot_HR_result(bin):
    path = '/Volumes/pulsar/WR/1671/txt/HR_BIN/'
    def get_info(filename):
        file=np.loadtxt(path+filename)
        HR=file[:,0]
        HR_low=file[:,1]
        HR_high=file[:,2]
        return[HR,HR_low,HR_high]
    temp1=get_info('HR_1671_I_bin20_result.txt')
    temp2 = get_info('HR_1671_G_bin20_result.txt')

    HR=temp1[0];HR_low=temp1[1];HR_high=temp1[2]
    HR_G = temp2[0];HR_low_G= temp2[1];HR_high_G = temp2[2]

    plt.xlabel('phase')
    plt.ylabel('HR')

    x=np.linspace(0,1-1./bin,bin)+0.5/bin
    plt.errorbar(x,HR,yerr=[HR-HR_low,HR_high-HR],fmt = 'o',
                 capsize = 3, elinewidth = 1, color = 'red', ecolor = 'red')
    plt.errorbar(x,HR_G,yerr=[HR_G-HR_low_G,HR_high_G-HR_G],fmt = 'o',
                 capsize = 3, elinewidth = 1, color = 'green', ecolor = 'green')
    plt.legend(['ACIS-I','ACIS-S'])
    plt.title('P=48016.90')
    plt.savefig(path+'HR_1671_bin20.eps')

    plt.show()

plot_HR_result(20)



