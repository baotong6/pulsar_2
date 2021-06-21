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
import useful_functions as func

font1 = {'family': 'Normal',
         'weight': 'normal',
         'size': 18, }

def phase_fold(data_file,bkg_file,epoch_file,p_test,bin,shift,label):
    src_evt=np.loadtxt(data_file)
    bkg_evt=np.loadtxt(bkg_file)
    epoch=np.loadtxt(epoch_file)
    use_id=epoch[:,2]
    (src_evt,bkg_evt)=func.filter_obs(src_evt,bkg_evt,use_id)
    time=src_evt[:,0];energy=src_evt[:,1]
    time_bkg= bkg_evt[:, 0];energy_bkg= bkg_evt[:, 1]
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
    for index in range(len(turns)):
        loc[int(turns[index]*bin)] += 1
    for index_b in range(len(turns_b)):
        loc_b[int(turns_b[index_b]*bin)]+=1
    x = np.array([(i / bin + 0.5 / bin) for i in range(bin)])

    bkg_y=np.concatenate((loc_b,loc_b))
    bkg_y /= bin
    bkg_y_err = np.array(poisson_conf_interval(bkg_y, interval='frequentist-confidence'))
    bkg_y_err[0]=bkg_y-bkg_y_err[0]
    bkg_y_err[1]=bkg_y_err[1]-bkg_y



    fig=plt.figure(1,(10,7.5))
    x2=np.concatenate((x,x+1))
    y2=np.concatenate((loc,loc))

    # correct_gap = T_in_perbin / (sum(T_in_perbin) / len(T_in_perbin))
    # y2 /= np.concatenate((correct_gap, correct_gap))
    y2_err=np.array(poisson_conf_interval(y2,interval='frequentist-confidence'))
    y2_err[0]=y2-y2_err[0]
    y2_err[1]=y2_err[1]-y2

    #label='F1'
    plt.title("#{0} P={1:.2f},C={2}".format(label,p_test,str(len(time))), fontsize = 18)
    plt.xlabel('phase',font1)
    plt.ylabel('counts/bin',font1)
    plt.tick_params(labelsize = 18)
    plt.ylim(0,(np.max(y2)+np.max(y2)**0.5)*1.05)
    plt.step(x2,y2,color='red')
    #plt.step(x2,y3,'--',color='green')
    print(np.size(y2))
    print(np.size(y2_err))
    plt.errorbar(x2 - 0.5 / bin, bkg_y, yerr = bkg_y_err, fmt = '.', capsize = 1, elinewidth = 0.5, ecolor = 'green')
    plt.errorbar(x2 - 0.5 / bin, y2, yerr = y2_err, fmt = '.', capsize = 1, elinewidth = 1, ecolor = 'red')

    # plt.savefig(pathout+'pfold_lc_{0}.eps'.format(label))
    plt.show()
    #plt.close()
if __name__=='__main__':
    path='/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep3/'
    id='210'
    data_file=path+'{0}.txt'.format(id)
    bkg_file=path+'{0}_bkg.txt'.format(id)
    epoch_file=path+'epoch_src_{0}.txt'.format(id)
    p_test=1614.96326
    phase_fold(data_file, bkg_file, epoch_file, p_test, bin=20, shift=0.4,label=id)