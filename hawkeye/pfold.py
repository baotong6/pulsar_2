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
from astropy.stats import poisson_conf_interval
import scipy
import hawkeye.timing_funcs as funcs

font1 = {'family': 'Normal',
         'weight': 'normal',
         'size': 18, }

def trans(t, p_test, shift):
    ti = t
    v = 1.0 / p_test
    turns = v * ti
    turns += shift
    # 初始相位
    for i in range(len(turns)):
        turns[i] = turns[i] - int(turns[i])
    return turns

def phase_fold(time,epoch_info,p_test,outpath,bin=20,net_percent=0.9,shift=0.0,label='test',save=False,show=True):
    turns=trans(time,p_test,shift)
    loc=np.zeros(bin)
    for index in turns:
        loc[int(index*bin)] += 1
    AM=1-min(loc)/max(loc)
    A0=AM/(2-AM)
    print('A0={0}'.format(A0))
    x = np.array([(i / bin + 0.5 / bin) for i in range(bin)])
    src_bkg = 1 - net_percent
    bkg_y = len(time) * src_bkg/bin
    b_1sigma = poisson_conf_interval(bkg_y, interval='frequentist-confidence').T
    bkg_y_low=b_1sigma[0];bkg_y_high=b_1sigma[1]
    fig=plt.figure(1,(10,7.5))
    ax1 = fig.add_subplot(111)
    bkg_x = [0, 2]
    plt.fill_between(bkg_x, bkg_y_low, bkg_y_high,facecolor = 'green', alpha = 0.5)
    x2 = np.concatenate((x, x + 1))
    y2 = np.concatenate((loc, loc))
    T_in_perbin = funcs.get_T_in_mbins(epoch_info, 2 * np.pi / p_test, bin, shift * 2 * np.pi)

    correct_gap = T_in_perbin / (sum(T_in_perbin) / len(T_in_perbin))
    print('correct_gap=',correct_gap)
    y2 /= np.concatenate((correct_gap, correct_gap))
    y2_err = np.array(poisson_conf_interval(y2, interval='frequentist-confidence'))
    y2_err[0] = y2 - y2_err[0]
    y2_err[1] = y2_err[1] - y2

    plt.title("#{0} P={1:.2f},C={2}".format(label, p_test, str(len(time))), fontsize=18)
    plt.xlabel('phase', font1)
    plt.ylabel('counts/bin', font1)
    plt.tick_params(labelsize=18)
    plt.ylim(0, (np.max(y2) + np.max(y2) ** 0.5) * 1.05)
    plt.step(np.concatenate(([0], x2)), np.concatenate(([y2[0]], y2)), color='red')
    plt.errorbar(x2 - 0.5 / bin, y2, yerr=y2_err, fmt='.', capsize=1, elinewidth=1, ecolor='red')

    ax2 = ax1.twinx()
    yhigh = (np.max(y2) + np.max(y2) ** 0.5) * 1.05 / np.mean(y2)
    ax2.set_ylabel('Normalized flux', font1)
    ax2.plot([0, 2], [1.0, 1.0], '--', color='green')
    ax2.set_ylim([0, yhigh])
    ax2.tick_params(labelsize=18)
    if save:plt.savefig(outpath + 'pfold_lc_{0}.eps'.format(label),bbox_inches='tight', pad_inches=0.0)
    if show:plt.show()
    else:plt.close()

