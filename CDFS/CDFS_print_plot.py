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

font1 = {'family': 'Normal',
         'weight': 'normal',
         'size': 16, }

figurepath='/Users/baotong/Desktop/aas/AGN_CDFS/figure/'

def plot_GL_result_of_CDFS():
    threshold=0.99
    plt.figure(1,(15,10))
    plt.ylim(0,1)
    def get_result_GL(k):
        path='/Users/baotong/Desktop/CDFS/result_all_obs_0.5_8_ep{0}/'.format(k)
        result_file_1h=np.loadtxt(path+'all_result_1h.txt')
        result_file_3h=np.loadtxt(path+'all_result_1h_3h.txt')
        id_1h=result_file_1h[:,0]
        prob_1h=result_file_1h[:,2]
        period_1h=result_file_1h[:,4]

        id_3h=result_file_3h[:,0]
        prob_3h=result_file_3h[:,2]
        period_3h=result_file_3h[:,4]


        prob=np.zeros(1055)
        period=np.zeros(1055)

        for i in range(len(1055)):
            if i in id_1h:
                if i in id_3h:
                    if prob_1h[np.where(id_1h)==i] >= prob_3h[np.where(id_3h==i)]:
                        prob[i]=prob_1h[np.where(id_1h)==i]
                        period[i]=period_3h[np.where(id_1h)==i]
                    else:
                        prob[i] = prob_3h[np.where(id_3h) == i]
                        period[i] = period_3h[np.where(id_3h) == i]
                else:
                    prob[i] = prob_1h[np.where(id_1h) == i]
                    period[i] = period_3h[np.where(id_1h) == i]
            else:continue

        prob_out=prob[np.where(prob>threshold)]
        period_out=period[np.where(prob>threshold)]
        prob_in=prob[np.where(prob<=threshold)]
        period_in=period[np.where(prob<=threshold)]

        label=220+k
        plt.subplot(label)
        plt.title('Epoch{0}'.format(k), font1)
        plt.scatter(period_out, prob_out, marker='^',s=80, color='purple')
        plt.scatter(period_in,prob_in,marker='.',s=50,color='black')
        plt.plot([0, 3000], [threshold, threshold], '--', color='black')
        plt.plot([707, 707], [0, threshold], '--', color='r')
        plt.plot([1000, 1000], [0, threshold], '--', color='green')
        plt.plot([707 * 2, 707 * 2], [0, threshold], '--', color='r')
        plt.plot([1000 * 2, 1000 * 2], [0, threshold], '--', color='green')
        plt.plot([707 * 3, 707 * 3], [0, threshold], '--', color='r')
        plt.plot([1000 * 3, 1000 * 3], [0, threshold], '--', color='green')
        plt.ylabel('Probability', font1)
        plt.tick_params(labelsize=16)
        plt.semilogx()

        if k>2:plt.xlabel('Period (s)',font1)
        return (prob,period)

    get_result_GL(1)
    get_result_GL(2)
    get_result_GL(3)
    get_result_GL(4)

    plt.savefig(figurepath+'GL_res.eps',bbox_inches='tight',pad_inches=0.0)
    plt.show()
plot_GL_result_of_CDFS()

def plot_LS_result_of_CDFS():
    threshold=0.9973
    plt.figure(1,(15,10))
    plt.ylim(0,1)
    def get_result_LS(k):
        path = '/Users/baotong/Desktop/CDFS/fig_LS_ep{0}_ovsamp_5_baluev/'.format(k)
        LS_result = np.loadtxt(path + 'LS_result_{0}.txt'.format(k))
        FP = LS_result[:, 1]
        prob=1-FP
        period = LS_result[:, 2]

        prob_out=prob[np.where(((period-200)>10)&(prob>threshold))]
        period_out=period[np.where(((period-200)>10)&(prob>threshold))]
        print(period_out)

        prob_in=prob[np.where(prob<=threshold)]
        period_in=period[np.where(prob<=threshold)]

        label=220+k
        plt.subplot(label)
        plt.title('Epoch{0}'.format(k), font1)

        plt.scatter(period_out, prob_out, marker='v',s=80, color='purple')
        plt.scatter(period_in,prob_in,marker='.',s=50,color='black')
        plt.plot([0, 10000], [threshold, threshold], '--', color='black')
        plt.plot([707, 707], [0, threshold], '--', color='r')
        plt.plot([1000, 1000], [0, threshold], '--', color='green')
        plt.plot([707 * 2, 707 * 2], [0, threshold], '--', color='r')
        plt.plot([1000 * 2, 1000 * 2], [0, threshold], '--', color='green')
        plt.plot([707 * 3, 707 * 3], [0, threshold], '--', color='r')
        plt.plot([1000 * 3, 1000 * 3], [0, threshold], '--', color='green')
        plt.ylabel('1-FAP', font1)
        plt.tick_params(labelsize=16)
        plt.semilogx()

        if k>2:plt.xlabel('Period (s)',font1)
        return (prob,period)
    get_result_LS(1)
    get_result_LS(2)
    get_result_LS(3)
    get_result_LS(4)

    plt.savefig(figurepath+'LS_res.eps',bbox_inches='tight',pad_inches=0.0)

    plt.show()

# plot_LS_result_of_CDFS()