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

font1 = {'family': 'Normal',
         'weight': 'normal',
         'size': 18, }

figurepath='/Users/baotong/Desktop/aas/AGN_CDFS_mod1/figure/'

def plot_obs_epoch(epoch_file,color='red',epoch_days='65 days'):
    epoch=np.loadtxt(epoch_file)
    Tstart=epoch[:,0];Tstop=epoch[:,1];exp=epoch[:,3]
    x=(Tstart+Tstop)/2
    error_x=exp
    x=x / 86400 + 2449352.5 - 2400000.5
    error_x=error_x/1000
    y=np.zeros(len(x))+error_x
    # error_y=error_x
    plt.xlabel('MJD', font1)
    plt.ylabel('Exposure time (ks)', font1)
    plt.tick_params(labelsize=15)
    plt.plot([x[0]-50,x[0]-50],[0,150],'--',color=color)
    plt.plot([x[-1]+50, x[-1]+50], [0, 150], '--', color=color)
    # plt.annotate(s='', xy=(1, 1), xytext=(0, 0), arrowprops=dict(arrowstyle='<->'))
    plt.annotate("", xy=(x[0]-70, 150),
                xytext=(x[-1]+70, 150), color=color,
                weight="bold",
                arrowprops=dict(arrowstyle="<->", connectionstyle="arc3", color=color))
    plt.text(x[-1]+50,150,epoch_days,color=color,fontsize=12)

    data1=plt.scatter(x, y, marker='*', s=100, color=color)
    # plt.semilogy()
    # plt.show()
    return data1


if __name__=='__main__':
    data1=[];
    epoch_days=['436 days','45 days','65 days','715 days']
    plt.figure(1,(9,6))
    color=['blue','orange','green','red']
    for i in range(4):
        epoch_file='/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep{0}/epoch_src_495.txt'.format(i+1)
        data1.append(plot_obs_epoch(epoch_file,color[i],epoch_days[i]))
    plt.legend(handles=data1,loc=3,bbox_to_anchor=(0.15,0.1),labels=['Epoch 1','Epoch 2','Epoch 3','Epoch 4'],fontsize=12)
    # ##第二个参数 bbox_to_anchor 被赋予的二元组中，num1 用于控制 legend 的左右移动，值越大越向右
    # ##边移动，num2 用于控制 legend 的上下移动，值越大，越向上移动。用于微调图例的位置。
    plt.xlim(xmax=56900)
    ##为了注释文字不出界
    plt.savefig(figurepath+'Epoch_obs.eps')
    plt.savefig(figurepath + 'Epoch_obs.pdf')
    plt.show()




