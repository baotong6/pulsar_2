'''
Author: baotong && baotong@smail.nju.edu.cn
Date: 2024-03-27 20:48:42
LastEditors: baotong && baotong@smail.nju.edu.cn
LastEditTime: 2024-03-28 10:02:13
FilePath: /pulsar/NSC/plot_HR67.py
Description: 

Copyright (c) 2024 by baotong, All Rights Reserved. 
'''
#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
from scipy.interpolate import lagrange
from scipy import optimize as op
from scipy.optimize import curve_fit
import pandas as pd
from astropy.stats import poisson_conf_interval
from astropy.coordinates import SkyCoord
import astropy.units as u
import hawkeye as hawk
from matplotlib.gridspec import GridSpec

ra_center=266.4166667;dec_center= -29.0077778
cat=fits.open('/Users/baotong/Desktop/period/'+'zhu18_3.fits')
ra = cat[1].data['_RAJ2000']
dec = cat[1].data['_DEJ2000']
srcID_list=np.arange(1,len(ra)+1,1)
target_coord = SkyCoord(ra=ra_center*u.deg, dec=dec_center*u.deg, frame='icrs')
catalog_coords =SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame='icrs')
distances = target_coord.separation(catalog_coords).to(u.arcsec).value

def gaussian(x, amplitude,mu, sigma):
    return amplitude * np.exp(-((x - mu) / sigma) ** 2 / 2)


def triple_gaussian(x, a1, mu1, sigma1, a2, mu2, sigma2, a3, mu3, sigma3):
    return (a1*np.exp(-(x-mu1)**2/(2*sigma1**2)) +
            a2*np.exp(-(x-mu2)**2/(2*sigma2**2)) +
            a3*np.exp(-(x-mu3)**2/(2*sigma3**2)))

def plot_HR67_hist(save=0,show=1):
    color_list = ['b']
    HR=[];HRl=[];HRu=[];F=[];Fl=[];Fu=[]
    file = pd.read_csv(f'/Users/baotong/Desktop/period/HR/S26_H67_evt/all_HR_sup.csv')
    pars=np.loadtxt(f'/Users/baotong/Desktop/period/HR/S26_H67_evt/all_inputpars.txt')
    S=pars[:,1]-pars[:,3]/pars[:,5];H=pars[:,2]-pars[:,4]/pars[:,6]
    seq=pars[:,0]
    seq=seq.astype('int')
    S[np.isnan(S)]=0;H[np.isnan(H)]=0
    HR=file['HR'];HRl=file['HRl'];HRu=file['HRu']
    SH=H+S
    SH_1sigma = poisson_conf_interval(SH, interval='frequentist-confidence').T
    SHl = SH_1sigma[:,0];SHu = SH_1sigma[:,1]
    ##==plot HR-histogram==##
    goodindex=np.where(SH>100)[0]
    outindex=np.intersect1d(goodindex,np.where(HR>-0.75)[0])
    print(seq[outindex])
    print(len(goodindex),len(outindex))
    # 绘制拟合曲线和原始直方图

    plt.figure(1,(7,6))
    temphist=plt.hist(HR[goodindex],bins=80,histtype='step',linewidth=3)
    histvalue=temphist[0]
    bin_centers=(temphist[1][:-1]+temphist[1][1:])/2
    # ----- 进行三个高斯函数的拟合 ---- #
    initial_guess = [20, -0.95, 0.05, 60, -0.67, 0.05, 20, -0.4, 0.05]  # 初始猜测值
    bounds = ([-np.inf, -0.95, 0, -np.inf, -0.8, 0, -np.inf, -0.4,0],
          [np.inf, -0.8, 0.2, np.inf, -0.5, 0.2, np.inf, -0.1,0.2])  # 将mu_1固定在10附近
    params, _ = curve_fit(triple_gaussian, bin_centers, histvalue, p0=initial_guess,bounds=bounds,maxfev=800000)
    bin_plot=np.linspace(bin_centers[0],bin_centers[-1],1000)
    plt.plot(bin_plot, triple_gaussian(bin_plot, *params), color='red', label='Fit',lw=3)
    plt.plot(bin_plot, gaussian(bin_plot, params[0], params[1], params[2]), '--', color='cyan',label='Component 1')
    plt.plot(bin_plot, gaussian(bin_plot, params[3], params[4], params[5]), '--', color='cyan',label='Component 2')
    plt.plot(bin_plot, gaussian(bin_plot, params[6], params[7], params[8]), '--', color='cyan',label='Component 3')


    print('params=',params)
    plt.xticks(size=18)
    plt.yticks(size=18)
    plt.xlabel('HR', hawk.font1)
    plt.ylabel('Number of sources per bin',hawk.font1)


    # 创建图形和网格布局
    fig = plt.figure(figsize=(12, 6))
    gs = GridSpec(1, 2, figure=fig)

    # 散点图子图
    ax_scatter = fig.add_subplot(gs[0, 0])
    scatter_good = ax_scatter.scatter(HR[goodindex], distances[goodindex], marker='+', s=50, label='Good Index')
    scatter_out = ax_scatter.scatter(HR[outindex], distances[outindex], marker='+', s=50, color='cyan', label='Out Index')
    ax_scatter.set_xlabel('HR')
    ax_scatter.set_ylabel('Distances')
    ax_scatter.legend()

    # 水平直方图子图
    ax_hist = fig.add_subplot(gs[0, 1], sharey=ax_scatter)
    ax_hist.hist(distances[goodindex], bins=20, alpha=0.5,histtype='step', orientation='horizontal', linestyle='--',label='Good Index')
    ax_hist.hist(distances[outindex], bins=20, alpha=0.5, histtype='step',orientation='horizontal', linestyle='--',color='cyan', label='Out Index')
    ax_hist.set_xlabel('Frequency')
    ax_hist.legend()

    # 调整子图之间的间距
    plt.tight_layout()

    # 显示图形
    plt.show()

# plt.figure(2)
# plt.scatter(HR[goodindex],distances[goodindex],marker='+',s=50)
# plt.scatter(HR[outindex], distances[outindex], marker='+', s=50,color='cyan')
    # for k in range(len(outindex)):
    #     plt.text(HR[outindex[k]], distances[outindex[k]], s=seq[outindex[k]])
    # plt.legend()
    # plt.semilogy()
# plt.show()

if __name__=='__main__':
    plot_HR67_hist(save=0,show=1)
