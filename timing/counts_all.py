# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import string
from scipy import interpolate
# path='/Users/baotong/Desktop/period/txt_all_obs_m03/'
# os.chdir(path)
# ID=np.linspace(1,3619,3619)
# cts=[]
# for src in ID:
#     src=int(src)
#     event=np.loadtxt(str(src)+'.txt')
#     cts.append(len(event))
# #print(np.sort(cts))
# cts=np.array(cts)
# cand_id=ID[np.where(cts>100)]
# print(cand_id)
#np.savetxt(path+'cand_id.txt',cand_id,fmt="%5d")

font1 = {'family': 'Normal',
         'weight': 'normal',
         'size': 18, }
def plot_obs():
    bins = np.logspace(np.log10(5000), np.log10(200000), 21)
    path_NSC='/Users/baotong/Desktop/period/txt/'
    path_LW='/Users/baotong/Desktop/period_LW/txt_all_obs/'
    ep_NSC_I=np.loadtxt(path_NSC+'SgrA_I_epoch.txt')[:,-1]
    ep_NSC_G=np.loadtxt(path_NSC+'SgrA_G_epoch.txt')[:,-1]
    ep_LW=np.loadtxt(path_LW + 'LW_epoch.txt')[:,-1]
    plt.semilogx()
    plt.xlabel('ACIS exposure time (s)',font1)
    plt.ylabel('N',font1)
    plt.tick_params(labelsize = 18)
    plt.hist(ep_NSC_I,bins=bins,histtype='step',lw=2,color='red')
    plt.hist(ep_NSC_G, bins=bins,histtype='step',lw=4,linestyle=':',color='green')
    plt.hist(ep_LW, bins=bins,histtype='step',lw=5,linestyle='-.',color='blue')
    plt.legend(['NSC_ACIS-I','NSC_ACIS-S','LW_ACIS-I'])
    plt.show()
#plot_obs()

def plot_obs_bin():
    path_NSC = '/Users/baotong/Desktop/period/txt/'
    path_LW = '/Users/baotong/Desktop/period_LW/txt_all_obs/'
    NSC_I = np.loadtxt(path_NSC + 'SgrA_I_epoch.txt')
    NSC_G = np.loadtxt(path_NSC + 'SgrA_G_epoch.txt')
    LW = np.loadtxt(path_LW + 'LW_epoch.txt')

    ep_NSC_I = NSC_I[:, -1]
    ep_NSC_G = NSC_G[:, -1]
    ep_LW = LW[:, -1]
    obs_time_I = (NSC_I[:, 0] + NSC_I[:, 1]) / 2
    time_I = obs_time_I / 86400 + 2449352.5 - 2400000.5
    obs_time_G = (NSC_G[:, 0] + NSC_G[:, 1]) / 2
    time_G = obs_time_G / 86400 + 2449352.5 - 2400000.5
    obs_time_LW = (LW[:, 0] + LW[:, 1]) / 2
    time_LW = obs_time_LW / 86400 + 2449352.5 - 2400000.5
    plt.semilogy()
    plt.bar(time_I,ep_NSC_I,width=50)
    plt.bar(time_G, ep_NSC_G, width=50)
    plt.bar(time_LW, ep_LW, width=50)
    plt.xlabel('Observation time (MJD)', font1)
    plt.ylabel('ACIS exposure time (s)', font1)
    plt.legend(['NSC_ACIS-I', 'NSC_ACIS-S', 'LW_ACIS-I'])
    plt.tick_params(labelsize=18)
    plt.show()
plot_obs_bin()

def plot_GL():
    font1 = {'family': 'Normal',
             'weight': 'normal',
             'size': 18,}
    N=1000
    ban=np.sort(np.random.random(9))
    y=[]
    y.append(int(N*ban[0]))
    for i in range(8):
        y.append(int(N*ban[i+1]-N*ban[i]))
    y.append(1000-int(N*ban[8]))
    print(y)
    x=np.linspace(0.05,0.95,10)
    y=np.random.random(10)*100+50
    y=np.array([30,30,50,90,300,300,90,50,30,30])


    func = interpolate.interp1d(x, y, kind='cubic')
    x_new = np.linspace(min(x),max(x),100)
    y_new = func(x_new)
    plt.bar(x,y,width=0.1,color='white',edgecolor='black')
    plt.plot(x_new,y_new)
    plt.xlabel('phase',font1)
    plt.ylabel('counts/bin',font1)
    plt.tick_params(labelsize=18)
    plt.show()
#plot_GL()