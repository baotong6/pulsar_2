# -*- coding: utf-8 -*-
"""
Created on Sun April 6 13:56:40 2023
@author: baotong

time: photon arrvial time
P: assumed period
m: number of bin in one period

"""

import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp
import functools
import datetime
import hawkeye as hawk
from NGC104.timing_comb import load_data

def f_window(t,x,r):
    f=(np.sign(r-np.abs(x-t))+1)/2
    return f


def photon_density_func(tlist,epoch_info=[],width=10000):
    tlist=tlist-tlist[0]
    x=np.arange(tlist[0]+width/2,tlist[-1]-width/2,width/200)
    y=np.zeros(len(x))
    for k in range(len(x)):
        y[k]=np.sum(f_window(tlist,x[k],width))
    plt.plot(x,y)
    plt.show()


def extract_var_list(time,epoch_info,P,m):
    phase=P*time
    T_tot=time[-1]-time[0]
    variance=[]

if __name__=='__main__':
    path_Tuc = '/Users/baotong/Desktop/period_M28/txt_all_obs_p90/'
    (src_evt_use, epoch_info_use)=load_data(dataname=106,ifpath=path_Tuc,ecf=90,ifobsID=[9133])
    time=src_evt_use[:,0]
    photon_density_func(time,width=2000)