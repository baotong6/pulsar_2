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

def load_data(dataname,ecf=90):
    # path_Tuc='/Users/baotong/Desktop/period_Tuc/txt_startover/txt_all_obs_p{0}/'.format(ecf)
    # path_Tuc = '/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep4/'
    path_Tuc = '/Users/baotong/Desktop/period_Tuc/txt_startover/txt_all_obs_p{0}/'.format(ecf)
    path_M31='/Users/baotong/Desktop/M31XRB/M31ACIS_txt/txt_all_obs_p90/'
    path = path_Tuc
    dataname = '{0}.txt'.format(dataname)
    epoch_file = path + 'epoch_src_' + dataname
    src_evt=np.loadtxt(path+dataname)
    epoch_info=np.loadtxt(epoch_file)
    if epoch_info.ndim==1:epoch_info=np.array([epoch_info])
    CR=hawk.plot_longT_V(src_evt=src_evt, bkg_file=None,epoch_info=epoch_info)
    CR/=ecf/100.
    useobsID=epoch_info[:,2][0:50].astype('int')
    print(useobsID)
    (useid, epoch_info_use)=hawk.choose_obs(epoch_info,flux_info=CR,
                                            flux_filter=10000,expT_filter=4000,
                                            if_flux_high=0, if_expT_high=1,obsID=[2735,2736,2737,2738])
    src_evt_use =hawk.filter_obs(src_evt, useid)
    print(useid)
    return (src_evt_use,epoch_info_use)

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
    (src_evt_use, epoch_info_use)=load_data(dataname=366,ecf=90)
    time=src_evt_use[:,0]
    photon_density_func(time,width=2000)