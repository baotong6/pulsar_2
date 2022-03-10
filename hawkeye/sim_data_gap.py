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
import hawkeye as hawk

def load_data(dataname,ecf=90,obsID=None):
    path_Tuc='/Users/baotong/Desktop/period_Tuc/txt_all_obs_p{0}/'.format(ecf)
    path = path_Tuc
    dataname = '{0}.txt'.format(dataname)
    epoch_file = path + 'epoch_src_' + dataname
    src_evt=np.loadtxt(path+dataname)
    epoch_info=np.loadtxt(epoch_file)
    CR=hawk.plot_longT_V(src_evt=src_evt, bkg_file=None,epoch_info=epoch_info)
    # print(epoch_info[:,2])
    CR/=ecf/100.
    (useid, epoch_info_use)=hawk.choose_obs(epoch_info,flux_info=CR,
                                            flux_filter=9,expT_filter=10000,
                                            if_flux_high=False, if_expT_high=True,obsID=obsID)
    src_evt_use =hawk.filter_obs(src_evt, useid)
    print(useid)
    return (src_evt_use,epoch_info_use)

def add_data_gap(evt1,evt2):
    time1=evt1[:,0];time2=evt2[:,0]
    gap=time2[0]-time1[-1]
    print(gap)
    CR1=len(time1)/(time1[-1]-time1[0])
    CR2=len(time2)/(time2[-1]-time2[0])
    counts_gap=gap*(CR1+CR2)/2
    fake_time=np.random.random(int(counts_gap))*gap+time1[-1]
    fake_time=np.sort(fake_time)
    fake_energy=np.random.random(int(counts_gap))*8000+500
    fake_energy=fake_energy.astype('int')
    fake_evt=np.column_stack((fake_time,fake_energy))
    return fake_evt

def main_process():
    (src_evt_use1, epoch_info_use1) = load_data(dataname='999', ecf=90,obsID=[953])
    (src_evt_use2, epoch_info_use2) = load_data(dataname='999', ecf=90,obsID=[955])
    fake_evt=add_data_gap(src_evt_use1,src_evt_use2)
    obsID=np.zeros(len(fake_evt))+99999
    fake_evt=np.column_stack((fake_evt,obsID))
    print(len(fake_evt))
    np.savetxt('/Users/baotong/Desktop/period_Tuc/txt_all_obs_p{0}/fake_gap_999.txt'.format(90),fake_evt,fmt='%.7f  %5.3f  %d')


if __name__=='__main__':
    main_process()