#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import re
import string
from datetime import datetime
from scipy.interpolate import lagrange
from scipy import interpolate
from scipy.optimize import curve_fit
from scipy.special import comb, perm
import pandas as pd
from astropy.stats import poisson_conf_interval
from astropy.timeseries import LombScargle
import scipy
def read_LSP(k,srcid,threshold=0.999,num_trials=10000):
    # smooth=100
    path='/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep{0}/simulation/'.format(k)
    T_exp = 11000154.981141508;dt=100
    freq = np.arange(1 / T_exp, 0.5 / dt, 1 / (5 * T_exp))
    freq = freq[np.where(freq > 1 / 20000.)]
    res=[];
    df=pd.read_table(path+'{0}_LS_simP.csv'.format(srcid),header=None)
    df=np.array(df)
    print('read success')
    for i in range(len(df)):
        temp=np.array(df[i][0].split(','))
        temp=temp.astype('float32')
        res.append(temp)
    res=np.array(res)
    res=res.T
    res_sort=np.sort(res)
    # print('local num={0}'.format(len(np.where(res[55100]>12.431478151451798)[0])))
    # print('global num={0}'.format(len(sum(np.where(res>12.431478151451798)))))
    maxP_trial=np.max(res,axis=0)
    np.savetxt(path + '{0}_LS_sim_noQPO/maxP_trial.txt'.format(srcid), maxP_trial, fmt='%10.2f')

    power_thres_ar = np.zeros(len(freq));
    power_thres_ar=res_sort[:,int(num_trials*threshold)]
    smoothP=np.column_stack((freq,power_thres_ar))
    np.savetxt(path+'{0}_LS_sim_noQPO/LSP_{1}.txt'.format(srcid,threshold),smoothP,fmt='%10.5f %10.2f')
    plt.step(freq,power_thres_ar)
    plt.savefig(path+'{0}_LS_sim_noQPO/LSP_{1}.eps'.format(srcid,threshold))
    return None
# read_LSP('3','19',threshold=0.999,num_trials=10000)
read_LSP('3','236',threshold=0.9999,num_trials=10000)
# read_LSP('3','19',threshold=0.9998,num_trials=10000)