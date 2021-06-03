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
from astropy.timeseries import LombScargle
from functools import reduce
# import scipy.signal.lombscargle as LombScargle
from astropy.stats import poisson_conf_interval
threshold=1-0.
def good_function_hhh(k):
    path='/Users/baotong/Desktop/CDFS/fig_LS_ep{0}_ovsamp_5_baluev/'.format(k)
    LS_res=np.loadtxt(path+'LS_result_{0}.txt'.format(k))
    cr = LS_res[:, 5]
    # LS_res=LS_res[np.where(cr>0)]
    src_id=LS_res[:,0];conf=LS_res[:,1];period=LS_res[:,2];cr = LS_res[:, 5]
    good_det_id=src_id[reduce(np.intersect1d, (np.where(conf < threshold), np.where(period > 360.),
                                               np.where(np.abs(period -707.)>5),np.where(np.abs(period -999.)>5),
                                               np.where(cr>4e-4),np.where(period<20000)))]
    # good_det_id = src_id[reduce(np.intersect1d, (np.where(conf < threshold), np.where(np.abs(period - 707.) > 5),np.where(np.abs(period - 999.) > 5)))]
    good_det_id=good_det_id.astype('int')
    return good_det_id
ID1=good_function_hhh(1)
ID2=good_function_hhh(2)
ID3=good_function_hhh(3)
ID4=good_function_hhh(4)
print(len(ID1))
print(len(ID2))
print(ID3)
print(len(ID4))

def sim_res_pop(k):
    threshold = 0.9972
    cts_rate=['2e-5','5e-5','1e-4','2e-4','3e-4','4e-4','5e-4']
    det_p=[]
    path = '/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep{0}/simulation/'.format(k)
    for i in range(len(cts_rate)):
        filename='trial_out_{0}_REJ1034+396_noQPO.txt'.format(cts_rate[i])
        sim_res=np.loadtxt(path+filename)
        conf=sim_res[:,0];period=sim_res[:,1]
        det_p=np.concatenate((det_p,period[np.where(conf<1-threshold)]))

    frac_400=len(np.where(det_p<500)[0])/len(det_p)
    return frac_400
# print(sim_res_pop(1))
