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

def get_info(srcid):
    src_cts_ep=[];bkg_cts_ep=[];exp_ep=[];CR_ep=[];net_cts_ep=[]
    for k in range(1,6):
        if k<5:
            path='/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep{0}/'.format(k)
        if k==5:
            path = '/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8/'

        src_evt=np.loadtxt(path+'{0}.txt'.format(srcid))
        bkg_evt=np.loadtxt(path+'{0}_bkg.txt'.format(srcid))

        epoch_file=np.loadtxt(path+'epoch_src_{0}.txt'.format(srcid))
        if len(epoch_file) == 0:
            src_cts_ep.append(0);
            bkg_cts_ep.append(0)
            net_cts_ep.append(0)
            exp_ep.append(0)
            CR_ep.append(0)
            continue
        if epoch_file.ndim == 1:
            epoch_file = np.array([epoch_file])
        src_cts_ep.append(len(src_evt));bkg_cts_ep.append(len(bkg_evt))
        net_cts_ep.append(len(src_evt)-len(bkg_evt)/12.)
        exp_ep.append(np.sum(epoch_file[:,3]))
        CR_ep.append((len(src_evt)-len(bkg_evt)/12.)/np.sum(epoch_file[:,3]))

    result=np.concatenate(([srcid],src_cts_ep,bkg_cts_ep,net_cts_ep,exp_ep,CR_ep))

    return result
def write_to_txt():
    path = '/Users/baotong/Desktop/CDFS/'
    source_id=np.arange(1,1055,1)
    all_result=get_info(source_id[0])
    print(get_info(source_id[2]))
    for i in range(1,len(source_id)):
        all_result=np.row_stack((all_result,get_info(source_id[i])))
    np.savetxt(path+'source_info.txt',all_result,fmt='%6d %8d %8d %8d %8d %8d %8d %8d %8d %8d %8d %10.2f %10.2f %10.2f %10.2f %10.2f'
                                                     '%12.2f %12.2f %12.2f %12.2f %12.2f %12f %12f %12f %12f %12f')

if __name__=='__main__':
    write_to_txt()
