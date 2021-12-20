#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.timeseries import LombScargle
import hawkeye.pfold as pfold

def plot_longT_V(src_evt,bkg_file,epoch_info,backscale=12.,iffold=False,p_test=None,shift=None,show=True):
    if epoch_info.ndim == 1:epoch_info=np.array([epoch_info])
    t_start = epoch_info[:, 0]
    t_end = epoch_info[:, 1]
    t_mid=(t_start+t_end)/2
    obsID = epoch_info[:, 2]
    expT = epoch_info[:, 3]
    cts=[];bkg_cts=[]
    if not bkg_file:
        for i in range(len(obsID)):
            cts.append(len(np.where(src_evt[:, 2] == obsID[i])[0]))
        cts = np.array(cts)
        CR = cts / expT
        CR_ERR = np.sqrt(CR * expT) / expT

    else:
        time_bkg = np.loadtxt(bkg_file)
        for i in range(len(obsID)):
            cts.append(len(np.where(src_evt[:,2]==obsID[i])[0]))
            bkg_cts.append(len(np.where(time_bkg[:, 2] == obsID[i])[0]))
        cts=np.array(cts);bkg_cts=np.array(bkg_cts)
        CR=(cts-bkg_cts/backscale)/expT
        CR_ERR=np.sqrt(CR*expT)/expT
    plt.figure(1)
    plt.semilogy()
    plt.errorbar(t_mid,CR,CR_ERR,fmt='o',capsize=3, elinewidth=1, ecolor='red')
    if show:
        plt.show()

    if iffold:
        plt.figure(2)
        turns=pfold.trans(t_mid,p_test=p_test,shift=shift)
        plt.errorbar(turns, CR, CR_ERR, fmt='o', capsize=3, elinewidth=1, ecolor='red')
        plt.errorbar(turns+1, CR, CR_ERR, fmt='o', capsize=3, elinewidth=1, ecolor='red')
        plt.show()

    return CR

def plot_singleobs_lc(lc):
    plt.title('T0={0}'.format(lc.time[0]))
    plt.step(lc.time-lc.time[0],lc.counts)
    plt.show()
