#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.timeseries import LombScargle
from astropy.stats import poisson_conf_interval
import hawkeye.pfold as pfold
from scipy import optimize as op
font1 = {'family': 'Normal',
         'weight': 'normal',
         'size': 18, }
plt.rc('legend',fontsize=14 )
def plot_longT_V(src_evt,bkg_file,epoch_info,backscale=12.,iffold=False,p_test=None,shift=None,show=False):
    if epoch_info.ndim == 1:epoch_info=np.array([epoch_info])
    t_start = epoch_info[:, 0]
    t_end = epoch_info[:, 1]
    t_mid=(t_start+t_end)/2
    obsID = epoch_info[:, 2]
    expT = epoch_info[:, 3]
    # expT=t_end-t_start
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
    for i in range(len(t_mid)):
        plt.text(t_mid[i],CR[i]*1.2,str(int(obsID[i])))
    if show:
        plt.show()
    else:plt.close()

    if iffold:
        plt.figure(2)
        turns=pfold.trans(t_mid,p_test=p_test,shift=shift)
        plt.errorbar(turns, CR, CR_ERR, fmt='o', capsize=3, elinewidth=1, ecolor='red')
        plt.errorbar(turns+1, CR, CR_ERR, fmt='o', capsize=3, elinewidth=1, ecolor='red')
        plt.show()

    return CR

def sin_temp(x,period,shift,A,B):
    return A*np.sin(2*np.pi/period*x+shift)+B

def curvefit_sin(x,y,yerr,period):
    param_bounds = ((period*0.95,0,8,8),(period*1.05, 0.1,10,10))
    popt, pcov = op.curve_fit(sin_temp, x,y,bounds=param_bounds)
    perr = np.sqrt(np.diag(pcov))
    return (popt,perr)

def plot_singleobs_lc(lc,period=None,ifsin=None,shift=0,figurepath=None,save=0,show=0,dataname=None):
    plt.figure(1,(15,6))
    plt.title(r'$T_0={0}$'.format(lc.time[0]),font1)
    x=lc.time-lc.time[0]
    y2=lc.counts
    y2_err = np.array(poisson_conf_interval(y2, interval='frequentist-confidence'))
    y2_err[0] = y2 - y2_err[0]
    y2_err[1] = y2_err[1] - y2
    plt.figure(1,(9,6))
    if ifsin:
        (popt,perr)=curvefit_sin(x,y2,0.5*(y2_err[0]+y2_err[1]),period)
        print(popt)
        # y1 = sin_temp(x,popt[0],popt[1],popt[2],popt[3])
        x1=np.linspace(x.min(),x.max(),10000)
        y1 = sin_temp(x1,period,6.28,8,8)
        plt.plot(x1,y1)
    # absolute_sigma = True, sigma = yerr,
    plt.errorbar(x, y2,yerr=y2_err, fmt='co', capsize=4, elinewidth=2, ecolor='red',color='green')
    plt.xlabel(r'Time-$T_0$ (second)',font1)
    plt.ylabel('Counts/bin',font1)
    plt.tick_params(labelsize=16)
    if save:plt.savefig(figurepath+f'{dataname}_lc_4longobs.eps',bbox_inches='tight', pad_inches=0.0)
    if show:plt.show()
    else:plt.close()
