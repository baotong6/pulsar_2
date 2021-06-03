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
import useful_functions as func
from stingray import Lightcurve, Powerspectrum, AveragedPowerspectrum

font1 = {'family': 'Normal',
         'weight': 'normal',
         'size': 18, }
def get_lc(data_file,bkg_file,epoch_file,bin,id1,id2):
    src_evt=np.loadtxt(data_file)
    bkg_evt=np.loadtxt(bkg_file)
    epoch=np.loadtxt(epoch_file)
    use_id = epoch[:, 2][id1:id2]
    TSTART=epoch[:,0][id1:id2]
    TSTOP=epoch[:,1][id1:id2]
    exptime=TSTOP[-1]-TSTART[0]
    (src_evt,bkg_evt)=func.filter_obs(src_evt,bkg_evt,use_id)
    time=src_evt[:,0];energy=src_evt[:,1]
    time_bkg= bkg_evt[:, 0];energy_bkg= bkg_evt[:, 1]

    lc_src = func.get_hist(time, len_bin=bin,tstart=TSTART[0],tstop=TSTOP[-1])
    lc_net=func.get_hist_withbkg(time, time_bkg,len_bin=bin,tstart=TSTART[0],tstop=TSTOP[-1])

    return [lc_src,lc_net,TSTART,TSTOP]
def plot_lc(data_file,bkg_file,epoch_file,p_test,bin,shift,label,id1,id2):
    [lc_src,lc_net,TSTART,TSTOP]=get_lc(data_file,bkg_file,epoch_file,bin,id1,id2)
    plt.errorbar(lc_src.time,lc_src.counts,yerr=np.sqrt(lc_src.counts))
    # plt.errorbar(lc_net.time, lc_net.counts, yerr=np.sqrt(lc_src.counts))

    exptime=TSTOP[-1]-TSTART[0]
    for i in range(int(exptime/p_test)):
        plt.plot([TSTART[0]+i*p_test+shift*p_test,TSTART[0]+i*p_test+shift*p_test],[0,10],'--',color='grey')
    plt.show()

def plot_pds(data_file,bkg_file,epoch_file,p_test,bin,shift,label,id1,id2):
    [lc_src, lc_net, TSTART, TSTOP] = get_lc(data_file, bkg_file, epoch_file, bin, id1, id2)
    ps = Powerspectrum(lc_src,norm='frac')
    fig, ax1 = plt.subplots(1, 1, figsize=(9, 6), sharex=True)
    ax1.loglog()
    ax1.step(ps.freq, ps.power, lw=2, color='blue')
    ax1.set_xlabel("Frequency (Hz)", fontproperties=font1)
    ax1.set_ylabel("Power ", fontproperties=font1)
    ax1.set_yscale('log')
    ax1.tick_params(axis='x', labelsize=16)
    ax1.tick_params(axis='y', labelsize=16)
    ax1.tick_params(which='major', width=1.5, length=7)
    ax1.tick_params(which='minor', width=1.5, length=4)
    for axis in ['top', 'bottom', 'left', 'right']:
        ax1.spines[axis].set_linewidth(1.5)
    plt.show()

    # avg_ps = AveragedPowerspectrum(lc_src, 20000,dt=lc_src.time[1]-lc_src.time[0],norm='leahy')
    # print("Number of segments: %d" % avg_ps.m)
    # fig, ax1 = plt.subplots(1, 1, figsize=(9, 6))
    # ax1.loglog()
    # ax1.step(avg_ps.freq, avg_ps.power, lw=2, color='blue')
    # ax1.set_xlabel("Frequency (Hz)", fontproperties=font1)
    # ax1.set_ylabel("Power ", fontproperties=font1)
    # ax1.set_yscale('log')
    # ax1.tick_params(axis='x', labelsize=16)
    # ax1.tick_params(axis='y', labelsize=16)
    # ax1.tick_params(which='major', width=1.5, length=7)
    # ax1.tick_params(which='minor', width=1.5, length=4)
    # for axis in ['top', 'bottom', 'left', 'right']:
    #     ax1.spines[axis].set_linewidth(1.5)
    # plt.show()



if __name__=='__main__':
    path='/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep4/'
    id='643'
    data_file=path+'{0}.txt'.format(id)
    bkg_file=path+'{0}_bkg.txt'.format(id)
    epoch_file=path+'epoch_src_{0}.txt'.format(id)
    p_test= 13273.02779663
    # plot_lc(data_file, bkg_file, epoch_file, p_test, bin=2000, shift=0.6,label=id,id1=6,id2=7)
    plot_pds(data_file, bkg_file, epoch_file, p_test, bin=30, shift=0.6,label=id,id1=27,id2=28)