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
# import scipy.signal.lombscargle as LombScargle
from astropy.stats import poisson_conf_interval
import stingray as sr
from stingray.events import EventList
from stingray.lightcurve import Lightcurve
from stingray import Lightcurve, Crossspectrum, sampledata,Powerspectrum,AveragedPowerspectrum
from stingray.simulator import simulator, models

figurepath = '/Users/baotong/Desktop/aas/AGN_CDFS_mod1/figure/'

font1 = {'family': 'Normal',
         'weight': 'normal',
         'size': 18, }
font2 = {'family': 'Normal',
         'weight': 'normal',
         'size': 16, }
def bending_po(x,p):
    """
    Parameters
    ----------
    :param x:numpy.ndarray
        non-zero frequencies
    :param p:
    p[0] = bending frequency,delta
    p[1] = alpha,即 power law index、
    p[2] = constant,gamma
    p[3] = normalization N or beta
    :return:
    """
    return p[3]*x**(-1)*(1+(x/p[0])**(p[1]-1))**(-1)+ p[2]

def standard_lorentzian(x,p):
    """
    standard Lorentzian function.

    Parameters
    ----------
    x: numpy.ndarray
        non-zero frequencies

    p: iterable
        p[0] = peak centeral frequency
        p[1] = Q=f/df
        p[2] = R （should represent fractional RMS)
        p[3] = 2, by default
    """
    if p[3]!=2: print("warning for the index")
    return p[2]**2*2*p[1]*p[0]/(np.pi*(p[0]**2+(2*p[1])**2*(x-p[0])**2))

def bendp_lorentz(x,p_1,p_2):

    return p_1[3]*x**(-1)*(1+(x/p_1[0])**(p_1[1]-1))**(-1)+ \
           p_1[2]+p_2[2]**2*2*p_2[1]*p_2[0]/(np.pi*(p_2[0]**2+(2*p_2[1])**2*(x-p_2[0])**2))

def generalized_lorentzian(x, p):
    """
    Generalized Lorentzian function.

    Parameters
    ----------
    x: numpy.ndarray
        non-zero frequencies

    p: iterable
        p[0] = peak centeral frequency
        p[1] = FWHM of the peak (gamma)
        p[2] = peak value at x=x0
        p[3] = power coefficient [n]
    Returns
    -------
    model: numpy.ndarray
        generalized lorentzian psd model
    """
    assert p[3] > 0., "The power coefficient should be greater than zero."
    return p[2] * (p[1] / 2)**p[3] * 1./(abs(x - p[0])**p[3] + (p[1] / 2)**p[3])
def smoothbknpo(x, p):
    """
    Smooth broken power law function.
    Parameters
    ----------

    x: numpy.ndarray
        non-zero frequencies

    p: iterable
        p[0] = normalization frequency
        p[1] = power law index for f --> zero
        p[2] = power law index for f --> infinity
        p[3] = break frequency
    Returns
    -------
    model: numpy.ndarray
        generalized smooth broken power law psd model
    """
    return p[0] * x**(-p[1]) / (1. + (x / p[3])**2)**(-(p[1] - p[2]) / 2)

def sim_evtlist(lc):
    num_evt_bin=np.random.poisson(lc.counts)
    evt_all=[]
    i=0
    while i < len(num_evt_bin):
        evt_all=np.concatenate((evt_all,(np.random.uniform(lc.time[i]-lc.dt/2,lc.time[i]+lc.dt/2,num_evt_bin[i]))))
        i+=1
    return evt_all

def filter_energy(time,energy,band):
    T=time
    E=energy
    i=0
    if len(time)!=len(energy):
        print('error')
        return None
    else:
        while i <len(E):
            if E[i]<band[0] or E[i]>band[1]:
                E=np.delete(E,i)
                T=np.delete(T,i)
            else:
                i+=1
    return (T,E)

def get_hist(t, len_bin,tstart=0,tstop=0):
    ###将输入的time信息，按照len_bin的长度输出为lc
    if tstart==0 and tstop==0:
        tstart=t[0]
        tstop=t[-1]
        tseg=tstop-tstart
    else:tseg=tstop-tstart

    t_test = t;dt=len_bin
    ev = EventList()
    ev.time = t_test
    lc_new = ev.to_lc(dt=dt, tstart=tstart, tseg=tseg)
    return lc_new

def get_hist_withbkg(t,t_bkg, len_bin,tstart=0,tstop=0):
    ###将输入的time信息，按照len_bin的长度输出为lc
    if tstart==0 and tstop==0:
        tstart=t[0]
        tstop=t[-1]
        tseg=tstop-tstart
    else:
        tseg = tstop - tstart
    t_test = t;t_bkg_test=t_bkg;dt=len_bin
    if type(t_bkg_test)==list:t_bkg_test=np.array(t_bkg_test)
    t_bkg_test=np.delete(t_bkg_test,t_bkg_test<0)
    ev = EventList();ev_bkg=EventList()
    ev.time=t_test;ev_bkg.time=t_bkg_test
    lc_new = ev.to_lc(dt=dt, tstart=tstart, tseg=tseg)
    lc_bkg = ev_bkg.to_lc(dt=dt, tstart=tstart, tseg=tseg)
    lc_out=lc_new
    lc_out.counts=lc_new.counts-(1/12.)*lc_bkg.counts
    # lc_out.counts = lc_new.counts
    # print(lc_bkg.counts)
    return lc_out

def get_LS(time, flux,freq):
    x = time
    y = flux
    LS = LombScargle(x, y,normalization = 'standard')
    # LS = LombScargle(x, y, normalization='psd')
    power = LS.power(freq)
    max_NormLSP=np.max(power)
    FP=LS.false_alarm_probability(power.max(),minimum_frequency = freq[0],maximum_frequency = freq[-1],method='baluev')
    FP_99 = LS.false_alarm_level(0.0027,minimum_frequency = freq[0], maximum_frequency = freq[-1],method='baluev')
    FP_95 = LS.false_alarm_level(0.05, minimum_frequency=freq[0],
                                 maximum_frequency=freq[-1], method='baluev')
    FP_68 = LS.false_alarm_level(0.32,minimum_frequency=freq[0],
                                 maximum_frequency=freq[-1], method='baluev')

    # if FP<0.01:print(dataname)
    # plt.title('Epoch {2}: XID={0},FAP={1}'.format(dataname,np.round(FP,4),k),font1)
    # plt.semilogx()
    plt.plot(freq, power)
    plt.semilogx()
    out_period=1./freq[np.where(power==np.max(power))][0]
    plt.plot([freq[0], freq[-1]], [FP_99, FP_99], '--')
    plt.plot([freq[0], freq[-1]], [FP_95, FP_95], '--')
    plt.plot([freq[0], freq[-1]], [FP_68, FP_68], '--')
    plt.text(freq[0], FP_99, '1-FAP 99%',font1)
    plt.text(freq[0], FP_95, '95%',font1)
    plt.text(freq[0], FP_68, '68%',font1)
    plt.xlabel('Frequency (Hz)',font1)
    plt.ylabel('Normalized LS Periodogram',font1)
    plt.tick_params(labelsize=16)
    # plt.show()
    # plt.savefig('/Users/baotong/Desktop/CDFS/fig_LS_ep{0}_ovsamp_5_baluev/{1}.eps'.format(k,dataname))
    plt.close()
    return [FP,out_period,max_NormLSP]

def get_T_in_mbins(epoch_file,w,m,fi):
    T=2*np.pi/w
    T_in_perbin = np.zeros(m)
    # 每个bin的总积分时间
    tbin = T/m
    # 每个bin的时间长度
    epoch_info = np.loadtxt(epoch_file)
    if epoch_info.ndim==2:
        t_start = epoch_info[:, 0]
        t_end = epoch_info[:, 1]
    if epoch_info.ndim == 1:
        t_start=np.array([epoch_info[0]])
        t_end = np.array([epoch_info[1]])

    N_bin_t_start=t_start/tbin+m*fi/(2*np.pi)
    N_bin_t_end=t_end/tbin+m*fi/(2*np.pi)
    intN_bin_t_start=np.floor(N_bin_t_start)+1
    intN_bin_t_end=np.floor(N_bin_t_end)
    intN_bin_t_start=intN_bin_t_start.astype(int)
    intN_bin_t_end=intN_bin_t_end.astype(int)
    for i in range(len(N_bin_t_start)):
        if intN_bin_t_end[i]>=intN_bin_t_start[i]:
            T_in_perbin+=int((intN_bin_t_end[i]-intN_bin_t_start[i])/m)*tbin
            #print(intN_bin_t_start[i]-1)
            T_in_perbin[np.mod(intN_bin_t_start[i],m)-1]+=(intN_bin_t_start[i]-N_bin_t_start[i])*tbin
            T_in_perbin[np.mod(intN_bin_t_end[i],m)]+=(N_bin_t_end[i]-intN_bin_t_end[i])*tbin
            rest=np.mod(intN_bin_t_end[i]-intN_bin_t_start[i],m)
            for k in range(rest):
                T_in_perbin[int(np.mod((intN_bin_t_start[i] + k), m))] += tbin
            #print(rest)
        else:
            T_in_perbin[np.mod(intN_bin_t_start[i],m)-1]+=(N_bin_t_end[i]-N_bin_t_start[i])*tbin
    return T_in_perbin

def filter_obs(src_evt,bkg_evt,useid):
    src_evt_use = src_evt[np.where(src_evt[:-1] == useid[0])[0]]
    bkg_evt_use = bkg_evt[np.where(bkg_evt[:-1] == useid[0])[0]]
    i=1
    while i < len(useid):
        id=useid[i]
        src_evt_use_temp=src_evt[np.where(src_evt[:-1]==id)[0]]
        bkg_evt_use_temp=bkg_evt[np.where(bkg_evt[:-1]==id)[0]]
        src_evt_use = np.concatenate((src_evt_use, src_evt_use_temp))
        bkg_evt_use = np.concatenate((bkg_evt_use, bkg_evt_use_temp))
        i+=1
    return (src_evt_use,bkg_evt_use)

def read_epoch(epoch_file):
    epoch=np.loadtxt(epoch_file)
    if epoch.ndim==1:
        epoch=np.array([epoch])
    # print(epoch)
    TSTART=epoch[:,0]
    TSTOP=epoch[:,1]
    OBSID=epoch[:,2]
    exptime=epoch[:,3]
    return (TSTART,TSTOP,OBSID,exptime)