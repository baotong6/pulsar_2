#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import string
from datetime import datetime
import random
from scipy import optimize
from scipy.optimize import fsolve
from itertools import chain
from astropy.timeseries import LombScargle
def get_time_series(T_START,N_cts,T_stop,period=300.,amp=0.7, model='sin',eqw=0.1):
    #N_cts为总光子数
    #T为总观测时间
    #period为周期
    #amp为model的amplitude,反应了时变幅度
    #model现只针对sin和eclipse函数
    def get_turns(t,period):
        v=1.0/period
        turns=t*v-int(t*v)
        return turns

    if model=='eclipse':
        i=0
        delta=amp
        w=2*np.pi/period
        T=T_stop-T_START
        lam_0=N_cts/T
        t=np.random.random(N_cts)*T+T_START
        t=np.sort(t)
        t_eclipse=t
        while i<len(t_eclipse):
            if (0.5-0.5*eqw)<get_turns(t_eclipse[i],period)<(0.5+0.5*eqw):
                rand=np.random.random(1)[0]
                if rand<amp:
                    t_eclipse=np.delete(t_eclipse,i)
                else:
                    i+=1
            else:i+=1
        return t_eclipse

    if model=='sin':
        delta=amp
        w=2*np.pi/period
        T=T_stop-T_START
        lam_0=N_cts/T
        t=np.random.random(N_cts)*T+T_START
        t=np.sort(t)
        t_sin=t
        i=0
        while i <len(t):
            def f(t1):
                return t1 - delta / w * np.cos(w * t1) - t[i]
            t_sin[i]=fsolve(f,t[i]+np.random.rand(1)[0])
            i+=1
        return t_sin

    if model=='const':
        T=T_stop-T_START
        lam_0=N_cts/T
        t=np.random.random(N_cts)*T+T_START
        t=np.sort(t)
        return t

def get_epoch_time_series(cts_rate,period,amp, model):
    t=[]
    path='/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep4/'
    epoch = np.loadtxt(path+'epoch_src_384.txt')
    # epoch = np.loadtxt('LW_epoch.txt')
    tstart=epoch[:,0]
    tstop=epoch[:,1]
    exp_time_epoch=epoch[:,-1]
    if model=='const':
        for i in range(len(tstart)):
            #VI=10
            cts = int((cts_rate * (tstop[i] - tstart[i]))*(np.random.random(1)[0]*3.16227766+0.316227766))
            t.append(get_time_series(tstart[i], cts, tstop[i], period, amp, model=model))
    else:
        for i in range(len(tstart)):
            cts = int(cts_rate * (tstop[i] - tstart[i]))
            t.append(get_time_series(tstart[i],cts,tstop[i],period,amp,model=model))
    t=list(chain.from_iterable(t))
    t = np.array(t)
    # t = t + np.random.rand(len(t)) * 3.2 - 1.6
    return t

