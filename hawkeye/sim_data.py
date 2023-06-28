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
        t=np.random.random(np.random.poisson(N_cts))*T+T_START
        t=np.sort(t)
        return t

def get_epoch_time_series(cts_rate,period,amp, model,epoch_info):
    t=[]
    if epoch_info.ndim==1:epoch_info=np.array([epoch_info])
    # epoch = np.loadtxt('LW_epoch.txt')
    tstart=epoch_info[:,0]
    tstop=epoch_info[:,1]
    exp_time_epoch=epoch_info[:,-1]
    obsid=epoch_info[:,2]
    obsidout=[]
    if model=='const':
        for i in range(len(tstart)):
            if type(cts_rate) == list or type(cts_rate) == np.ndarray:
            #VI=10
                cts = int((cts_rate[i] * (tstop[i] - tstart[i])))
            else:
                cts = int((cts_rate * (tstop[i] - tstart[i])))
            tempt=get_time_series(tstart[i], cts, tstop[i], period, amp, model=model)
            t.append(tempt)
            obsidout.append([int(obsid[i]) for k in range(len(tempt))])
    else:
        for i in range(len(tstart)):
            cts = int(cts_rate * (tstop[i] - tstart[i]))
            t.append(get_time_series(tstart[i],cts,tstop[i],period,amp,model=model))
    t=list(chain.from_iterable(t));obsidout=list(chain.from_iterable(obsidout))
    t = np.array(t);obsidout=np.array(obsidout)
    t = t + np.random.rand(len(t)) * 3.2 - 1.6
    return t,obsidout

