# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp
import functools
import datetime
from astropy.io import fits
import sys
import os
import string
#import correct as correct
import random
from scipy import optimize
from sympy import *
from scipy.optimize import fsolve
from itertools import chain
def get_time_series(T_START,N_cts,T_stop,period,amp, model,eqw=0.1):
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
    epoch = np.loadtxt('CDFS_epoch_ep4.txt')
    # epoch = np.loadtxt('LW_epoch.txt')
    tstart=epoch[:,0]
    tstop=epoch[:,1]
    exp_time_epoch=epoch[:,-1]
    if model=='const':
        for i in range(len(tstart)):
            cts = int((cts_rate * (tstop[i] - tstart[i])))
            t.append(get_time_series(tstart[i], cts, tstop[i], period, amp, model=model))
    else:
        for i in range(len(tstart)):
            cts = int(cts_rate * (tstop[i] - tstart[i]))
            t.append(get_time_series(tstart[i],cts,tstop[i],period,amp,model=model))
    t=list(chain.from_iterable(t))
    t = np.array(t)
    t = t + np.random.rand(len(t)) * 3.2 - 1.6
    return t

def get_sim_time_series(dataID,mode):
    #CRT是cts-rate;CTS是光子数
    time_orig=np.loadtxt(path+str(dataID)+'.txt')
    epoch=np.loadtxt(path+'epoch_src_'+str(dataID)+'.txt')
    tstart=epoch[:,0]
    tstop=epoch[:,1]
    exp_time_epoch=epoch[:,-1]
    t=[];S=[]
    for i in range(len(epoch[:,2])):
        index=np.where(time_orig==epoch[:,2][i])
        if mode=='CRT':
            cts_num=exp_time_epoch[i]*len(time_orig)/np.sum(exp_time_epoch)
        elif mode=='CTS':
            cts_num = len(time_orig[index])
        CTS_pois = np.random.poisson(lam=int(cts_num), size=1)[0]
        t.append(get_time_series(tstart[i], CTS_pois, tstop[i], model='const'))
        #S.append(CTS_pois/exp_time_epoch[i])
    # S=np.array(S)
    # S_del = np.delete(S, np.where(S <1e-15))
    # VI=np.max(S_del)/np.min(S_del)
    t=list(chain.from_iterable(t))
    t = np.array(t)
    t = t + np.random.rand(len(t)) * 3.2 - 1.6
    return t
def get_LS(time, flux,freq,dataname,k):
    x = time
    y = flux
    # dy=np.sqrt(y)
    # plt.scatter(x,y)
    # plt.show()

    # LS = LombScargle(x, y, dy = 1, normalization = 'standard', fit_mean = True,
    #                  center_data = True).power(freq, method = 'cython')
    LS = LombScargle(x, y,normalization = 'standard')
    # LS = LombScargle(x, y, dy, normalization='psd')
    power = LS.power(freq)
    FP=LS.false_alarm_probability(power.max(),minimum_frequency = freq[0],maximum_frequency = freq[-1],method='baluev')
    FP_99 = LS.false_alarm_level(0.0027,minimum_frequency = freq[0], maximum_frequency = freq[-1],method='baluev')
    FP_90 = LS.false_alarm_level(0.05, minimum_frequency=freq[0],
                                 maximum_frequency=freq[-1], method='baluev')
    FP_68 = LS.false_alarm_level(0.32,minimum_frequency=freq[0],
                                 maximum_frequency=freq[-1], method='baluev')

    if FP<0.01:print(dataname)
    plt.title('{0},FP={1}'.format(dataname,FP))
    # FP_99 = LombScargle(x, y, dy=1, normalization='standard',
    #                     fit_mean=True, center_data=True).false_alarm_level(
    #     0.01,
    #     minimum_frequency=freq[0], maximum_frequency=freq[-1])
    #
    # FP_90 = LombScargle(x, y, dy = 1, normalization = 'standard',
    #                     fit_mean = True, center_data = True).false_alarm_level(
    #     0.1,
    #     minimum_frequency = freq[0], maximum_frequency = freq[-1])
    # FP_68 = LombScargle(x, y, dy = 1, normalization = 'standard',
    #                     fit_mean = True, center_data = True).false_alarm_level(
    #     0.32,
    #     minimum_frequency = freq[0], maximum_frequency = freq[-1])

    # plt.semilogx()
    plt.plot(freq, power)
    # print(1. / freq[np.where(power == np.max(power))])
    if FP<0.01:print(1./freq[np.where(power==np.max(power))])
    out_period=1./freq[np.where(power==np.max(power))]
    plt.plot([freq[0], freq[-1]], [FP_99, FP_99], '--')
    plt.plot([freq[0], freq[-1]], [FP_90, FP_90], '--')
    plt.plot([freq[0], freq[-1]], [FP_68, FP_68], '--')
    # plt.show()
    # plt.savefig('/Users/baotong/Desktop/CDFS/fig_LS_ep{0}_ovsamp_5_baluev/{1}.eps'.format(k,dataname))
    plt.close()
    return [FP,out_period]

def plot_sim_LS_res():
    get_epoch_time_series(cts_rate=0.0001,)