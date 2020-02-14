# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import string
#import correct as correct
from datetime import datetime
import random
from scipy import optimize
from sympy import *
from scipy.optimize import fsolve
from itertools import chain

def get_time_series(T_START,lam_0,T_stop,period,amp, model):
    #N_cts为总光子数
    #T为总观测时间
    #period为周期
    #amp为model的amplitude,反应了时变强度
    #model现只针对sin函数
    def get_turns(t,period):
        v=1.0/period
        turns=t*v-int(t*v)

        return turns

    if model=='eclipse':
        t0 = np.random.random(1)[0] * (1. / lam_0) + T_START
        t = [t0]
        delta = amp
        while t[-1] < T_stop:
            if 0.45<get_turns(t[-1],period)<0.55:
                lam = lam_0 * (1 -delta)
                temp = np.float(-1 / lam * np.log(np.random.random(1)[0]))
                t.append(t[-1] + temp)

            else:
                lam = lam_0
                temp = np.float(-1 / lam * np.log(np.random.random(1)[0]))
                t.append(t[-1] + temp)
        t = np.array(t)
        return t

    if model=='sin':
        #lam_0=N_cts/T
        t0=np.random.random(1)[0]*(1./lam_0)+T_START
        t=[t0]
        # dt=-1/lam*np.log(np.random.random(N_cts))
        while t[-1]<T_stop:
            delta = amp
            w = 2*np.pi*1.0 / period
            def f(t1):
                return t1 - delta / w * np.cos(w * t1) - t[-1]
            t_n=fsolve(f,t[-1])
            lam=lam_0*(1+delta*np.sin(w*t_n))
            #lam_simp=lam_0*(1+delta*np.sin(w*t[-1]))
            temp=np.float(-1/lam*np.log(np.random.random(1)[0]))
            #temp = -1 / lam * np.log(np.random.random(1)[0])
            t.append(t[-1]+temp)

        t=np.array(t)
    return t

def get_a_from_b(T_all,T_start,T_stop):
    #T_all为总的time series
    T_abs_1=np.abs(T_all-T_start)
    T_abs_2=np.abs(T_all-T_stop)
    index_1=np.where(T_abs_1==np.min(T_abs_1))[0][0]
    print(index_1)
    if T_all[index_1]<T_start:
        index_1+=1
    index_2=np.where(T_abs_2==np.min(T_abs_2))[0][0]
    if T_all[index_2]>T_stop:
        index_2-=1
    if index_1>=index_2:
        return []
    else:
        return T_all[index_1:index_2+1]

def get_epoch_ts_md(epoch_file,cts_rate,period,amp,model):
    t=[]
    #a = get_time_series(54268974.41, 2644, 428094744.28, period = 5540., amp = 0.5, model = 'sin')
    epoch = np.loadtxt(epoch_file)
    tstart = epoch[:, 0]
    tstop = epoch[:, 1]
    exp_time_epoch = epoch[:, -1]
    a = get_time_series(tstart[0], cts_rate, tstop[-1], period = period, amp = amp, model = model)
    for i in range(len(tstart)):
        t.append(get_a_from_b(a,tstart[i],tstop[i]))
        # t.append(get_time_series(tstart[i], int(cts_rate * exp_time_epoch[i]) + 1, exp_time_epoch[i], period, amp,
        #                          model = 'sin'))
    t = list(chain.from_iterable(t))
    t = np.array(t)
    return t


t=get_epoch_ts_md('SgrA_I_epoch.txt',cts_rate=0.000617,period=554.,amp=0.99999,model='eclipse')
np.savetxt('test_e1.txt',t)
    # T=1e6
# period=55400.
# t=get_time_series(4e7,10000,T,period,0.5,'sin')
# print(t)
#
# #t=get_epoch_time_series(0.000617,554.,0.9,'sin')
#print(t)
# plt.hist(t,bins=int(bin),histtype = 'step')
# plt.show()