#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from datetime import datetime
from scipy.interpolate import lagrange
import pandas as pd
from scipy import optimize
# from sympy import *
from scipy.optimize import fsolve
from itertools import chain
from astropy.stats import poisson_conf_interval
import scipy
import hawkeye as hawk
import rednoise as rednoise
import os 

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
def make_nustar_epoch_file(srcname,exp):
    tstart=np.random.random(1)[0]*100;id=1
    epoch_info=[];exptemp=0
    while exptemp<exp:
        if exptemp+3200<exp:
            epoch_info.append([tstart,tstart+3200,id,3200])
            exptemp+=3200
        else:
            epoch_info.append([tstart,tstart+exp-exptemp,id,exp-exptemp])
            exptemp+=exp-exptemp
        tstart+=95*60;id+=1
    epoch_info=np.array(epoch_info)
    return epoch_info
def get_epoch_time_series(cts_rate,period,amp, model,epochname=None):
    t=[]
    path='/Users/baotong/nustar/simulation/'
    epoch = np.loadtxt(path+epochname)
    # epoch = np.loadtxt('LW_epoch.txt')
    tstart=epoch[:,0]
    tstop=epoch[:,1]
    exp_time_epoch=epoch[:,-1]
    if model=='const':
        for i in range(len(tstart)):
            cts = np.random.poisson(int((cts_rate * (tstop[i] - tstart[i]))))
            t.append(get_time_series(tstart[i], cts, tstop[i], period, amp, model=model))
    else:
        for i in range(len(tstart)):
            cts = np.random.poisson(int(cts_rate * (tstop[i] - tstart[i])))
            t.append(get_time_series(tstart[i],cts,tstop[i],period,amp,model=model))
    t=list(chain.from_iterable(t))
    t = np.array(t)
    t = t + np.random.rand(len(t)) * 3.2 - 1.6
    return t

def read_plot_simres():
    path='/Users/baotong/nustar/simulation/result_nustar_ft0.1/'
    res1=np.loadtxt(path+'result_src1_{0}.txt'.format(0))
    for i in range(1,99):
        if os.path.exists(path+'result_src1_{0}.txt'.format(i)):
            res_single1=np.loadtxt(path+'result_src1_{0}.txt'.format(i))
            res1=np.row_stack((res1,res_single1))
    res2=np.loadtxt(path+'result_src2_{0}.txt'.format(0))
    for i in range(1,99):
        if os.path.exists(path+'result_src2_{0}.txt'.format(i)):
            res_single2=np.loadtxt(path+'result_src2_{0}.txt'.format(i))
            res2=np.row_stack((res2,res_single2))

    res3=np.loadtxt(path+'result_src3_{0}.txt'.format(0))
    for i in range(1,99):
        if os.path.exists(path+'result_src3_{0}.txt'.format(i)):
            res_single3=np.loadtxt(path+'result_src3_{0}.txt'.format(i))
            res3=np.row_stack((res3,res_single3))
    period1=np.array(res1[:,3]);P1=res1[:,2]
    period2=np.array(res2[:,3]);P2=res2[:,2]
    period3=np.array(res3[:,3]);P3=res3[:,2]
    period1=period1[np.where(P1>0.99)]
    period2=period2[np.where(P2>0.99)]
    period3=period3[np.where(P3>0.99)]
    print(len(period1))
    plt.figure(1,(10,8))
    plt.hist(np.abs(2*np.pi/period1/354.6-1),bins=4,histtype='step',lw=2,label='Source 1,DR={0}%'.format(len(np.where(P1>0.99)[0])))
    plt.hist(np.abs(2*np.pi/period2/4093.3-1),bins=10,histtype='step',lw=2,label='Source 2,DR={0}%'.format(len(np.where(P2>0.99)[0])))
    plt.hist(np.abs(2*np.pi/period3/607.5-1),bins=5,histtype='step',lw=2,label='Source 3,DR={0}%'.format(len(np.where(P3>0.99)[0])))
    print(np.mean(np.abs(2*np.pi/period1-354.6)))
    print(np.mean(np.abs(2*np.pi/period2-4093.3)))
    print(np.mean(np.abs(2*np.pi/period3-607.5)))
    print(np.max(np.abs(2*np.pi/period1-354.6)))
    print(np.max(np.abs(2*np.pi/period2-4093.3)))
    print(np.max(np.abs(2*np.pi/period3-607.5)))
    plt.semilogy()
    plt.semilogx()
    plt.xlabel(r'|$P_{det}/P_{true}-1$|',hawk.font1)
    plt.tick_params(labelsize=16)
    plt.ylabel('Number of simulations',hawk.font1)
    plt.legend(fontsize=12)
    plt.savefig(path+'res.pdf',bbox_inches='tight', pad_inches=0.05)
    plt.show()
    # res=np.column_stack((srcid,res))
    # np.savetxt(path+ 'all_result_20ks.txt', res, fmt='%10d %10.2f %10.5f %10.10f %10.5f %10d %10.10f %10.10f %10d')

if __name__=='__main__':
    # epoch=make_nustar_epoch_file(1,70000)
    # np.savetxt('/Users/baotong/nustar/simulation/epoch_src_1.txt',epoch,fmt='%10d %10d %5d %10d')
    read_plot_simres()
    # for i in range(100):
    #     t=get_epoch_time_series(cts_rate=400/50000,period=4093.3,amp=0.348, model='sin',epochname='epoch_src_2.txt')
    #     energy=np.zeros(len(t))+10000
    #     evt=np.column_stack((t,energy))
    #     np.savetxt('/Users/baotong/nustar/simulation/evt_ft0.1_2/src2_{0}.txt'.format(i),evt,fmt='%10.5f %10d')
