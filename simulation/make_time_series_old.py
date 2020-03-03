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
import random
from scipy import optimize
from sympy import *
from scipy.optimize import fsolve
from itertools import chain
def get_time_series(T_START,N_cts,T_stop,period,amp,model,eqw=0.1):
    #N_cts为总光子数
    #T为总观测时间
    #period为周期
    #amp为model的amplitude,反应了时变强度
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


def get_epoch_time_series(cts_rate,period,amp, model):
    t=[]
    # path='/Users/baotong/Desktop/period/'
    # epoch = np.loadtxt(path + 'txt_all_obs_I/' + 'SgrA_I_epoch.txt')
    epoch = np.loadtxt('LW_epoch.txt')
    tstart=epoch[:,0]
    tstop=epoch[:,1]
    exp_time_epoch=epoch[:,-1]
    for i in range(len(tstart)):
        cts = int(cts_rate * (tstop[i] - tstart[i]))
        t.append(get_time_series(tstart[i],cts,tstop[i],period,amp,model=model))
    t=list(chain.from_iterable(t))
    t = np.array(t)
    return t

# t=get_epoch_time_series(cts_rate=0.000617,period=554.,amp=0.99, model='sin')
# np.savetxt('old.txt',t)
# print(t)


# All_T=1618692.94
# cts_rate=0.0006

def get_Z2(time,freq):
    #time=   np.loadtxt(path + 'txt/' + dataname + '.txt')[:,0]
    #time = np.loadtxt(path + 'txt/' + dataname + '_s.txt')[:, 0]
    #time = np.loadtxt(path + 'txt_G/' + dataname + '.txt')[:, 0]
    #time=np.loadtxt(path + 'txt/evt_list10956.txt')[:,0]
    N=len(time)
    def turns(t,f):
        ti=t-t[0]
        #print ti
        v=f
        #p_test = 1.0/5.500550055005501e-06
        p=1.0/v
        # pdot=-vdot/(v*v)
        # vddot = 2.0 * pdot * pdot / (p * p * p)
        freq = v # + ti * vdot + ti * ti / 2.0 * vddot
        preq = 1.0 / freq
        turns=v*ti
        INT_turns=np.trunc(turns)
        turns=turns-INT_turns
        turns = 2.0*np.pi*turns #+ vdot * ti * ti / 2.0 + vddot * ti * ti * ti / 6.0
        #print turns
        #turns=list(turns)
        #turns=np.array(turns[0])
        #初始相位
        return turns
    #print time
    Z2=[]
    #plt.hist(time-time[0],bins=100,histtype='step')
    #plt.show()

    for fi in freq:
        Z2.append((2.0 / N)*(sum(np.cos(turns(time,fi)))**2+sum(np.sin(turns(time,fi))**2)))
    cts=len(time)
    return [Z2,cts]

def plot_amp():
    period=554.
    amp=1.
    cts_rate_temp=0.000617
    border = 10000
    vary = np.array([i for i in range(0, border)])
    freq = 1 / 3000. + vary * 1e-7
    y=[]
    period_detect=[]
    x=np.linspace(1,20,20)
    t = get_epoch_time_series(cts_rate = i * cts_rate_temp, period = period, amp = amp, model = 'sin')
    plt.hist(t,bins=20)
    plt.show()

    for i in x:
        #t=get_time_series(T_START=0,N_cts = int(i*100),T=1e6,period=period,amp=0.5,model='sin')
        t=get_epoch_time_series(cts_rate=i*cts_rate_temp,period=period,amp=amp, model='sin')
        Z2=get_Z2(t,freq)[0]
        y.append(np.max(Z2))
        period_detect.append(1./freq[np.where(Z2==np.max(Z2))])
    p99=23.00
    All_T=1618692.94
    cts=x*cts_rate_temp*All_T

    plt.figure(1,(8,6))
    plt.semilogy()
    plt.subplot(211)
    plt.title('amp={0}'.format(amp))
    plt.plot(cts,np.zeros(len(cts))+p99,'--',color='r')
    plt.xlabel('cts')
    plt.ylabel('Z2')
    plt.scatter(cts,y)
    plt.subplot(212)
    plt.plot(cts,np.zeros(len(cts))+period,'--',color='r')
    plt.scatter(cts,period_detect)
    plt.xlabel('cts')
    plt.ylabel('detect_p')
    plt.savefig('/Users/baotong/Desktop/li_pr/11-15/50k/amp_0.eps')
    plt.show()

#plot_amp()
def plot_Z2():
    period = 5540.
    amp = 0.99
    cts_rate_temp = 0.0006171
    border = 5000
    vary = np.array([i for i in range(0, border)])
    freq = 1 /7000. + vary * 1e-8

    # border=10000
    # vary=np.array([i for i in range(0,border)])
    # freq=1.e-4+vary*1.e-6
    p99_Z2 = 27.6244
    p90_Z2 = 22.9232

    t = get_epoch_time_series(cts_rate =cts_rate_temp, period = period, amp = amp, model = 'sin')
    np.savetxt('old.txt',t)
    #t = np.loadtxt('test_1.txt')
    Z2=get_Z2(t,freq)[0]
    cts=get_Z2(t,freq)[1]
    Z2=np.array(Z2)
    plt.figure(1,(7,7))
    plt.title('sim'+',cts={0}'.format(cts))
    plt.semilogx(freq,[p99_Z2 for i in range(len(Z2))],'--',color='black')
    plt.semilogx(freq,[p90_Z2 for i in range(len(Z2))],'--',color='red')
    plt.step(freq,Z2,color='black')
    plt.text(0.0005,p99_Z2+1,"99%")
    plt.text(0.0005,p90_Z2+1,"90%")
    plt.show()
#plot_Z2()