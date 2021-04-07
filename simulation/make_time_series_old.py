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
# from sympy import *
from scipy.optimize import fsolve
from itertools import chain
from astropy.timeseries import LombScargle

def get_T_in_mbins(epoch_file,w,m,fi):
    T=2*np.pi/w
    T_in_perbin = np.zeros(m)
    # 每个bin的总积分时间
    tbin = T/m
    # 每个bin的时间长度
    epoch_info = np.loadtxt(epoch_file)
    t_start = epoch_info[:, 0]
    t_end = epoch_info[:, 1]
    ID = epoch_info[:, 2]
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

def pfold(time,p_test,shift,bin=20):
    def trans(t,p_test,shift):
      ti =t
      v = 1.0 /p_test
      turns = v * ti
      turns += shift
      # 初始相位
      for i in range(len(turns)):
        turns[i] = turns[i] - int(turns[i])
      return turns
    epoch_file='CDFS_epoch.txt'
    T_in_perbin = get_T_in_mbins(epoch_file, 2 * np.pi / p_test, bin, shift * 2 * np.pi)
    turns=trans(time,p_test,shift)
    loc=np.zeros(bin)
    for index in turns:
      loc[int(index*bin)] += 1

    x = np.array([(i / bin + 0.5 / bin) for i in range(bin)])

    x2=np.concatenate((x,x+1))
    y2=np.concatenate((loc,loc))

    correct_gap = T_in_perbin / (sum(T_in_perbin) / len(T_in_perbin))
    y2 /= np.concatenate((correct_gap, correct_gap))
    plt.figure(1,(8,8))

    plt.xlabel('phase')
    plt.ylabel('counts/bin')
    plt.step(x2,y2,color='red')
    plt.show()

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
            cts = int((cts_rate * (tstop[i] - tstart[i]))*(np.random.random(1)[0]*3.16227766+0.316227766))
            t.append(get_time_series(tstart[i], cts, tstop[i], period, amp, model=model))
    else:
        for i in range(len(tstart)):
            cts = int(cts_rate * (tstop[i] - tstart[i]))
            t.append(get_time_series(tstart[i],cts,tstop[i],period,amp,model=model))
    t=list(chain.from_iterable(t))
    t = np.array(t)
    t = t + np.random.rand(len(t)) * 3.2 - 1.6
    return t

# pfold(time_s,893.,0.5)

def get_hist(t, len_bin):
    ###将输入的time信息，按照len_bin的长度输出为lc
    t_test = t-t[0]
    a = [0 for i in range(int(t_test[-1] / len_bin) + 1)]
    for i in range(len(t_test)):
        a[int(t_test[i] / len_bin)] += 1
    a = np.array(a)
    return a
def get_LS(time, flux,freq):
    def get_hist(t, len_bin):
        ###将输入的time信息，按照len_bin的长度输出为lc
        t_test = t - t[0]
        a = [0 for i in range(int(t_test[-1] / len_bin) + 1)]
        for i in range(len(t_test)):
            a[int(t_test[i] / len_bin)] += 1
        a = np.array(a)

    x = np.arange(bin_len / 2., (time[-1] - time[0]) + bin_len / 2., bin_len)
    y = flux
    LS = LombScargle(x, y,normalization = 'standard')
    # LS = LombScargle(x, y, dy, normalization='psd')
    power = LS.power(freq)
    FP=LS.false_alarm_probability(power.max(),samples_per_peak=5, nyquist_factor=5,minimum_frequency = freq[0], maximum_frequency = freq[-1],method='baluev')
    FP_99 = LS.false_alarm_level(0.01, samples_per_peak=10, nyquist_factor=5,minimum_frequency = freq[0], maximum_frequency = freq[-1],method='naive')
    FP_90 = LS.false_alarm_level(0.1, samples_per_peak=10, nyquist_factor=5, minimum_frequency=freq[0],
                                 maximum_frequency=freq[-1], method='naive')
    FP_68 = LS.false_alarm_level(0.32, samples_per_peak=10, nyquist_factor=5, minimum_frequency=freq[0],
                                 maximum_frequency=freq[-1], method='naive')

    plt.title('{0},FP={1}'.format('test',FP))

    plt.plot(freq, power)
    print(1./freq[np.where(power==np.max(power))])
    plt.plot([freq[0], freq[-1]], [FP_99, FP_99], '--')
    plt.plot([freq[0], freq[-1]], [FP_90, FP_90], '--')
    plt.plot([freq[0], freq[-1]], [FP_68, FP_68], '--')
    plt.show()

border = 500000
vary = np.array([i for i in range(0, border)])
freq = 1 /50000 + vary * 1.e-8
bin_len=100.
time_s = get_epoch_time_series(1e-4, 893, 0.7, 'sin')
print(len(time_s))
flux = get_hist(time_s, bin_len)
get_LS(time_s,flux,freq)

def get_sim_time_series(dataID):
    path='/Users/baotong/Desktop/period/txt_all_obs_I/'
    time_orig=np.loadtxt(path+str(dataID)+'.txt')
    epoch=np.loadtxt(path+'epoch_src_'+str(dataID)+'.txt')
    tstart=epoch[:,0]
    tstop=epoch[:,1]
    exp_time_epoch=epoch[:,-1]
    t=[];S=[]
    for i in range(len(epoch[:,2])):
        index=np.where(time_orig==epoch[:,2][i])
        #cts_num=len(time_orig[index])
        cts_num=exp_time_epoch[i]*len(time_orig)/np.sum(exp_time_epoch)
        CTS_pois = np.random.poisson(lam=int(cts_num), size=1)[0]
        t.append(get_time_series(tstart[i], CTS_pois, tstop[i], model='const'))
        S.append(CTS_pois/exp_time_epoch[i])
    S=np.array(S)
    S_del = np.delete(S, np.where(S <1e-15))
    VI=np.max(S_del)/np.min(S_del)
    t=list(chain.from_iterable(t))
    t = np.array(t)
    t = t + np.random.rand(len(t)) * 3.2 - 1.6
    return [t,VI]
# pfold(get_sim_time_series(2238)[0],44529.,0.5)
# def plot_test():
#     VI_plt=[]
#     for i in range(200):
#         [t,VI]=get_sim_time_series('2338')
#         VI_plt.append(VI)
#     VI_plt=np.array(VI_plt)
#     plt.hist(VI_plt,bins=20)
#     plt.show()
# plot_test()
def get_VI_hist(VI_input=10):
    VI=[]
    epoch = np.loadtxt('SgrA_I_epoch.txt')
    # epoch = np.loadtxt('LW_epoch.txt')
    expT=epoch[:,-1]
    print(expT)
    R=5e-4
    for i in range(6400):
        S=[]
        for j in range(49):
            CTS= np.random.poisson(lam=int(R*expT[j]), size=1)
            S.append((CTS+0.0)/expT[j])
        S=np.array(S)
        S=np.delete(S,np.where(S==0.0))
        #print(np.min(S))
        # S=[]
        # for i in range(49):
        #     S.append(np.random.random(1)[0] * VI_input**0.5 + 1./VI_input**0.5)
        VI.append(np.max(S)/np.min(S))
    print(VI)
    ax=plt.hist(VI,bins=15,histtype='step')
    print(ax)
    plt.show()
#get_VI_hist()

def get_var_W(P,w_range):
    TOT_N=100
    ni = 10;
    m = 12;
    MEAN_OUT=[]
    Om1w = np.zeros(m)
    for i in range(TOT_N):
        t=get_epoch_time_series(cts_rate=0.0001,period=P,amp=0.7, model='sin')
        N=len(t)
        fa = np.zeros(m_max)
        for m in range(0, m_max):  # precompute factors for each m
            fa[m] = compute_factor(N, m + 1, v)

        fbin = precompute_binmult(N)
        Tlist=t;m=10;w=2*np.pi/w_range;fi=0.
        n = compute_bin(Tlist, m, w, fi)  # find bin histogram for a given bin number m, frequency w and fi
        f = 0
        for i in range(0, m):
            # f=f+np.sum(np.log(np.arange(2,n[i]+1)))
            f = f + fbin[n[i]]

        y = np.exp(f + factor)
        MEAN_OUT.append(f)
        O_out.append(Om1w)


    MEAN_OUT=np.array(MEAN_OUT)
    O_out=np.array(O_out)
    var=np.var(MEAN_OUT)
    return [var,O_out]

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