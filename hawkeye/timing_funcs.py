#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from stingray.events import EventList
import scipy.stats as stats
import random


def get_T_in_mbins(epoch_info,w,m,fi):
    T=2*np.pi/w
    T_in_perbin = np.zeros(m)
    # 每个bin的总积分时间
    tbin = T/m
    # 每个bin的时间长度
    if epoch_info.ndim==1:epoch_info=np.array([epoch_info])
    t_start=epoch_info[:,0];t_end = epoch_info[:, 1]

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
        else:
            T_in_perbin[np.mod(intN_bin_t_start[i],m)-1]+=(N_bin_t_end[i]-N_bin_t_start[i])*tbin
    return T_in_perbin

def choose_obs(epoch_info,flux_info=None,flux_filter=0,expT_filter=0,if_flux_high=True,if_expT_high=True,obsID=[]):
    if len(obsID)>0:
        filter=[]
        for i in range(len(obsID)):
            index=np.where(epoch_info[:,2].astype('int')==obsID[i])[0][0]
            filter.append(index)
        epoch_info_use = epoch_info[filter]
        useid = epoch_info_use[:, 2]
        return (useid,epoch_info_use)
    if len(epoch_info)!=len(flux_info):
        print('ERROR! Flux info!')
        return None
    if if_expT_high:
        if if_flux_high:
            filter=np.where((epoch_info[:, 3] > expT_filter)&(flux_info > flux_filter))
        else:
            filter=np.where((epoch_info[:, 3] > expT_filter)&(flux_info < flux_filter))
    else:
        if if_flux_high:
            filter=np.where((epoch_info[:, 3] < expT_filter)&(flux_info > flux_filter))
        else:
            filter=np.where((epoch_info[:, 3] < expT_filter)&(flux_info < flux_filter))

    epoch_info_use = epoch_info[filter]
    useid = epoch_info_use[:, 2]

    return (useid,epoch_info_use)

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

    return T

def filter_obs(src_evt,useid,bkg_evt=None):
    src_evt_use = src_evt[np.where(src_evt[:,-1] == useid[0])]
    i=1
    while i < len(useid):
        id=useid[i]
        src_evt_use_temp=src_evt[np.where(src_evt[:,-1]==id)]
        src_evt_use = np.concatenate((src_evt_use, src_evt_use_temp))
        i+=1
    if not bkg_evt:
        return src_evt_use
    if bkg_evt:
        bkg_evt_use = bkg_evt[np.where(bkg_evt[:-1] == useid[0])[0]]
        i = 1
        while i < len(useid):
            bkg_evt_use_temp = bkg_evt[np.where(bkg_evt[:-1] == id)[0]]
            bkg_evt_use = np.concatenate((bkg_evt_use, bkg_evt_use_temp))
            i+=1
        return (src_evt_use, bkg_evt_use)

def filter_time_t1t2(time,t1,t2):
    index1=np.where(time>(time[0]+t1))
    index2=np.where(time<(time[0]+t2))
    indexall=np.intersect1d(index1,index2)

    return time[indexall]


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

def get_hist_withbkg(t,t_bkg, len_bin,bkgscale=0.,tstart=0,tstop=0):
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
    lc_out.counts=lc_new.counts-bkgscale*lc_bkg.counts
    # lc_out.counts = lc_new.counts
    # print(lc_bkg.counts)
    return lc_out

def read_epoch(epoch_file):
    epoch=np.loadtxt(epoch_file)
    if epoch.ndim==1:
        epoch=np.array([epoch])
    # print(epoch)
    TSTART=epoch[:,0]
    TSTOP=epoch[:,1]
    OBSID=epoch[:,2]
    exptime=epoch[:,3]
    return (TSTART, TSTOP, OBSID, exptime)


def get_Z2(time,freq):
    N=len(time)
    def turns(t,f):
        ti=t-t[0]
        v=f
        p=1.0/v
        # pdot=-vdot/(v*v)
        # vddot = 2.0 * pdot * pdot / (p * p * p)
        freq = v # + ti * vdot + ti * ti / 2.0 * vddot
        turns=v*ti
        INT_turns=np.trunc(turns)
        turns=turns-INT_turns
        turns = 2.0*np.pi*turns #+ vdot * ti * ti / 2.0 + vddot * ti * ti * ti / 6.0
        return turns
    #print time
    Z2=[]

    for fi in freq:
        Z2.append((2.0 / N)*(sum(np.cos(turns(time,fi)))**2+sum(np.sin(turns(time,fi))**2)))
    cts=len(time)
    return (Z2,cts)

def plot_Z2(time,freq):
    T_tot=1e5
    freq=np.arange(1/T_tot,0.5/100,1/(5*T_tot))
    freq=freq[np.where(freq>1/10000)]
    v=len(freq)
    [Z2,cts]=get_Z2(time,freq)
    Z2=np.array(Z2)
    (p99_Z2,p90_Z2)=get_Z2_thres(v)
    plt.figure(1,(7,7))
    plt.semilogx(freq,[p99_Z2 for i in range(len(Z2))],'--',color='black')
    plt.semilogx(freq,[p90_Z2 for i in range(len(Z2))],'--',color='red')
    plt.plot([1/1754.38596,1/1754.38596],[0,50],'--')
    plt.step(freq,Z2,color='black')
    plt.text(0.0005,p99_Z2+1,"99%")
    plt.text(0.0005,p90_Z2+1,"90%")
    plt.show()

def get_Z2_thres(v):
    # 采样频率数
    v = 3e7
    ##(1-p)**N=0.99
    # 置信度99%
    p99 = 1 - 0.99 ** (1.0 / v)
    p90 = 1 - 0.90 ** (1.0 / v)
    x1 = []
    y1 = []
    x2 = []
    y2 = []
    for x in np.linspace(1, 1000, 5000):
        if stats.chi2.cdf(x, 2) != 0:
            y1.append(np.abs((1 - stats.chi2.cdf(x, 2)) - p99))
            y2.append(np.abs((1 - stats.chi2.cdf(x, 2)) - p90))
            # print y1
            x1.append(x)
            x2.append(x)
        else:
            break
    y1 = np.array(y1)
    x1 = np.array(x1)
    y2 = np.array(y2)
    x2 = np.array(x2)
    p99_Z2=x1[np.where(y1 == min(y1))][0]
    p90_Z2=x1[np.where(y2 == min(y2))][0]
    # print min(y1)
    # print(x1[np.where(y1 == min(y1))])
    # print("%e" % (1 - stats.chi2.cdf(x1[np.where(y1 == min(y1))], 2)))
    # print(p99)
    # print min(y2)
    # print(x1[np.where(y2 == min(y2))])
    # print("%e" % (1 - stats.chi2.cdf(x1[np.where(y2 == min(y2))], 2)))
    # print(p90)
    return (p99_Z2,p90_Z2)