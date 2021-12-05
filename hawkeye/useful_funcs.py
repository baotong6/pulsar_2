#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
from stingray.events import EventList

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
            #print(rest)
        else:
            T_in_perbin[np.mod(intN_bin_t_start[i],m)-1]+=(N_bin_t_end[i]-N_bin_t_start[i])*tbin
    return T_in_perbin

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
    src_evt_use = src_evt[np.where(src_evt[:-1] == useid[0])[0]]
    i=1
    while i < len(useid):
        id=useid[i]
        src_evt_use_temp=src_evt[np.where(src_evt[:-1]==id)[0]]
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