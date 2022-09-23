#!/bin/bash
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 19 18:13:40 2022
@author: baotong
from real data to simulate lightcurve, eventlist
by best-fit psd or const, plus sine light curve
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
import stingray as sr
from stingray.events import EventList
from stingray.lightcurve import Lightcurve
from stingray import Lightcurve, Crossspectrum, sampledata,Powerspectrum,AveragedPowerspectrum
from stingray.simulator import simulator, models
import rednoise
from rednoise import useful_functions as func
from scipy import integrate
import hawkeye as hawk
import rednoise

def sinfunc(t,w,A,fi,B=0):
    y=A*np.sin(w*t+fi)+B
    return y


def simulate_srcinfo(srcid):
    lenbin=100.0
    figurepath='/Users/baotong/Desktop/aas/pXS_Tuc/figure/rednoise/'
    srcpsdinfo=np.loadtxt(figurepath+str(srcid)+'_2738.txt')
    if srcpsdinfo[1]<0:srcpsdinfo[1]=1e-3
    (src_evt_use, epoch_info_use) = rednoise.load_data(srcid,ifobsid=[2738])
    w = rednoise.make_freq_range(dt=lenbin, epoch_info=epoch_info_use)
    w=w[w>1/50000.]
    lc_use = func.get_hist(src_evt_use[:,0], len_bin=100, tstart=epoch_info_use[:,0][0], tstop=epoch_info_use[:,1][-1])
    psd_org=Powerspectrum(lc_use)
    (frms,frms_err)=psd_org.compute_rms(min_freq=w[0], max_freq=w[-1], white_noise_offset=0)
    print('frms=',frms)
    (frms,frms_err)=sr.excess_variance(lc_use,normalization='fvar')
    print('frms=', frms)
    frms=(lc_use.counts.std()/lc_use.counts.mean())
    print('frms=', frms)
    # epoch_info_use=np.array(epoch_info_use)
    # psd_model = build_psd(x=w, p=srcpsdinfo, type='smoothbkpl')
    psd_model = rednoise.build_psd(x=w, p=srcpsdinfo[:-1], type='breakp')
    CR=len(src_evt_use)/epoch_info_use[:,3]
    lc = rednoise.make_lc_from_psd(psd=psd_model, cts_rate=CR * lenbin, dt=lenbin, epoch_info=epoch_info_use)
    # print('counts={0}'.format(np.sum(lc.counts)))
    lc_evt = Lightcurve(time=lc.time, counts=lc.counts, dt=lc.dt, gti=lc.gti)
    lc_evt.counts = np.random.poisson(lc_evt.counts)
    psd=rednoise.plot_psd(lc,show=0)
    plt.close()

    # rednoise.plot_psd(lc_evt,show=1)
    (src_evt_use, epoch_info_use) = rednoise.load_data(srcid,ifobsid=[953,955,2736,2738,16527,15747,16529,17420,15748,16528])
    CR = len(src_evt_use) / np.sum(epoch_info_use[:,3])
    print(CR)


    for i in range(len(epoch_info_use)):
        t_src = src_evt_use[:, 0][np.where(src_evt_use[:, 2] == epoch_info_use[i][2])]
        epoch_temp=np.array([epoch_info_use[i]])
        if len(t_src)<2 or epoch_info_use[i][3]<lenbin:
            continue
        lc_new = rednoise.make_lc_from_psd(psd=psd_model, cts_rate=CR * lenbin, dt=lenbin, epoch_info=epoch_temp,frms=0.24)
        evt = EventList()
        EventList.simulate_times(evt,lc=lc_new)
        # lc_evt = Lightcurve(time=lc.time, counts=lc.counts, dt=lc.dt, gti=lc.gti)
        # lc_evt.counts = np.random.poisson(lc_evt.counts)
        # if i == 0:
        #     lc_all = lc_evt
        # else:
        #     lc_all = lc_all.join(lc_evt)
        if i==0:
            evt_all=evt
        else:
            evt_all=evt_all.join(evt)

    return evt_all

def sim_onesrc_psdandsin():
    srcid=312


