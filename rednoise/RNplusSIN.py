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

def sinfunc(x, w,fi,A, B):
    return A * np.sin(w*x+fi) + B

def get_turns(t,period):
    v=1.0/period
    turns=t*v-int(t*v)
    return turns


def sim_onesrc_psdandsin(srcid,period=None,ifobsid=None,amp=None):
    dt=1000
    figurepath='/Users/baotong/Desktop/aas/pXS_Tuc/figure/rednoise/'
    srcpsdinfo=np.loadtxt(figurepath+str(srcid)+'_2738.txt')
    if srcpsdinfo[1]<0:srcpsdinfo[1]=1e-3
    (src_evt_use, epoch_info_use) = rednoise.load_data(srcid,ifobsid=[2738])
    w = rednoise.powerlw2evt.make_freq_range(dt=dt, epoch_info=epoch_info_use)
    w=w[w>1/50000.]

    psd_model = rednoise.powerlw2evt.build_psd(x=w, p=srcpsdinfo[:-1], type='breakp')

    (src_evt_use, epoch_info_use) = rednoise.load_data(srcid,ifobsid=ifobsid)
    CR = len(src_evt_use) / np.sum(epoch_info_use[:,3])
    for i in range(len(epoch_info_use)):
        t_src = src_evt_use[:, 0][np.where(src_evt_use[:, 2] == epoch_info_use[i][2])]
        epoch_temp=np.array([epoch_info_use[i]])
        if len(t_src)<2 or epoch_info_use[i][3]<dt:
            # print(len(t_src),exptime[i])
            continue
        lc_new = rednoise.powerlw2evt.make_lc_from_psd(psd=psd_model, cts_rate=CR * dt, dt=dt, epoch_info=epoch_temp,frms=0.2)
        lc_evt = Lightcurve(time=lc_new.time, counts=lc_new.counts, dt=lc_new.dt, gti=lc_new.gti)
        addsin = sinfunc(lc_evt.time, 2 * np.pi / period, fi=0, A=CR*amp * dt, B=0)
        lc_evt.counts = lc_evt.counts + addsin
        lc_evt.counts[np.where(lc_evt.counts < 0)] = 0
        ##这里假设amp=0.5，则high/low~3
        lc_evt.counts = np.random.poisson(lc_evt.counts)
        evt = EventList()
        EventList.simulate_times(evt,lc=lc_evt)
        if i == 0:
            lc_all = lc_evt
        else:
            lc_all = lc_all.join(lc_evt)

        if i == 0:
            evt_all = evt
        else:
            evt_all = evt_all.join(evt)

    return (evt_all, lc_all)

def sim_onesrc_psdandeclipse(srcid,period,ifobsid=None,width=0.1,amp=0.9):
    dt=1000
    figurepath='/Users/baotong/Desktop/aas/pXS_Tuc/figure/rednoise/'
    srcpsdinfo=np.loadtxt(figurepath+str(srcid)+'_2738.txt')
    if srcpsdinfo[1]<0:srcpsdinfo[1]=1e-3
    (src_evt_use, epoch_info_use) = rednoise.load_data(srcid,ifobsid=[2738])
    w = rednoise.powerlw2evt.make_freq_range(dt=dt, epoch_info=epoch_info_use)
    w=w[w>1/50000.]
    psd_model = rednoise.powerlw2evt.build_psd(x=w, p=srcpsdinfo[:-1], type='breakp')
    (src_evt_use, epoch_info_use) = rednoise.load_data(srcid,ifobsid=ifobsid)
    CR = len(src_evt_use) / np.sum(epoch_info_use[:,3])
    for i in range(len(epoch_info_use)):
        t_src = src_evt_use[:, 0][np.where(src_evt_use[:, 2] == epoch_info_use[i][2])]
        epoch_temp=np.array([epoch_info_use[i]])
        if len(t_src)<2 or epoch_info_use[i][3]<dt:
            # print(len(t_src),exptime[i])
            continue
        lc_new = rednoise.powerlw2evt.make_lc_from_psd(psd=psd_model, cts_rate=CR * dt, dt=dt, epoch_info=epoch_temp,frms=0.2)
        lc_evt = Lightcurve(time=lc_new.time, counts=lc_new.counts, dt=lc_new.dt, gti=lc_new.gti)
        lc_evt.counts[np.where(lc_evt.counts < 0)] = 0
        lc_evt.counts = np.random.poisson(lc_evt.counts)
        evt = EventList()
        EventList.simulate_times(evt,lc=lc_evt)
        if i == 0:
            lc_all = lc_evt
        else:
            lc_all = lc_all.join(lc_evt)
        if i == 0:
            evt_all = evt
        else:
            evt_all = evt_all.join(evt)
    t_eclipse=evt_all.time

    i=0
    while i < len(t_eclipse):
        if (0.5 - 0.5 * width) < get_turns(t_eclipse[i], period) < (0.5 + 0.5 * width):
            rand = np.random.random(1)[0]
            if rand < amp:
                t_eclipse = np.delete(t_eclipse, i)
            else:
                i += 1
        else:   i += 1
    evt_all=EventList()
    evt_all.time=t_eclipse
    return evt_all
if __name__ == '__main__':
    (evt_all, lc_all)=sim_onesrc_psdandsin(srcid='414',period=8646.78,amp=0.5)
    t_eclipse=sim_onesrc_psdandeclipse(srcid='366',period=10311.94,width=0.1,amp=0.9)
    print(len(t_eclipse))
    # hawk.get_LS(time=lc_all.time,flux=lc_all.counts,freq=np.arange(1/10000.,1/3000.,1e-7),show=1)
    # plt.show()
    # dt=5000
    # # first test ##
    # time=np.arange(0,100000,dt)
    # lc=Lightcurve(time=time,counts=np.zeros(len(time))+0.0038*dt)
    #
    # addsin=sinfunc(time,2*np.pi/48780.49,fi=0,A=0.0038*0.2*dt,B=0)
    # lc.counts=lc.counts+addsin
    # lc.counts = np.random.poisson(lc.counts)
    # lc.plot()
    # plt.show()
    #
    # evt = EventList()
    # EventList.simulate_times(evt, lc=lc)
    # lc_new=evt.to_lc(dt=dt)
    # lc_new.plot()
    # plt.show()