#!/bin/bash
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 19 18:13:40 2022
@author: baotong
from real data to simulate lightcurve, eventlist
by best-fit psd or const
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.optimize import curve_fit
import astropy.units as u
import astropy.constants as c
import stingray as sr
from stingray.events import EventList
from stingray.lightcurve import Lightcurve
from stingray import Lightcurve, Crossspectrum, sampledata,Powerspectrum,AveragedPowerspectrum
from stingray.simulator import simulator, models
import rednoise
from rednoise import useful_functions as func
import hawkeye as hawk
from scipy import integrate

font1=func.font1
font2=func.font2

def build_psd(x,p,x_2=1,p_2=1,type='bendp'):
    ##
    if type=='bendp':
        psd=func.bending_po(x,p)
    elif type=='lorentz':
        psd=func.generalized_lorentzian(x,p)
    elif type=='bendp+lorentz':
        psd=func.bending_po(x,p)+func.standard_lorentzian(x_2,p_2)
    elif type=='breakp':
        psd=func.break_po(x,p)
    elif type=='smoothbkpl':
        psd=func.smoothbkpl(x,p)
    return psd

def make_freq_range(dt,epoch_info):
    if type(epoch_info)==np.ndarray:epoch_info=np.array(epoch_info)
    else:epoch_info=np.loadtxt(epoch_info)
    if epoch_info.ndim == 1:
        epoch_info=np.array([epoch_info])
    TSTART=epoch_info[:,0];TSTOP=epoch_info[:,1];OBSID=epoch_info[:,2];exptime=epoch_info[:,3]
    T_tot=TSTOP[-1]-TSTART[0]
    w = np.arange(2 / T_tot, 0.5 / dt, 1 / T_tot)
    return w

def make_lc_from_psd(psd,cts_rate,dt,epoch_info,frms=1):
    if type(epoch_info)==np.ndarray:epoch_info=np.array(epoch_info)
    else:epoch_info=np.loadtxt(epoch_info)
    if epoch_info.ndim == 1:
        epoch_info=np.array([epoch_info])
    TSTART=epoch_info[:,0];TSTOP=epoch_info[:,1];OBSID=epoch_info[:,2];exptime=epoch_info[:,3]
    T_tot=TSTOP[-1]-TSTART[0]
    num_bins=int(T_tot/dt)
    sim = simulator.Simulator(N=num_bins + 1, mean=cts_rate,dt=dt,rms=frms)
    # lc=sim.simulate(1)
    lc=sim.simulate(psd)
    lc.time=lc.time+TSTART[0]
    lc.counts[np.where(lc.counts < 0)] = 0
    lc.gti=[[lc.time[0],lc.time[-1]]]
    return lc

'''
def test_simple():
    lenbin=1
    T_exp=113468.73
    CR=5e-2
    epoch_89='/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep3/epoch_src_test_single.txt'
    (TSTART, TSTOP, OBSID, exptime) = func.read_epoch(epoch_89)
    epoch_89=(TSTART, TSTOP, OBSID, exptime)
    w=make_freq_range(dt=lenbin,epoch_info=epoch_89)

    p_1=[10,1e-3,1,2,1/4]
    psd_model=build_psd(x=w,p=p_1,type='smoothbkpl')
    frms = integrate.quad(func.smoothbkpl, w[0], w[-1], args=(p_1))[0]
    frms = np.sqrt(frms)
    print(frms)
    plt.loglog()
    plt.plot(w, psd_model)
    plt.xlabel('Frequency (Hz)', font2)
    plt.ylabel('Power', font2)
    plt.tick_params(labelsize=16)
    plt.show()
    lc = make_lc_from_psd(psd=psd_model, cts_rate=CR * lenbin, dt=lenbin, epoch_file=epoch_89, frms=0.22)
    ps_org = rednoise.plot_psd(lc)
    print('counts={0}'.format(np.sum(lc.counts)))
    lc_evt = Lightcurve(time=lc.time, counts=lc.counts, dt=lc.dt, gti=lc.gti)
    lc_evt.counts = np.random.poisson(lc_evt.counts)
    # ev_all = EventList()
    # ev_all.time = func.sim_evtlist(lc)
    # lc_evt = ev_all.to_lc(dt=lenbin, tstart=ev_all.time[0] - 0.5 * lenbin,
    #                       tseg=ev_all.time[-1] - ev_all.time[0])
    ps_real = rednoise.plot_psd(lc_evt)
    freq = np.arange(1 / T_exp, 0.5 / lenbin, 1 / (5 * T_exp))
    freq = freq[np.where(freq > 1 / 10000.)]
    # temp = func.get_LS(lc_evt.time, lc_evt.counts, freq=freq)
    # plt.subplot(121)
    plt.xlabel('Time', font2)
    plt.ylabel('Counts/bin', font2)
    plt.tick_params(labelsize=16)
    plt.plot(lc.time, lc.counts, color='red')
    # plt.subplot(122)
    plt.show()
    plt.plot(lc_evt.time, lc_evt.counts, color='green')
    print('counts={0}'.format(np.sum(lc_evt.counts)))
    plt.xlabel('Time', font2)
    plt.ylabel('Counts/bin', font2)
    plt.tick_params(labelsize=16)
    plt.show()
    evt = EventList()
    EventList.simulate_times(evt, lc=lc)
    print('counts={0}'.format(len(evt.time)))

    rednoise.optimize_psdmodel(ps_org, whitenoise=100, show=1, label='test', save=0, figurepath=None, outfigname=None)
    return None
'''

def simulate_srcinfo(srcid,ifobsid=None,path_provide=None):
    lenbin=100.0
    figurepath='/Users/baotong/Desktop/aas/pXS_Tuc/figure/rednoise/'
    srcpsdinfo=np.loadtxt(figurepath+str(srcid)+'_2738.txt')
    if srcpsdinfo[1]<0:srcpsdinfo[1]=1e-3
    (src_evt_use, epoch_info_use) = rednoise.singleobs_psd.load_data(srcid,ecf=90,ifobsid=ifobsid,ifexpT=10000,path_provide=path_provide)
    w = make_freq_range(dt=lenbin, epoch_info=epoch_info_use)
    # w=w[w>1/50000.]
    lc_use = func.get_hist(src_evt_use[:,0], len_bin=100, tstart=epoch_info_use[:,0][0], tstop=epoch_info_use[:,1][-1])
    psd_org=Powerspectrum(lc_use)
    # (frms,frms_err)=psd_org.compute_rms(min_freq=w[0], max_freq=w[-1], white_noise_offset=0)
    # print('frms=',frms)
    # (frms,frms_err)=sr.excess_variance(lc_use,normalization='fvar')
    # print('frms=', frms)
    # frms=(lc_use.counts.std()/lc_use.counts.mean())
    # print('frms=', frms)
    # epoch_info_use=np.array(epoch_info_use)
    # psd_model = build_psd(x=w, p=srcpsdinfo[:-1], type='breakp')
    # psd_model=psd_model+srcpsdinfo[-1]  ## add poisson noise
    # psd_model=np.zeros(len(psd_model))+srcpsdinfo[-1]/1000  # pure poisson noise
    # print(psd_model)
    # CR=len(src_evt_use)/epoch_info_use[:,3]
    # lc = make_lc_from_psd(psd=psd_model, cts_rate=CR * lenbin, dt=lenbin, epoch_info=epoch_info_use)
    # print('counts={0}'.format(np.sum(lc.counts)))
    # psd=rednoise.plot_psd(lc,show=0)
    # plt.close()

    # rednoise.plot_psd(lc_evt,show=1)
    (src_evt_use, epoch_info_use) = rednoise.singleobs_psd.load_data(srcid,ecf=90,ifobsid=None,ifexpT=10000,path_provide=path_provide)
    CR = len(src_evt_use) / np.sum(epoch_info_use[:,3])
    # print(CR)
    evt_all = EventList()
    for i in range(len(epoch_info_use)):
        t_src = src_evt_use[:, 0][np.where(src_evt_use[:, 2] == epoch_info_use[i][2])]
        w = make_freq_range(dt=lenbin, epoch_info=epoch_info_use[i])
        psd_model = build_psd(x=w, p=srcpsdinfo[:-1], type='breakp')
        epoch_temp=np.array([epoch_info_use[i]])
        if len(t_src)<2 or epoch_info_use[i][3]<lenbin:
            # print(len(t_src),exptime[i])
            continue
        CR = len(t_src) / epoch_info_use[i][3]
        lc_new = make_lc_from_psd(psd=psd_model, cts_rate=0.5*CR * lenbin, dt=lenbin, epoch_info=epoch_temp,frms=0.1)
        lc_new.counts += 0.5*CR * lenbin
        lc_new.counts[np.where(lc_new.counts < 0)] = 0
        # plt.plot(lc_new.time,lc_new.counts)
        # plt.plot([lc_new.time[0],lc_new.time[-1]],[CR*lenbin,CR*lenbin],'--')
        # plt.show()
        # evt = EventList()
        # EventList.simulate_times(evt,lc=lc_new)
        # lc_evt=evt.to_lc(dt=lc.dt)
        # plt.plot(lc_evt.time,lc_evt.counts)
        # plt.show()
        # evt = EventList()
        # EventList.simulate_times(evt,lc=lc_new)
        lc_evt = Lightcurve(time=lc_new.time, counts=lc_new.counts)
        lc_evt.counts = np.random.poisson(lc_evt.counts)
        evt = EventList()
        # EventList.simulate_times(evt,lc=lc_evt)
        evt.time = func.sim_evtlist(lc_evt) + epoch_info_use[i][0]
        if i == 0:
            lc_all = lc_evt
        else:
            lc_all = lc_all.join(lc_evt)
        if i==0:
            evt_all=evt
        else:
            evt_all=evt_all.join(evt)

    return (evt_all,lc_all)

def simulate_const(src_evt,epoch_info,model='const'):
    ## src_evt [time,energy,obsid]
    tstart=epoch_info[:,0]
    tstop=epoch_info[:,1]
    obsid=epoch_info[:,2]
    exp_time_epoch=epoch_info[:,3]
    t_all=np.zeros(0)
    if model=='const':
        for i in range(len(exp_time_epoch)):
            N_cts=len(src_evt[np.where(src_evt[:,2]==obsid[i])])
            T=tstop[i]-tstart[i]
            N_cts_noise=np.random.poisson(N_cts)
            t=np.random.random(N_cts_noise)*T+tstart[i]
            t=np.sort(t)
            t_all=np.concatenate((t_all,t))

    return t_all


if __name__=='__main__':
    (evt,lc)=simulate_srcinfo(srcid=423)
    print(len(evt.time))
    print(np.sum(lc.counts))
    (src_evt_use, epoch_info_use) = rednoise.singleobs_psd.load_data(srcid=423, ecf=90, ifobsid=None, ifexpT=10000,
                                                                     path_provide=path_provide)
    hawk.phase_fold(time=evt.time,epoch_info=epoch_info_use,net_percent=net_p,p_test=18273.,outpath=figurepath,bin=20,shift=0.,
                    label=dataname,text='Seq.{}'.format(dataname),save=0,show=1)