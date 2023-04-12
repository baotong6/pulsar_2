#!/bin/bash
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 11 9:13:40 2023
@author: baotong
from real data to simulate lightcurve, eventlist
by best-fit psd or const
test everything
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
from rednoise import powerlw2evt as p2evt
from NGC104.timing_comb import load_data
from rednoise import Vaughan as Vaughan
import hawkeye as hawk

def make_freq_range(dt,epoch_info):
    if type(epoch_info)==np.ndarray:epoch_info=np.array(epoch_info)
    else:epoch_info=np.loadtxt(epoch_info)
    if epoch_info.ndim == 1:
        epoch_info=np.array([epoch_info])
    TSTART=epoch_info[:,0];TSTOP=epoch_info[:,1];OBSID=epoch_info[:,2];exptime=epoch_info[:,3]
    T_tot=TSTOP[-1]-TSTART[0]
    w = np.arange(2 / T_tot, 0.5 / dt, 1 / T_tot)
    return w

def bestpsd(lc,epoch_info):
    (result_mu,psd)=Vaughan.apply_Vaughan(lc,epoch_info=epoch_info,model=Vaughan.powerlaw)
    freq=psd.freq
    psd_sim=Vaughan.powerlaw(freq,result_mu)
    # plt.step(freq,psd_sim)
    # plt.loglog()
    # plt.show()
    return (psd_sim,result_mu)


def make_lc_from_psd(psd,cts_rate,dt,epoch_info,poisson=True,frms=0):
    if type(epoch_info)==np.ndarray:epoch_info=np.array(epoch_info)
    else:epoch_info=np.loadtxt(epoch_info)
    if epoch_info.ndim == 1:
        epoch_info=np.array([epoch_info])
    TSTART=epoch_info[:,0];TSTOP=epoch_info[:,1];OBSID=epoch_info[:,2];exptime=epoch_info[:,3]
    T_tot=TSTOP[-1]-TSTART[0]
    num_bins=int(T_tot/dt)
    sim = simulator.Simulator(N=num_bins + 1, mean=cts_rate,dt=dt,rms=frms,poisson=poisson)
    # lc=sim.simulate(1)
    lc=sim.simulate(psd)
    lc.time=lc.time+TSTART[0]
    # print('caution',np.where(lc.counts<0))
    lc.counts[np.where(lc.counts < 0)] = 0
    lc.gti=[[lc.time[0],lc.time[-1]]]
    return lc

def singleobs2simevt(path,dataname,dt=500,chosen_obs=None,ifoutobs=[],randseed=1):
    ## from one obs(usually the longest one obs) to simulate singleobs_evt
    ## if you need only just some obs to be simulated, make ifoutobs not []
    (src_evt_all, epoch_info_all) = load_data(dataname=dataname, ifpath=path, ifobsID=ifoutobs)
    if not chosen_obs:
        print('ID for chosen_obs must be given!')
        print('Use the longest obs')
        chosen_obs=[epoch_info_all[:,2][np.argmax(epoch_info_all[:,3])]]

    (src_evt_use, epoch_info_use) = load_data(dataname=dataname, ifpath=path, ifobsID=chosen_obs)
    time = src_evt_use[:, 0]
    dt = dt
    lc_cho = hawk.get_hist(time, len_bin=dt, tstart=epoch_info_use[:, 0][0], tstop=epoch_info_use[:, 1][-1])
    CR_cho = np.sum(lc_cho.counts) / np.sum(epoch_info_use[:, 3])
    frac_rms = np.sqrt(np.var(lc_cho.counts) * lc_cho.counts.size / (lc_cho.counts.size - 1) / np.mean(lc_cho.counts) ** 2)
    (psd_sim,result_mu) = bestpsd(lc_cho, epoch_info=epoch_info_use)

    evt_all = EventList()
    for i in range(len(epoch_info_all)):
        t_src = src_evt_all[:, 0][np.where(src_evt_all[:, 2] == epoch_info_all[i][2])]
        epoch_temp = np.array([epoch_info_all[i]])
        w = make_freq_range(dt=dt, epoch_info=epoch_temp)
        psd_model = Vaughan.powerlaw(x=w, p=result_mu)
        psd_model*=CR_cho**2 ## abs powerspectrum
        if len(t_src) < 2 or epoch_info_all[i][3] < dt:
            continue
        epoch_temp = np.array([epoch_info_all[i]])
        if len(t_src) < 2 or epoch_info_all[i][3] < dt:
            # print(len(t_src),exptime[i])
            continue
        CR = len(t_src) / epoch_info_all[i][3]
        lc_sim = make_lc_from_psd(psd_model, CR, dt, epoch_info=epoch_temp, frms=frac_rms, poisson=False)
        lc_sim.counts = lc_sim.counts * dt
        lc_sim.counts[np.where(lc_sim.counts < 0)] = 0
        np.random.seed(randseed)
        # lc_evt = Lightcurve(time=lc_sim.time, counts=lc_sim.counts)
        # lc_evt.counts = np.random.poisson(lc_evt.counts)
        evt = EventList()
        EventList.simulate_times(evt,lc=lc_sim)
        # evt.time = func.sim_evtlist(lc_sim) + epoch_info_use[i][0]
        lc_out = evt.to_lc(dt=dt)
        if i == 0:
            lc_all = lc_out
        else:
            lc_all = lc_all.join(lc_out)
        if i == 0:
            evt_all = evt
        else:
            evt_all = evt_all.join(evt)

    return (evt_all, lc_all)

def GL_simevt(simN=100):
    ##s模拟simN多少组
    path_M31 = '/Users/baotong/Desktop/M31XRB/M31ACIS_txt/txt_all_obs_p90/'  ##数据的路径
    dataname='2'  ##源的编号
    ifoutobs=[14197] #ifoutobs是，你要选的这个源，模拟他的哪些观测的数据。如果是所有的，想办法在别的地方读取一下，别手动敲这么多

    (src_evt_all, epoch_info_all) = load_data(dataname=dataname, ifpath=path_M31, ifobsID=ifoutobs)
    w_range=2*np.pi*np.arange(1./10000,1./3000,1.e-7)
    ## w_range针对你的信号来改，这个没关系，只要和我们table里这个源GL流程的range一致即可
    Prob,wpeak,wmean,mopt,wconf_lo,wconf_hi,simcounts = np.zeros(simN),np.zeros(simN),np.zeros(simN),\
                                                        np.zeros(simN),np.zeros(simN),np.zeros(simN),np.zeros(simN)
    for i in range(simN):
        (sim_evt_all, sim_lc_all) = singleobs2simevt(path=path_M31, dataname=dataname, dt=100, chosen_obs=None,
                                                     ifoutobs=ifoutobs,randseed=i)
        ## dt最好取大一点，以200为宜，再小的话rms会很大；如果是对短周期来做，就只能取小点。但是无所谓，因为短周期会被泊松噪声dominate
        time=sim_evt_all.time
        GL_R=hawk.GL_algorithm_single.compute_GL(time,epoch_info=epoch_info_all,w_range=w_range,m_max=20,parallel=True)
        Prob[i]=GL_R[1]
        wpeak[i]=GL_R[5]
        wmean[i]=GL_R[6]
        mopt[i]=GL_R[2]
        wconf_lo[i]=GL_R[7][0]
        wconf_hi[i]=GL_R[7][1]
        simcounts[i]=len(time)
    sim_result=np.column_stack((Prob,2*np.pi/wpeak,2*np.pi/wmean,mopt,
                                2*np.pi/wconf_lo-2*np.pi/wpeak,2*np.pi/wpeak-2*np.pi/wconf_hi,simcounts))
    np.savetxt('/Users/baotong/Desktop/M31XRB/M31ACIS_txt/rednoise/{0}.txt'.format(dataname),sim_result,
               fmt='%10.5f %10.2f %10.2f %10d %10.5f %10.5f %10d')
    ##改这里的path来存你的sim结果
if __name__=='__main__':
    GL_simevt(10)