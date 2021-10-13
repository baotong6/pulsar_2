import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import pandas as pd
import sys
import os
from scipy import interpolate
from scipy.optimize import curve_fit
import astropy.units as u
import astropy.constants as c
from scipy import interpolate
import stingray as sr
from stingray.events import EventList
from stingray.lightcurve import Lightcurve
from stingray import Lightcurve, Crossspectrum, sampledata,Powerspectrum,AveragedPowerspectrum
from stingray.simulator import simulator, models
import matplotlib.font_manager as font_manager
from astropy.timeseries import LombScargle
import warnings
from functools import reduce
import csv
from CDFS.CDFS_startover import useful_functions as func
from CDFS.CDFS_startover import sim_psd as sim

font1=func.font1

def build_psd(x,p,x_2=1,p_2=1,type='bendp'):
    ##
    if type=='bendp':
        psd=func.bending_po(x,p)
    elif type=='lorentz':
        psd=func.generalized_lorentzian(x,p)
    elif type=='bendp+lorentz':
        psd=func.bending_po(x,p)+func.standard_lorentzian(x_2,p_2)
    return psd

def make_freq_range(dt,epoch_file):
    (TSTART, TSTOP, OBSID, exptime)=epoch_file
    T_tot=TSTOP[-1]-TSTART[0]
    w = np.arange(1 / T_tot, 0.5 / dt, 1 / T_tot)
    return w

def make_lc_from_psd(psd,cts_rate,dt,epoch_file):
    (TSTART, TSTOP, OBSID, exptime)=epoch_file
    T_tot=TSTOP[-1]-TSTART[0]
    num_bins=int(T_tot/dt)
    sim = simulator.Simulator(N=num_bins + 1, mean=cts_rate,dt=dt)
    lc=sim.simulate(psd)
    lc.time=lc.time+TSTART[0]
    # plt.plot(lc.time,lc.counts)
    # plt.show()
    # lc.counts += cts_rate
    lc.counts[np.where(lc.counts < 0)] = 0
    lc.gti=[[lc.time[0],lc.time[-1]]]
    return lc

def plot_psd(lc):
    ps = Powerspectrum(lc,norm='frac')
    fig, ax1 = plt.subplots(1, 1, figsize=(9, 6), sharex=True)
    ax1.loglog()
    ax1.step(ps.freq, ps.power, lw=2, color='blue')
    ax1.set_xlabel("Frequency (Hz)", fontproperties=font1)
    ax1.set_ylabel("Power ", fontproperties=font1)
    ax1.tick_params(axis='x', labelsize=16)
    ax1.tick_params(axis='y', labelsize=16)
    ax1.tick_params(which='major', width=1.5, length=7)
    ax1.tick_params(which='minor', width=1.5, length=4)
    for axis in ['top', 'bottom', 'left', 'right']:
        ax1.spines[axis].set_linewidth(1.5)
    plt.show()
    return ps

def test_somefunc():
    lenbin=100
    T_exp=59000
    CR=5e-4
    epoch_89='/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep3/epoch_src_test_single.txt'
    (TSTART, TSTOP, OBSID, exptime) = func.read_epoch(epoch_89)
    epoch_89=(TSTART, TSTOP, OBSID, exptime)
    w=make_freq_range(dt=lenbin,epoch_file=epoch_89)
    # psd_model=build_psd(x=w,p=[2.3e-3, 3.4, 0., 4.3e-4],x_2=w,
    #                     p_2=[4.01e-4, 4.01e-4 / 16, 100, 2],type='bendp+lorentz')
    psd_model=build_psd(x=w,p=[4.3e-4, 3.4, 0., 2e-2],x_2=w,
                        p_2=[2.67e-4, 16, 0.15, 2],type='bendp+lorentz')

    plt.loglog()
    plt.plot(w,psd_model)
    plt.show()
    lc=make_lc_from_psd(psd=psd_model,cts_rate=CR*lenbin,dt=lenbin,epoch_file=epoch_89)
    ps_org=plot_psd(lc)
    print('counts={0}'.format(np.sum(lc.counts)))
    lc_evt=Lightcurve(time=lc.time,counts=lc.counts,dt=lc.dt,gti=lc.gti)
    lc_evt.counts=np.random.poisson(lc_evt.counts)
    # ev_all = EventList()
    # ev_all.time = func.sim_evtlist(lc)
    # lc_evt = ev_all.to_lc(dt=lenbin, tstart=ev_all.time[0] - 0.5 * lenbin,
    #                       tseg=ev_all.time[-1] - ev_all.time[0])
    ps_real = plot_psd(lc_evt)
    freq = np.arange(1 / T_exp, 0.5 / lenbin, 1 / (5 * T_exp))
    freq = freq[np.where(freq > 1 / 10000.)]
    temp = func.get_LS(lc_evt.time, lc_evt.counts, freq=freq)

    plt.subplot(211)
    plt.plot(lc.time,lc.counts,color='red')
    plt.subplot(212)
    plt.plot(lc_evt.time,lc_evt.counts-lc.counts,color='green')
    print('counts={0}'.format(np.sum(lc_evt.counts)))
    plt.xlabel('Time')
    plt.ylabel('Counts/bin')
    plt.show()
    evt=EventList()
    EventList.simulate_times(evt,lc=lc)
    print('counts={0}'.format(len(evt.time)))

def sim_bunch_lc(CR,dt,period,epoch_file,outpath,num_trials=100):
    lenbin=dt
    # epoch_89 = '/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep3/epoch_src_89.txt'
    w = make_freq_range(dt=lenbin, epoch_file=epoch_file)
    qpo_f=1./period

    psd_model=build_psd(x=w,p=[4.3e-4, 3.4, 0., 2.3e-3],x_2=w,
                        p_2=[qpo_f, 16, 0.05, 2],type='bendp+lorentz')
    # psd_model = build_psd(x=w, p=[4.3e-4, 3.4, 0., 2.3e-3],
    #                       type='bendp')
    k_trial=0
    (TSTART, TSTOP, OBSID, exptime) = epoch_file
    # with open(outpath+'CR_{0}_P_{1}_A002_R15_fake4.txt'.format("%.0e" %CR ,"%.0d" %period),'a+') as f:
    # with open(outpath + 'CR_{0}_P_{1}_REJ1034_fake3.txt'.format("%.0e" % CR, "%.0d" % period), 'a+') as f:
    with open(outpath+'CR_{0}_P_{1}_REJ1034.txt'.format("%.0e" %CR ,"%.0d" %period),'a+') as f:
        while k_trial < num_trials:
            lc = make_lc_from_psd(psd=psd_model, cts_rate=CR* lenbin, dt=lenbin, epoch_file=epoch_file)
            print('trial'+':  '+str(k_trial))
            lc_evt = Lightcurve(time=lc.time, counts=lc.counts, dt=lc.dt, gti=lc.gti)
            lc_evt.counts = np.random.poisson(lc_evt.counts)
            T_tot=TSTOP[-1]-TSTART[0]
            for i in range(len(OBSID)):
                lc_cut = lc_evt.truncate(start=TSTART[i], stop=TSTOP[i], method='time')
                if i == 0:lc_long=lc_cut
                else:
                    lc_long=lc_long.join(lc_cut)
            print('counts={0}'.format(np.sum(lc_long.counts)))
            T_exp = lc_long.time[-1] - lc_long.time[0]
            freq = np.arange(1 / T_exp, 0.5 / lenbin, 1 / (5 * T_exp))
            freq = freq[np.where(freq > 1 / 20000.)]
            # print(T_tot)
            # print(T_exp)
            temp = func.get_LS(lc_long.time, lc_long.counts, freq=freq)
            f.writelines((str(temp[0]) + '        ' + str(temp[1]) + '        ' + str(temp[2]) + '  ' + '\n'))
            k_trial+=1
    f.close()

def sim_bunch_lc_evt(CR,dt,period,epoch_file,outpath,num_trials=100):
    lenbin=dt
    # epoch_89 = '/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep3/epoch_src_89.txt'
    w = make_freq_range(dt=lenbin, epoch_file=epoch_file)
    qpo_f=1./period
    (TSTART, TSTOP, OBSID, exptime) = epoch_file
    if type(TSTART)==type(1.2):
        print('Single exposure')
        TSTART=[TSTART];TSTOP=[TSTOP];OBSID=[OBSID];exptime=[exptime]

    with open(outpath+'CR_{0}_P_{1}_REJ1034_simpleway.txt'.format("%.0e"%CR,str(int(period))),'a+') as f:
        k_trial = 0
        while k_trial < num_trials:
            ev_all = EventList()
            # lc = make_lc_from_psd(psd=psd_model, cts_rate=CR / 2 * lenbin, dt=lenbin, epoch_file=epoch_file)
            print('trial'+':  '+str(k_trial))
            # lc_evt = Lightcurve(time=lc.time, counts=lc.counts, dt=lc.dt, gti=lc.gti)
            # lc_evt.counts = np.random.poisson(lc_evt.counts)
            T_tot=TSTOP[-1]-TSTART[0]
            for i in range(len(exptime)):
                cts_rate = 0.5 * CR *lenbin # 实际的cts-rate应为这个的2倍
                num_bins = int(exptime[i] / dt)
                sim = simulator.Simulator(N=num_bins, mean=cts_rate, dt=lenbin)
                w = np.arange(1 / exptime[i], 0.5 / lenbin, 1 / exptime[i])
                psd_model = build_psd(x=w, p=[2.3e-3, 3.4, 0., 4.3e-4], x_2=w, p_2=[qpo_f, qpo_f / 16, 100, 2],
                                      type='bendp+lorentz')
                # psd_model = build_psd(x=w, p=[2.3e-3, 3.4, 0., 4.3e-4],
                #                       type='bendp')
                spectrum=psd_model
                lc_cut = sim.simulate(spectrum)
                lc_cut.counts += cts_rate
                lc_cut.counts[np.where(lc_cut.counts < 0)] = 0
                lc_cut.counts=np.random.poisson(lc_cut.counts)
                # print(np.sum(lc_cut.counts))
                ev = EventList()
                ev.time = func.sim_evtlist(lc_cut) + TSTART[i]
                ev_all = ev_all.join(ev)
            lc_long = ev_all.to_lc(dt=dt, tstart=ev_all.time[0] - 0.5 * dt, tseg=ev_all.time[-1] - ev_all.time[0])
            print('counts={0}'.format(np.sum(lc_long.counts)))
            T_exp = lc_long.time[-1] - lc_long.time[0]
            freq = np.arange(1 / T_exp, 0.5 / lenbin, 1 / (5 * T_exp))
            freq = freq[np.where(freq > 1 / 20000.)]
            # print(T_tot)
            # print(T_exp)
            temp = func.get_LS(lc_long.time, lc_long.counts, freq=freq)
            f.writelines((str(temp[0]) + '        ' + str(temp[1]) + '        ' + str(temp[2]) + '  ' + '\n'))
            k_trial+=1
    f.close()

def run_many_sim():
    epoch1 = '/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep1/CDFS_epoch_ep1.txt'
    epoch2 = '/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep2/CDFS_epoch_ep2.txt'
    epoch3 = '/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep3/CDFS_epoch_ep3.txt'
    epoch4 = '/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep4/CDFS_epoch_ep4.txt'
    epoch_all=[epoch1,epoch2,epoch3,epoch4]
    CR_all=np.array([5e-4,6e-4,7e-4,8e-4,9e-4,1e-3])
    # period_all=[1800,3600,7200]
    period_all=[3600,5400,7200]

    ep=2;i=5;j=2

    (TSTART, TSTOP, OBSID, exptime) = func.read_epoch(epoch_all[ep])
    # id=5  ##前id次观测
    # TSTART=TSTART[0:id];TSTOP=TSTOP[0:id];OBSID=OBSID[0:id];exptime=exptime[0:id]

    sim_bunch_lc(CR=CR_all[i], dt=100, period=period_all[j],
                 epoch_file=(TSTART, TSTOP, OBSID, exptime),
                 outpath='/Users/baotong/Desktop/CDFS/simulation/EP{0}/'.format(int(ep + 1)),
                 num_trials=1000)


    # ep=1;i=5;j=1
    #
    # (TSTART, TSTOP, OBSID, exptime) = func.read_epoch(epoch_all[ep])
    # # id=5  ##前id次观测
    # # TSTART=TSTART[0:id];TSTOP=TSTOP[0:id];OBSID=OBSID[0:id];exptime=exptime[0:id]
    #
    # sim_bunch_lc(CR=CR_all[i], dt=100, period=period_all[j],
    #              epoch_file=(TSTART, TSTOP, OBSID, exptime),
    #              outpath='/Users/baotong/Desktop/CDFS/simulation/EP{0}/'.format(int(ep + 1)),
    #              num_trials=1000)
    # ep=1;i=5;j=2
    #
    # (TSTART, TSTOP, OBSID, exptime) = func.read_epoch(epoch_all[ep])
    # # id=5  ##前id次观测
    # # TSTART=TSTART[0:id];TSTOP=TSTOP[0:id];OBSID=OBSID[0:id];exptime=exptime[0:id]
    #
    # sim_bunch_lc(CR=CR_all[i], dt=100, period=period_all[j],
    #              epoch_file=(TSTART, TSTOP, OBSID, exptime),
    #              outpath='/Users/baotong/Desktop/CDFS/simulation/EP{0}/'.format(int(ep + 1)),
    #              num_trials=1000)

    ##内存有限，别用循环一次性做啦##
    # for ep in range(len(epoch_all)):
    #     (TSTART, TSTOP, OBSID, exptime) = func.read_epoch(epoch_all[ep])
    # #     sim_bunch_lc(CR=CR_all[i], dt=100, period=period_all[j],
    # #                  epoch_file=(TSTART, TSTOP, OBSID, exptime),
    # #                  outpath='/Users/baotong/Desktop/CDFS/simulation/EP{0}/'.format(int(ep + 1)),
    # #                  num_trials=1000)
    #     for i in range(len(CR_all)):
    #         for j in range(len(period_all)):
    #             sim_bunch_lc(CR=CR_all[i],dt=100,period=period_all[j],
    #                              epoch_file=(TSTART, TSTOP, OBSID, exptime),
    #                              outpath='/Users/baotong/Desktop/CDFS/simulation/EP{0}/'.format(int(ep+1)),
    #                              num_trials=1000)
    ## ---------------------##
if __name__=='__main__':
    # run_many_sim()
    # test_somefunc()
    run_many_sim()