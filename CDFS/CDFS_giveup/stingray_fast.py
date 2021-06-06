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
# font_prop = font_manager.FontProperties(size=16)
import warnings
from functools import reduce
import csv

font1 = {'family': 'Normal',
         'weight': 'normal',
         'size': 16, }

def filter_obs(src_evt,bkg_evt,useid):
    src_evt_use = src_evt[np.where(src_evt[:-1] == useid[0])[0]]
    bkg_evt_use = bkg_evt[np.where(bkg_evt[:-1] == useid[0])[0]]
    i=1
    while i < len(useid):
        id=useid[i]
        src_evt_use_temp=src_evt[np.where(src_evt[:-1]==id)[0]]
        bkg_evt_use_temp=bkg_evt[np.where(bkg_evt[:-1]==id)[0]]
        src_evt_use = np.concatenate((src_evt_use, src_evt_use_temp))
        bkg_evt_use = np.concatenate((bkg_evt_use, bkg_evt_use_temp))
        i+=1
    return (src_evt_use,bkg_evt_use)

def get_hist(t, len_bin):
    ###将输入的time信息，按照len_bin的长度输出为lc
    t_test = t-t[0]
    a = [0 for i in range(int(t_test[-1] / len_bin) + 1)]
    for i in range(len(t_test)):
        a[int(t_test[i] / len_bin)] += 1
    a = np.array(a)
    return a


def bending_po(x,p):
    """
    Parameters
    ----------
    :param x:numpy.ndarray
        non-zero frequencies
    :param p:
    p[0] = bending frequency
    p[1] = alpha,即 power law index、
    p[2] = constant
    p[3] = normalization N
    :return:
    """
    return p[3]*x**(-1)*(1+(x/p[0])**(p[1]-1))**(-1)+ p[2]

def generalized_lorentzian(x, p):
    """
    Generalized Lorentzian function.

    Parameters
    ----------
    x: numpy.ndarray
        non-zero frequencies

    p: iterable
        p[0] = peak centeral frequency
        p[1] = FWHM of the peak (gamma)
        p[2] = peak value at x=x0
        p[3] = power coefficient [n]
    Returns
    -------
    model: numpy.ndarray
        generalized lorentzian psd model
    """
    assert p[3] > 0., "The power coefficient should be greater than zero."
    return p[2] * (p[1] / 2)**p[3] * 1./(abs(x - p[0])**p[3] + (p[1] / 2)**p[3])
def smoothbknpo(x, p):
    """
    Smooth broken power law function.
    Parameters
    ----------

    x: numpy.ndarray
        non-zero frequencies

    p: iterable
        p[0] = normalization frequency
        p[1] = power law index for f --> zero
        p[2] = power law index for f --> infinity
        p[3] = break frequency
    Returns
    -------
    model: numpy.ndarray
        generalized smooth broken power law psd model
    """
    return p[0] * x**(-p[1]) / (1. + (x / p[3])**2)**(-(p[1] - p[2]) / 2)
def get_LS(time, flux,freq):
    x = time
    y = flux
    LS = LombScargle(x, y,normalization = 'standard')
    power = LS.power(freq)
    max_NormLSP=np.max(power)

    LS_power = LombScargle(x, y, normalization='psd')
    LSP=LS_power.power(freq)
    # power_freq=np.max(LSP[np.where((1/freq-2492.73)<200)])
    max_LSP=np.max(LSP)

    FP = LS.false_alarm_probability(power.max(), minimum_frequency=freq[0], maximum_frequency=freq[-1], method='baluev')

    plt.plot(freq, power)
    # plt.plot([1/2492,1/2492.],[0,np.max(power)],'--',linewidth=1)
    plt.semilogx()
    print(1. / freq[np.where(power == np.max(power))])
    out_period=1./freq[np.where(power==np.max(power))][0]
    # plt.show()
    plt.close()
    # return [FP,out_period]
    return [FP,out_period,max_LSP,max_NormLSP]

def sim_evtlist(lc):
    num_evt_bin=np.random.poisson(lc.counts)
    evt_all=[]
    i=0
    while i < len(num_evt_bin):
        evt_all=np.concatenate((evt_all,(np.random.uniform(lc.time[i]-lc.dt/2,lc.time[i]+lc.dt/2,num_evt_bin[i]))))
        i+=1
    return evt_all
def slow_sim_fixpds_qpo_XMM(src_index,num_trials=2):
    k_trial =0
    dt=100
    path='/Users/baotong/Desktop/CDFS/xmm_CDFS/xmm_txt/'
    epoch_file = np.loadtxt(path + 'epoch6_xmmobs.txt')

    tstart = epoch_file[:, 0];tstop = epoch_file[:, 1]
    ID=epoch_file[:,2];exptime = epoch_file[:, 3]

    # evt_file_origin =np.loadtxt(path + 'XID{0}_{1}_all_obs_10sec.txt'.format(source_id, det))
    # bkgevt_file_origin = np.loadtxt(path + 'bkg_XID{0}_{1}_all_obs_10sec.txt'.format(source_id, det))
    # (evt_file,bkgevt_file)=filter_obs(evt_file_origin,bkgevt_file_origin,ID)

    T_tot=tstop[-1]-tstart[0]
    # cts_rate = len(evt_file)/ (2 * np.sum(exptime)) * dt  # 实际的cts-rate应为这个的2倍
    cts_rate=0.001*0.67/2*dt ##这里cts_rate单位都是per bin
    bkg_cts_rate=0.001*0.33

    num_bins = int(T_tot / dt)
    freq=np.arange(1/T_tot,0.5/dt,1/(5*T_tot))
    freq=freq[np.where(freq > 1 / 20000.)]

    sim = simulator.Simulator(N=num_bins+1, mean=cts_rate, dt=dt)
    w = np.arange(1 / T_tot, 0.5 / dt, 1 / T_tot)
    spectrum = bending_po(w, [2.3e-3, 3.4, 0.40, 4.3e-4]) + generalized_lorentzian(w, [4.01e-4, 4.01e-4 / 16, 200, 2])
    # spectrum = bending_po(w, [2.3e-3, 3.4, 0.40, 4.3e-4])

    lc = sim.simulate(spectrum)
    lc.time=lc.time+tstart[0]
    # lc.gti=[[lc.time[0],lc.time[-1]]]
    lc.counts += cts_rate
    lc.counts[np.where(lc.counts < 0)] = 0
    with open(path + 'simulation/ep6_1e-3_SN3_LSP_qpo.txt', 'a+') as f:
        while k_trial < num_trials:
            ev_all = EventList();ev_all_bkg=EventList()
            for i in range(len(ID)):

                lc_cut = lc.truncate(start=tstart[i], stop=tstop[i],method='time')
                ev = EventList()
                ev.time = sim_evtlist(lc_cut)
                ev_all = ev_all.join(ev)
                index_b = (tstop[i]-tstart[i])*bkg_cts_rate
                temp = index_b
                cts_num = np.random.poisson(temp)
                ev_bkg = EventList()
                ev_bkg.time = np.random.uniform(tstart[i] - dt / 2, tstop[i] + dt / 2, cts_num)
                ev_all_bkg = ev_all_bkg.join(ev_bkg)

            lc_src = ev_all.to_lc(dt=dt, tstart=tstart[0] , tseg=T_tot)
            lc_bkg= ev_all_bkg.to_lc(dt=dt, tstart=tstart[0], tseg=T_tot)
            print(np.sum(lc_src))
            print(np.sum(lc_bkg))
            # T_exp=lc_new.time[-1]-lc_new.time[0]
            # lc_new=lc_src-lc_bkg
            lc_new=lc_src+lc_bkg

            temp = get_LS(lc_new.time, lc_new.counts, freq=freq)
            f.writelines((str(temp[0]) + '        ' + str(temp[1]) + '        ' + str(temp[2]) + '        ' + str(temp[3])+ '\n'))
            k_trial += 1
    f.close()

def slow_sim_fixpds_qpo_chandra(src_index,num_trials=2,dt=100):
    k_trial =0
    path='/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep3/'
    epoch_file=np.loadtxt(path+'epoch_src_{0}.txt'.format(src_index))
    # epoch_file=epoch_file[23:]

    tstart = epoch_file[:, 0];tstop = epoch_file[:, 1]
    ID=epoch_file[:,2];exptime = epoch_file[:, 3]
    evt_file_origin = np.loadtxt(path + '{0}.txt'.format(src_index))
    bkgevt_file_origin = np.loadtxt(path + '{0}_bkg.txt'.format(src_index))
    # (evt_file,bkgevt_file)=filter_obs(evt_file_origin,bkgevt_file_origin,ID)
    evt_file=evt_file_origin;bkgevt_file=bkgevt_file_origin
    T_tot=tstop[-1]-tstart[0]
    cts_rate = len(evt_file)/ (2 * np.sum(exptime)) * dt  # 实际的cts-rate应为这个的2倍
    bkg_rate=len(bkgevt_file_origin)/(2 * np.sum(exptime)) * dt
    num_bins = int(T_tot / dt)
    freq=np.arange(1/T_tot,0.5/dt,1/(5*T_tot))
    freq=freq[np.where(freq > 1 / 20000.)]

    sim = simulator.Simulator(N=num_bins + 1, mean=cts_rate, dt=dt)
    w = np.arange(1 / T_tot, 0.5 / dt, 1 / T_tot)
    spectrum = bending_po(w, [2.3e-3, 3.4, 0.40, 4.3e-4]) + generalized_lorentzian(w,
                                                                                   [4.01e-4, 4.01e-4 / 16, 200,
                                                                                    2])
    # spectrum = bending_po(w, [2.3e-3, 3.4, 0.40, 4.3e-4])
    lc = sim.simulate(spectrum)
    # lc=sim.simulate(pow)
    lc.time = lc.time + tstart[0]
    lc.gti=[[lc.time[0],lc.time[-1]]]
    lc.counts += cts_rate
    lc.counts[np.where(lc.counts < 0)] = 0
    ps = Powerspectrum(lc, norm='frac')
    fig, ax1 = plt.subplots(1, 1, figsize=(9, 6), sharex=True)
    ax1.loglog()
    ax1.step(ps.freq, ps.power, lw=2, color='blue')
    ax1.plot([1 / 950.7, 1 / 950.7], [0, np.max(ps.power)], '--', linewidth=1)
    ax1.set_ylabel("Frequency (Hz)", fontproperties=font1)
    ax1.set_ylabel("Power (raw)", fontproperties=font1)
    ax1.set_yscale('log')
    ax1.tick_params(axis='x', labelsize=16)
    ax1.tick_params(axis='y', labelsize=16)
    ax1.tick_params(which='major', width=1.5, length=7)
    ax1.tick_params(which='minor', width=1.5, length=4)
    for axis in ['top', 'bottom', 'left', 'right']:
        ax1.spines[axis].set_linewidth(1.5)
    plt.show()

    with open(path + 'simulation/{0}_LSP_FP_qpo_nofix.txt'.format(src_index), 'a+') as f:
        while k_trial < num_trials:


            ev_all = EventList();ev_all_bkg=EventList()
            for i in range(len(ID)):

                lc_cut = lc.truncate(start=tstart[i], stop=tstop[i],method='time')
                ev = EventList()
                ev.time = sim_evtlist(lc_cut)
                ev_all = ev_all.join(ev)
                index_b = len(np.where(bkgevt_file[:, 2] == ID[i])[0])
                temp = index_b/12.
                cts_num = np.random.poisson(temp)
                ev_bkg = EventList()
                ev_bkg.time = np.random.uniform(tstart[i] - dt / 2, tstop[i] + dt / 2, cts_num)
                ev_all_bkg = ev_all_bkg.join(ev_bkg)

            lc_src = ev_all.to_lc(dt=dt, tstart=tstart[0] , tseg=T_tot)
            lc_bkg= ev_all_bkg.to_lc(dt=dt, tstart=tstart[0], tseg=T_tot)
            print(np.sum(lc_src))
            print(np.sum(lc_bkg))
            # T_exp=lc_new.time[-1]-lc_new.time[0]
            lc_new=lc_src
            lc_new.counts-=lc_bkg.counts

            temp = get_LS(lc_new.time, lc_new.counts, freq=freq)
            # f.writelines((str(temp[0]) + '        ' + str(temp[1]) +'            '+ str(temp[2])+'             '+ str(temp[3])+'\n'))
            k_trial += 1
    f.close()

def slow_sim_fixpds_qpo_chandra_epoch23(src_index,num_trials=2,dt=100):
    k_trial =0
    def make_lc(evtfile,bkgevtfile,epochfile,dt):
        epoch_file=np.loadtxt(epochfile)
        evt_file=np.loadtxt(evtfile)
        bkgevt_file=np.loadtxt(bkgevtfile)

        tstart = epoch_file[:, 0];
        tstop = epoch_file[:, 1]
        ID = epoch_file[:, 2];
        exptime = epoch_file[:, 3]
        T_tot = tstop[-1] - tstart[0]
        cts_rate = len(evt_file) / (2 * np.sum(exptime)) * dt  # 实际的cts-rate应为这个的2倍
        num_bins = int(T_tot / dt)
        sim = simulator.Simulator(N=num_bins + 1, mean=cts_rate, dt=dt)
        w = np.arange(1 / T_tot, 0.5 / dt, 1 / (T_tot))
        spectrum = bending_po(w, [2.3e-3, 3.4, 0.40, 4.3e-4]) + generalized_lorentzian(w,
                                                                                       [4.01e-4, 4.01e-4 / 16, 200,
                                                                                        2])
        lc = sim.simulate(spectrum)
        # lc=sim.simulate(pow)
        lc.time = lc.time + tstart[0]
        # lc.gti=[[lc.time[0],lc.time[-1]]]
        lc.counts += cts_rate
        lc.counts[np.where(lc.counts < 0)] = 0

        return lc

    path2='/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep2/'
    epoch_file2=path2+'epoch_src_{0}.txt'.format(src_index)
    path3='/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep3/'
    epoch_file3 = path3 + 'epoch_src_{0}.txt'.format(src_index)
    lc2=make_lc(path2+'{0}.txt'.format(src_index),path2+'{0}_bkg.txt'.format(src_index),epoch_file2,dt)
    lc3= make_lc(path3+ '{0}.txt'.format(src_index), path3 + '{0}_bkg.txt'.format(src_index), epoch_file3, dt)
    lc=lc3.join(lc2)
    bkgevt_file=np.loadtxt(path2+'{0}_bkg_ep23.txt'.format(src_index))
    epoch_file=np.loadtxt(path2+'epoch_src_{0}_ep23.txt'.format(src_index))
    tstart = epoch_file[:, 0];
    tstop = epoch_file[:, 1]
    ID = epoch_file[:, 2];
    exptime = epoch_file[:, 3]
    T_tot = tstop[-1] - tstart[0]
    with open(path2 + 'simulation/{0}_LSP_FP_qpo_ep23.txt'.format(src_index), 'a+') as f:
        while k_trial < num_trials:

            ev_all = EventList();ev_all_bkg=EventList()
            for i in range(len(ID)):

                lc_cut = lc.truncate(start=tstart[i], stop=tstop[i],method='time')
                ev = EventList()
                ev.time = sim_evtlist(lc_cut)
                ev_all = ev_all.join(ev)

                index_b = len(np.where(bkgevt_file[:, 2] == ID[i])[0])
                temp = index_b/12.
                cts_num = np.random.poisson(temp)
                ev_bkg = EventList()
                ev_bkg.time = np.random.uniform(tstart[i] - dt / 2, tstop[i] + dt / 2, cts_num)
                ev_all_bkg = ev_all_bkg.join(ev_bkg)

            lc_src = ev_all.to_lc(dt=dt, tstart=tstart[0] , tseg=T_tot)
            lc_bkg= ev_all_bkg.to_lc(dt=dt, tstart=tstart[0], tseg=T_tot)

            # T_exp=lc_new.time[-1]-lc_new.time[0]
            lc_new_ep2_src = lc_src.truncate(start=tstart[0], stop=tstop[1],method='time')
            lc_new_ep2_bkg = lc_bkg.truncate(start=tstart[0], stop=tstop[1], method='time')
            lc_new_ep2_src.counts-=lc_new_ep2_bkg.counts

            lc_new_ep3_src = lc_src.truncate(start=tstart[2], stop=tstop[-1],method='time')
            lc_new_ep3_bkg = lc_bkg.truncate(start=tstart[2], stop=tstop[-1],method='time')
            lc_new_ep3_src.counts-=lc_new_ep3_bkg.counts

            print(np.sum(lc_new_ep2_src.counts))
            print(np.sum(lc_new_ep3_src.counts))

            freq2 = np.arange(1 / (tstop[1]-tstart[0]), 0.5 / dt, 1 / (5 * (tstop[1]-tstart[0])))
            freq2 = freq2[np.where(freq2 > 1 / 20000.)]

            freq3 = np.arange(1 / (tstop[-1]-tstart[2]), 0.5 / dt, 1 / (5 * (tstop[-1]-tstart[2])))
            freq3 = freq3[np.where(freq3 > 1 / 20000.)]

            temp2 = get_LS(lc_new_ep2_src.time, lc_new_ep2_src.counts, freq=freq2)
            temp3 = get_LS(lc_new_ep3_src.time, lc_new_ep3_src.counts, freq=freq3)
            f.writelines((str(temp2[0]) + '                 ' + str(temp2[1]) +'               '+ str(temp2[2])+'                 '+ str(temp2[3])+
            '               '+str(temp3[0]) + '                 ' + str(temp3[1]) +'               '+ str(temp3[2])+'                 '+ str(temp3[3])+'\n'))
            k_trial += 1
    f.close()

def plot_scatter_backnoback(src_index,period=2493.7655860349128,threshold=0.01):
    # path = '/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep3/'
    path='/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep2/'
    noqpo=np.loadtxt(path + 'simulation/19_LSP_FP_noqpo.txt')
    qpo=np.loadtxt(path + 'simulation/19_LSP_FP_qpo.txt')
    # const=np.loadtxt(path + 'simulation/{0}_LSP_FP_const.txt'.format(src_index))
    # plt.scatter(noback[:,1],noback[:,0],marker='o',s=10,color='green')
    plt.scatter(noqpo[:, 1], noqpo[:, 0], marker='x', s=30, color='red')
    plt.scatter(qpo[:, 1], qpo[:, 0], marker='o', s=30, color='green')
    # plt.scatter(const[:,1],const[:,0],marker='v',s=30,color='black')
    plt.xlabel('Period',font1)
    plt.ylabel('1-FAP',font1)
    plt.tick_params(labelsize=16)
    qpo_detect=qpo[np.where(qpo[:,0]<threshold)[0]]
    FWHM=32*period/(17*15)

    plt.plot([period-FWHM,period-FWHM],[1e-7,1],'--','grey')
    plt.plot([period+FWHM,period+FWHM],[1e-7,1],'--','grey')
    true_P=qpo_detect[np.where(np.abs(qpo_detect[:,1]-period)<FWHM)[0]]
    true_N=qpo_detect[np.where(np.abs(qpo_detect[:,1]-period)>FWHM)[0]]
    print('TP={0}'.format(len(true_P)))
    print('TN={0}'.format(len(true_N)))

    noqpo_detect=noqpo[np.where(noqpo[:,0]<threshold)[0]]
    false_P=noqpo_detect[np.where(np.abs(noqpo_detect[:,1]-period)<FWHM)[0]]
    false_N=noqpo_detect[np.where((np.abs(noqpo_detect[:,1]-period)>FWHM)&(noqpo_detect[:,1]>200.))[0]]

    print('FP={0}'.format(len(false_P)))
    print('FN={0}'.format(len(false_N)))


    # const_detect=const[np.where(const[:,0]<threshold)[0]]
    # false_P_const=const_detect[np.where(np.abs(const_detect[:,1]-period)<FWHM)[0]]
    # false_N_const=const_detect[np.where(np.abs(const_detect[:,1]-period)>FWHM)[0]]
    #
    # print('FP_const={0}'.format(len(false_P_const)))
    # print('FN_const={0}'.format(len(false_N_const)))
    # plt.ylim(1e-7,1)
    plt.semilogx()
    plt.semilogy()
    # plt.legend(['noback','back'])
    plt.plot([200,20000],[threshold,threshold],'--')
    plt.show()

def slow_sim_fixpds(k,src_index,num_trials=2):
    k_trial =0
    dt=100
    FP = [];
    period = [];
    cts_num=[];
    peakP=[];
    power_P=[]
    path = '/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep{0}/'.format(k)
    epoch_file = np.loadtxt(path + 'CDFS_epoch_ep{0}.txt'.format(k))
    tstart = epoch_file[:, 0];tstop = epoch_file[:, 1]
    ID=epoch_file[:,2];exptime = epoch_file[:, 3]
    evt_file = np.loadtxt(path + '{0}.txt'.format(src_index))
    bkgevt_file = np.loadtxt(path + '{0}_bkg.txt'.format(src_index))

    T_tot=tstop[-1]-tstart[0]
    cts_rate = len(evt_file)/ (2 * np.sum(exptime)) * dt  # 实际的cts-rate应为这个的2倍
    num_bins = int(T_tot / dt)

    sim = simulator.Simulator(N=num_bins+1, mean=cts_rate, dt=dt)
    w = np.arange(1 / T_tot, 0.5 / dt, 1 / T_tot)
    spectrum = bending_po(w, [2.3e-3, 3.4, 0.40, 4.3e-4]) + generalized_lorentzian(w, [2.7e-4, 2.7e-4 / 16, 200, 2])
    lc = sim.simulate(spectrum)
    lc.time=lc.time+tstart[0]
    # lc.gti=[[lc.time[0],lc.time[-1]]]
    lc.counts += cts_rate
    lc.counts[np.where(lc.counts < 0)] = 0

    ps = Powerspectrum(lc, norm='leahy')
    fig, ax1 = plt.subplots(1, 1, figsize=(9, 6), sharex=True)
    ax1.loglog()
    ax1.step(ps.freq, ps.power, lw=2, color='blue')
    ax1.plot([1 / 950.7, 1 / 950.7], [0, np.max(ps.power)], '--', linewidth=1)
    ax1.set_ylabel("Frequency (Hz)", fontproperties=font1)
    ax1.set_ylabel("Power (raw)", fontproperties=font1)
    ax1.set_yscale('log')
    ax1.tick_params(axis='x', labelsize=16)
    ax1.tick_params(axis='y', labelsize=16)
    ax1.tick_params(which='major', width=1.5, length=7)
    ax1.tick_params(which='minor', width=1.5, length=4)
    for axis in ['top', 'bottom', 'left', 'right']:
        ax1.spines[axis].set_linewidth(1.5)
    plt.show()
    T_exp = T_tot
    freq=np.arange(1/T_exp,0.5/dt,1/(5*T_exp))
    freq=freq[np.where(freq > 1 / 20000.)]
    with open (path+'/simulation/{0}_max_LSP_fixpds.txt'.format(src_index),'a+') as f:
        res_period=[];res_power=[]
        while k_trial<num_trials:
            ev_all = EventList()
            for i in range(len(exptime)):
                lc_cut = lc.truncate(start=tstart[i], stop=tstop[i],method='time')
                ev = EventList()
                ev.time = sim_evtlist(lc_cut)
                ev_all = ev_all.join(ev)

            lc_new = ev_all.to_lc(dt=dt, tstart=ev_all.time[0] - 0.5 * dt, tseg=ev_all.time[-1] - ev_all.time[0])

            temp = get_LS(lc_new.time, lc_new.counts, freq=freq, trial=k_trial)
            # res_period.append(temp[0]);res_power.append(temp[-1])
            f.writelines((str(temp[0])+'        '+str(temp[1])+'\n'))
            k_trial += 1
    f.close()
    # res=np.column_stack((res_power,res_period))
    # np.savetxt(path+'/simulation/{0}_max_LSP_fixpds.txt'.format(src_index),res,fmt='%10.2f %10.2f')
# slow_sim_fixpds('3','89',num_trials=1)

def slow_sim_const(k,src_index,num_trials=2):
    k_trial =0
    dt=100
    FP = [];
    period = [];
    cts_num=[];
    peakP=[];
    power_P=[]
    path = '/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep{0}/'.format(k)
    epoch_file = np.loadtxt(path + 'CDFS_epoch_ep{0}.txt'.format(k))
    tstart = epoch_file[:, 0];tstop = epoch_file[:, 1]
    ID=epoch_file[:,2];exptime = epoch_file[:, 3]
    evt_file = np.loadtxt(path + '{0}.txt'.format(src_index))
    bkgevt_file = np.loadtxt(path + '{0}_bkg.txt'.format(src_index))

    T_tot=tstop[-1]-tstart[0]
    cts_rate = len(evt_file)/ (2 * np.sum(exptime)) * dt  # 实际的cts-rate应为这个的2倍
    num_bins = int(T_tot / dt)
    T_exp = T_tot
    freq=np.arange(1/T_exp,0.5/dt,1/(5*T_exp))
    freq=freq[np.where(freq > 1 / 20000.)]

    with open(path + '/simulation/{0}_LSP_FP_const.txt'.format(src_index), 'a+') as f:
        while k_trial < num_trials:
            ev_all = EventList()
            for i in range(len(ID)):
                index = len(np.where(evt_file[:, 2] == ID[i])[0])
                index_b = len(np.where(bkgevt_file[:, 2] == ID[i])[0])
                temp = index - index_b / 12.
                if temp < 0: temp = 0
                cts_num = np.random.poisson(temp)
                ev = EventList()
                ev.time = np.random.uniform(tstart[i] - dt / 2, tstop[i] + dt / 2, cts_num)
                ev_all = ev_all.join(ev)
            lc_new = ev_all.to_lc(dt=dt, tstart=ev_all.time[0] - 0.5 * dt, tseg=ev_all.time[-1] - ev_all.time[0])
            # T_exp=lc_new.time[-1]-lc_new.time[0]
            print(np.sum(lc_new))
            temp = get_LS(lc_new.time, lc_new.counts, freq=freq)
            f.writelines((str(temp[0]) + '        ' + str(temp[1]) + '\n'))
            k_trial += 1
    f.close()
    # res=np.column_stack((res_power,res_period))
    # np.savetxt(path+'/simulation/{0}_max_LSP_fixpds.txt'.format(src_index),res,fmt='%10.2f %10.2f')
# slow_sim_const('3','19',num_trials=5000)
if __name__ == "__main__":
    # slow_sim_fixpds_qpo_chandra_epoch23('19', num_trials=10000, dt=100)
    # plot_scatter_backnoback('19', period=1 / 4.01e-4, threshold=0.13587)
    # slow_sim_fixpds_qpo_XMM('19', num_trials=1000)
    slow_sim_fixpds_qpo_chandra('19', num_trials=1)
