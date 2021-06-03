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

warnings.filterwarnings('ignore')
font1 = {'family': 'Normal',
         'weight': 'normal',
         'size': 16, }
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
def get_LS(time, flux,freq,trial=1):
    path='/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep3/simulation/19_LS_sim_noQPO/'
    x = time
    y = flux
    # dy=np.sqrt(y)
    # plt.scatter(x,y)
    # plt.show()

    # LS = LombScargle(x, y, dy = 1, normalization = 'standard', fit_mean = True,
    #                  center_data = True).power(freq, method = 'cython')
    LS = LombScargle(x, y,normalization = 'standard')
    # LS = LombScargle(x, y, dy, normalization='psd')
    power = LS.power(freq)

    # print('freq_num={0}'.format(len(freq)))
    FP=LS.false_alarm_probability(power.max(),minimum_frequency = freq[0], maximum_frequency = freq[-1],method='baluev')
    # FP_99 = LS.false_alarm_level(0.0027, minimum_frequency = freq[0], maximum_frequency = freq[-1],method='baluev')
    # FP_90 = LS.false_alarm_level(0.05,  minimum_frequency=freq[0],
    #                              maximum_frequency=freq[-1], method='baluev')
    # FP_68 = LS.false_alarm_level(0.32, minimum_frequency=freq[0],
    #                              maximum_frequency=freq[-1], method='baluev')
    # if FP<0.01:print(dataname)
    # plt.title('FP={0}'.format(FP))
    # plt.semilogx()
    # plt.plot(freq, power)
    # print(1./freq[np.where(power==np.max(power))])
    # plt.plot([freq[0], freq[-1]], [FP_99, FP_99], '--'
    # plt.plot([freq[0], freq[-1]], [FP_90, FP_90], '--')
    # plt.plot([freq[0], freq[-1]], [FP_68, FP_68], '--')
    # plt.show()
    res=1e2*power
    res=np.round(res,2)
    # res=res[::smooth]
    # np.savetxt(path+'LS_simres_{0}.txt'.format(trial+1),res,fmt='%10.5f')

    return [1. / freq[np.where(power == np.max(power))],np.max(power),res]

def sim_evtlist(lc):
    num_evt_bin=np.random.poisson(lc.counts)
    evt_all=[]
    i=0
    while i < len(num_evt_bin):
        evt_all=np.concatenate((evt_all,(np.random.uniform(lc.time[i]-lc.dt/2,lc.time[i]+lc.dt/2,num_evt_bin[i]))))
        i+=1
    return evt_all
def get_lc_byspec_single():
    ## stingray做simulation的时候powerspectrum的normlization的对lc似乎毫无影响；反正模拟的结果是这么显示的``````
    exptime=1e6
    dt = 10
    cts_rate = 1e-4* dt  # 实际得到的lc中的cts-rate应为这个的2倍，换言之这里应该是你要的cts-rate除以2
    num_bins = int(exptime / dt)
    sim = simulator.Simulator(N=num_bins, mean=cts_rate, dt=dt)
    w = np.arange(1 / exptime, 0.5 / dt, 1 / exptime)
    # w = np.fft.rfftfreq(sim.N, d=sim.dt)[1:]
    spectrum=generalized_lorentzian(w,[1/500.,1/4000.,100,2])
    # spectrum = smoothbknpo(w, [cts_rate * dt * 0.1, 3, 0, 1e-3])
    # spectrum = bending_po(w, [2.3e-3, 3.4, 0.40, 4.3e-4]) ## 这是RE J1034+396 的参数
    spectrum= bending_po(w, [2.3e-3, 3.4, 0.40, 4.3e-4])+generalized_lorentzian(w,[2.7e-4,2.7e-4/16,200,2])
    # spectrum = smoothbknpo(w, [cts_rate*dt*1, 2, 0, 1e-3])  + generalized_lorentzian(w, [1 / 1000., 1 / 10000., cts_rate*dt*100, 2])
    print(np.mean(spectrum))
    lc = sim.simulate(spectrum)
    lc.counts += cts_rate
    lc.counts[np.where(lc.counts < 0)] = 0
    print(np.mean(lc.counts))
    print(np.sum(lc.counts))
    print(np.var(lc.counts))
    plt.plot(lc.time,lc.counts)
    plt.show()
    ps = Powerspectrum(lc, norm='abs')
    fig, ax1 = plt.subplots(1, 1, figsize=(9, 6), sharex=True)
    ax1.plot(ps.freq, ps.power, lw=2, color='blue')
    ax1.set_ylabel("Frequency (Hz)", fontproperties=font_prop)
    ax1.set_ylabel("Power (raw)", fontproperties=font_prop)
    # ax1.set_yscale('log')
    ax1.tick_params(axis='x', labelsize=16)
    ax1.tick_params(axis='y', labelsize=16)
    ax1.tick_params(which='major', width=1.5, length=7)
    ax1.tick_params(which='minor', width=1.5, length=4)
    for axis in ['top', 'bottom', 'left', 'right']:
        ax1.spines[axis].set_linewidth(1.5)
    plt.loglog()
    plt.show()

    T_exp =  exptime
    freq = np.arange(1 / T_exp, 0.5 / dt, 1 / T_exp)
    freq = freq[np.where(freq > 1 / 20000.)]
    # temp = get_LS(lc.time, lc.counts, freq=freq)
# get_lc_byspec_single()

cr=[2e-5,5e-5,1e-4,2e-4,3e-4,4e-4,5e-4]
cr_str=['2e-5','5e-5','1e-4','2e-4','3e-4','4e-4','5e-4']
def get_lc_byspec(k,j,num_trials=1000):
    k_trial=0
    FP = [];
    period = [];
    while k_trial <num_trials:
        path = '/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep{0}/'.format(k)
        epoch_file=np.loadtxt(path+'CDFS_epoch_ep{0}.txt'.format(k))
        # path = '/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8/'
        # epoch_file = np.loadtxt(path + 'CDFS_epoch.txt')
        tstart=epoch_file[:,0];tstop=epoch_file[:,1];exptime=epoch_file[:,3]
        ev_all = EventList()
        for i in range(len(exptime)):
            dt=100
            cts_rate=cr[j] *dt  #实际的cts-rate应为这个的2倍
            num_bins=int(exptime[i]/dt)
            sim = simulator.Simulator(N=num_bins, mean=cts_rate, dt=dt)
            w = np.arange(1 / exptime[i], 0.5 / dt, 1 / exptime[i])
            # w = np.fft.rfftfreq(sim.N, d=sim.dt)[1:]
            # spectrum = smoothbknpo(w, [0.01, 2, 1e-2, 1e-3])
            spectrum = bending_po(w, [2.3e-3, 3.4, 0.40, 4.3e-4]) ## 这是RE J1034+396 的参数
            # spectrum = bending_po(w, [2.3e-3, 3.4, 0.40, 4.3e-4]) + generalized_lorentzian(w, [6.68e-4 ,6.68e-4 /16,200,2])
            # spectrum = bending_po(w, [2.3e-3, 3.4, 0.40, 4.3e-4]) + generalized_lorentzian(w,[2.7e-4, 2.7e-4 / 16, 200,2])
            # spectrum =powerlaw + generalized_lorentzian(w, [1 / 1000., 1 / 10000., 0.5*np.max(powerlaw), 2])
            lc = sim.simulate(spectrum)
            lc.counts += cts_rate
            lc.counts[np.where(lc.counts<0)]=0
            # plt.plot(w,spectrum)
            # plt.loglog()
            # plt.xlabel("Frequency (Hz)", fontproperties=font_prop)
            # plt.ylabel("Power (abs)", fontproperties=font_prop)
            # plt.show()
            ps=Powerspectrum(lc,norm='abs')
            # fig, ax1 = plt.subplots(1, 1, figsize=(9, 6), sharex=True)
            # ax1.plot(ps.freq, ps.power, lw=2, color='blue')
            # ax1.set_ylabel("Frequency (Hz)", fontproperties=font_prop)
            # ax1.set_ylabel("Power (raw)", fontproperties=font_prop)
            # # ax1.set_yscale('log')
            # ax1.tick_params(axis='x', labelsize=16)
            # ax1.tick_params(axis='y', labelsize=16)
            # ax1.tick_params(which='major', width=1.5, length=7)
            # ax1.tick_params(which='minor', width=1.5, length=4)
            # for axis in ['top', 'bottom', 'left', 'right']:
            #     ax1.spines[axis].set_linewidth(1.5)
            # plt.loglog()
            # plt.show()
            ev = EventList()
            # ev.simulate_times(use_spline=False,lc=lc,bin_time=dt)
            # ev.time+=tstart[i]
            ev.time=sim_evtlist(lc)+tstart[i]
            # print(ev.time)RX J1301.9+2747
            ev_all=ev_all.join(ev)
            # lc_temp = ev.to_lc(dt=dt, tstart=ev.time[0]-0.5*dt, tseg=ev.time[-1]-ev.time[0])
            # ps_temp=Powerspectrum(lc_temp,norm='abs')
            # fig2, ax2 = plt.subplots(1, 1, figsize=(9, 6), sharex=True)
            # ax2.plot(ps_temp.freq, ps_temp.power, lw=2, color='blue')
            # ax2.set_ylabel("Frequency (Hz)", fontproperties=font_prop)
            # ax2.set_ylabel("Power (raw)", fontproperties=font_prop)
            # # ax1.set_yscale('log')
            # ax2.tick_params(axis='x', labelsize=16)
            # ax2.tick_params(axis='y', labelsize=16)
            # ax2.tick_params(which='major', width=1.5, length=7)
            # ax2.tick_params(which='minor', width=1.5, length=4)
            # for axis in ['top', 'bottom', 'left', 'right']:
            #     ax2.spines[axis].set_linewidth(1.5)
            # plt.loglog()
            # plt.show()
            # 
            # plt.figure(2)
            # plt.plot(lc.time,lc.counts,c='green')
            # plt.plot(lc_temp.time,lc_temp.counts,c='r')
            # plt.show()
            # plt.figure(3)
            # get_LS(lc.time,lc.counts,w)
        # print('cts={0}'.format(len(ev_all.time)))
        lc_new = ev_all.to_lc(dt=dt, tstart=ev_all.time[0]-0.5*dt, tseg=ev_all.time[-1]-ev_all.time[0])
        # plt.figure(2)
        # plt.plot(lc_new.time,lc_new.counts,c='r')
        # plt.show()
        # plt.figure(3)
        T_exp=lc_new.time[-1]-lc_new.time[0]
        freq=np.arange(1/T_exp,0.5/dt,1/(5*T_exp))
        freq=freq[np.where(freq > 1 / 10000.)]
        # print(len(freq))
        temp=get_LS(lc_new.time, lc_new.counts, freq=freq)
        FP.append(temp[0]);period.append(temp[1])
        k_trial+=1
    result=np.column_stack((FP,period))
    np.savetxt(path+'simulation/'+'trial_out_{0}_REJ1034+396_noQPO.txt'.format(cr_str[j]),result,fmt="%10.5f %10.5f")
    return ev_all
# get_lc_byspec_single()
def get_lc_byspec_1hr(k,j,num_trials=1000):
    k_trial=0
    FP = [];
    period = [];
    while k_trial <num_trials:
        path = '/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep{0}/'.format(k)
        epoch_file=np.loadtxt(path+'CDFS_epoch_ep{0}.txt'.format(k))
        # path = '/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8/'
        # epoch_file = np.loadtxt(path + 'CDFS_epoch.txt')
        tstart=epoch_file[:,0];tstop=epoch_file[:,1];exptime=epoch_file[:,3]
        ev_all = EventList()
        for i in range(len(exptime)):
            dt=100
            cts_rate=cr[j] *dt  #实际的cts-rate应为这个的2倍
            num_bins=int(exptime[i]/dt)
            sim = simulator.Simulator(N=num_bins, mean=cts_rate, dt=dt)
            w = np.arange(1 / exptime[i], 0.5 / dt, 1 / exptime[i])
            # w = np.fft.rfftfreq(sim.N, d=sim.dt)[1:]
            # spectrum = smoothbknpo(w, [0.01, 2, 1e-2, 1e-3])
            # spectrum = bending_po(w, [2.3e-3, 3.4, 0.40, 4.3e-4]) ## 这是RE J1034+396 的参数
            # spectrum = bending_po(w, [2.3e-3, 3.4, 0.40, 4.3e-4]) + generalized_lorentzian(w, [6.68e-4 ,6.68e-4 /16,200,2])
            spectrum = bending_po(w, [2.3e-3, 3.4, 0.40, 4.3e-4]) + generalized_lorentzian(w,[2.7e-4, 2.7e-4 / 16, 200,2])
            # spectrum =powerlaw + generalized_lorentzian(w, [1 / 1000., 1 / 10000., 0.5*np.max(powerlaw), 2])
            lc = sim.simulate(spectrum)
            lc.counts += cts_rate
            lc.counts[np.where(lc.counts<0)]=0
            # plt.plot(w,spectrum)
            # plt.loglog()
            # plt.xlabel("Frequency (Hz)", fontproperties=font_prop)
            # plt.ylabel("Power (abs)", fontproperties=font_prop)
            # plt.show()
            ps=Powerspectrum(lc,norm='abs')
            # fig, ax1 = plt.subplots(1, 1, figsize=(9, 6), sharex=True)
            # ax1.plot(ps.freq, ps.power, lw=2, color='blue')
            # ax1.set_ylabel("Frequency (Hz)", fontproperties=font_prop)
            # ax1.set_ylabel("Power (raw)", fontproperties=font_prop)
            # # ax1.set_yscale('log')
            # ax1.tick_params(axis='x', labelsize=16)
            # ax1.tick_params(axis='y', labelsize=16)
            # ax1.tick_params(which='major', width=1.5, length=7)
            # ax1.tick_params(which='minor', width=1.5, length=4)
            # for axis in ['top', 'bottom', 'left', 'right']:
            #     ax1.spines[axis].set_linewidth(1.5)
            # plt.loglog()
            # plt.show()
            ev = EventList()
            # ev.simulate_times(use_spline=False,lc=lc,bin_time=dt)
            # ev.time+=tstart[i]
            ev.time=sim_evtlist(lc)+tstart[i]
            # print(ev.time)RX J1301.9+2747
            ev_all=ev_all.join(ev)
            # lc_temp = ev.to_lc(dt=dt, tstart=ev.time[0]-0.5*dt, tseg=ev.time[-1]-ev.time[0])
            # ps_temp=Powerspectrum(lc_temp,norm='abs')
            # fig2, ax2 = plt.subplots(1, 1, figsize=(9, 6), sharex=True)
            # ax2.plot(ps_temp.freq, ps_temp.power, lw=2, color='blue')
            # ax2.set_ylabel("Frequency (Hz)", fontproperties=font_prop)
            # ax2.set_ylabel("Power (raw)", fontproperties=font_prop)
            # # ax1.set_yscale('log')
            # ax2.tick_params(axis='x', labelsize=16)
            # ax2.tick_params(axis='y', labelsize=16)
            # ax2.tick_params(which='major', width=1.5, length=7)
            # ax2.tick_params(which='minor', width=1.5, length=4)
            # for axis in ['top', 'bottom', 'left', 'right']:
            #     ax2.spines[axis].set_linewidth(1.5)
            # plt.loglog()
            # plt.show()
            #
            # plt.figure(2)
            # plt.plot(lc.time,lc.counts,c='green')
            # plt.plot(lc_temp.time,lc_temp.counts,c='r')
            # plt.show()
            # plt.figure(3)
            # get_LS(lc.time,lc.counts,w)
        # print('cts={0}'.format(len(ev_all.time)))
        lc_new = ev_all.to_lc(dt=dt, tstart=ev_all.time[0]-0.5*dt, tseg=ev_all.time[-1]-ev_all.time[0])
        # plt.figure(2)
        # plt.plot(lc_new.time,lc_new.counts,c='r')
        # plt.show()
        # plt.figure(3)
        T_exp=lc_new.time[-1]-lc_new.time[0]
        freq = np.arange(1 / T_exp, 0.5 / dt, 1 / (5 * T_exp))
        freq=freq[np.where(freq > 1 / 20000.)]
        # print(len(freq))
        temp=get_LS(lc_new.time, lc_new.counts, freq=freq)
        FP.append(temp[0]);period.append(temp[1])
        k_trial+=1
    result=np.column_stack((FP,period))
    np.savetxt(path+'simulation/'+'trial_out_1hr_{0}_REJ1034+396.txt'.format(cr_str[j]),result,fmt="%10.5f %10.5f")
    return ev_all
def get_lc_onesource(k,src_index,num_trials=2):
    k_trial =0
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
    for i in range(len(ID)):
        index = len(np.where(evt_file[:,2] == ID[i])[0])
        index_b=len(np.where(bkgevt_file[:,2] == ID[i])[0])
        cts_num.append(index-index_b/12.)
    dt = 100
    T_exp = 11000154.981141508
    freq=np.arange(1/T_exp,0.5/dt,1/(5*T_exp))
    freq=freq[np.where(freq > 1 / 20000.)]
    if os.path.exists(path+'/simulation/{0}_LS_simP.csv'.format(src_index)):
        print('caution! file exists')
        return None
    with open(path + '/simulation/{0}_LS_simP.csv'.format(src_index), 'a+') as csvfile:
        header = freq
        header = header.astype('str')
        writer = csv.writer(csvfile)
        while k_trial <num_trials:
            ev_all = EventList()
            for i in range(len(exptime)):
                cts_rate = cts_num[i]/(2*exptime[i]) * dt  # 实际的cts-rate应为这个的2倍
                num_bins = int(exptime[i] / dt)
                sim = simulator.Simulator(N=num_bins, mean=cts_rate, dt=dt)
                w = np.arange(1 / exptime[i], 0.5 / dt, 1 / exptime[i])
                spectrum = bending_po(w, [2.3e-3, 3.4, 0.40, 4.3e-4])
                # spectrum = bending_po(w, [2.3e-3, 3.4, 0.40, 4.3e-4]) + generalized_lorentzian(w, [1.0518215e-3,1.0518215e-3/16,200,2])
                lc = sim.simulate(spectrum)
                lc.counts += cts_rate
                lc.counts[np.where(lc.counts<0)]=0
                ps = Powerspectrum(lc, norm='abs')
                ev = EventList()
                ev.time = sim_evtlist(lc) + tstart[i]
                ev_all = ev_all.join(ev)
            # print(len(ev_all.time))
            lc_new = ev_all.to_lc(dt=dt, tstart=ev_all.time[0]-0.5*dt, tseg=ev_all.time[-1]-ev_all.time[0])
            # T_exp=lc_new.time[-1]-lc_new.time[0]
            temp=get_LS(lc_new.time, lc_new.counts, freq=freq,trial=k_trial)
            writer.writerows([temp[-1]])
            k_trial+=1

def get_lc_onesource_fixpds(k,src_index,num_trials=2):
    k_trial =0
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
    for i in range(len(ID)):
        index = len(np.where(evt_file[:,2] == ID[i])[0])
        index_b=len(np.where(bkgevt_file[:,2] == ID[i])[0])
        cts_num.append(index-index_b/12.)
    dt = 100
    T_exp = 11000154.981141508
    freq=np.arange(1/T_exp,0.5/dt,1/(5*T_exp))
    freq=freq[np.where(freq > 1 / 20000.)]
    if os.path.exists(path+'/simulation/{0}_LS_simP_fixpds.csv'.format(src_index)):
        print('caution! file exists')
        return None
    with open(path + '/simulation/{0}_LS_simP_fixpds.csv'.format(src_index), 'a+') as csvfile:
        header = freq
        header = header.astype('str')
        writer = csv.writer(csvfile)
        for i in range(len(exptime)):
            cts_rate = cts_num[i] / (2 * exptime[i]) * dt  # 实际的cts-rate应为这个的2倍
            num_bins = int(exptime[i] / dt)
            sim = simulator.Simulator(N=num_bins, mean=cts_rate, dt=dt)
            w = np.arange(1 / exptime[i], 0.5 / dt, 1 / exptime[i])
            spectrum = bending_po(w, [2.3e-3, 3.4, 0.40, 4.3e-4])
            # spectrum = bending_po(w, [2.3e-3, 3.4, 0.40, 4.3e-4]) + generalized_lorentzian(w, [1.0518215e-3,1.0518215e-3/16,200,2])
            lc = sim.simulate(spectrum)
            lc.counts += cts_rate
            lc.counts[np.where(lc.counts < 0)] = 0
            lc.time+=tstart[i]
            if i==0: lc_all=lc
            else: lc_all=lc_all.join(lc)
        print('run')
        print(lc_all.time)
        while k_trial<num_trials:
            ev_all = EventList()
            ev_all.time = sim_evtlist(lc_all) + tstart[0]
            # ev_all = ev_all.join(ev)
            lc_new = ev_all.to_lc(dt=dt, tstart=ev_all.time[0] - 0.5 * dt, tseg=ev_all.time[-1] - ev_all.time[0])
            temp = get_LS(lc_new.time, lc_new.counts, freq=freq, trial=k_trial)
            writer.writerows([temp[-1]])
            k_trial+=1
# get_lc_onesource_fixpds(3,'236',num_trials=10000)
def get_lc_onesource_const(k,src_index,num_trials=100):
    dt=100
    k_trial =0
    FP = [];
    period = [];
    cts_num=[];
    peakP=[];
    power_P=[]
    dt = 100
    T_exp = 11000154.981141508
    freq=np.arange(1/T_exp,0.5/dt,1/(5*T_exp))
    freq=freq[np.where(freq > 1 / 20000.)]
    path = '/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep{0}/'.format(k)
    epoch_file = np.loadtxt(path + 'CDFS_epoch_ep{0}.txt'.format(k))
    tstart = epoch_file[:, 0];tstop = epoch_file[:, 1]
    ID=epoch_file[:,2];exptime = epoch_file[:, 3]
    evt_file = np.loadtxt(path + '{0}.txt'.format(src_index))
    bkgevt_file = np.loadtxt(path + '{0}_bkg.txt'.format(src_index))
    if os.path.exists(path+'/simulation/{0}_LS_simP_const.csv'.format(src_index)):
        print('caution! file exists')
        return None
    with open(path + '/simulation/{0}_LS_simP_const.csv'.format(src_index), 'a+') as csvfile:
        header = freq
        header = header.astype('str')
        writer = csv.writer(csvfile)
        while k_trial <num_trials:
            ev_all = EventList()
            for i in range(len(ID)):
                index = len(np.where(evt_file[:,2] == ID[i])[0])
                index_b=len(np.where(bkgevt_file[:,2] == ID[i])[0])
                temp=index-index_b/12.
                if temp<0:temp=0
                cts_num=np.random.poisson(temp)
                ev=EventList()
                ev.time=np.random.uniform(tstart[i]-dt/2,tstop[i]+dt/2,cts_num)
                ev_all = ev_all.join(ev)
            lc_new = ev_all.to_lc(dt=dt, tstart=ev_all.time[0] - 0.5 * dt, tseg=ev_all.time[-1] - ev_all.time[0])
            # T_exp=lc_new.time[-1]-lc_new.time[0]
            temp = get_LS(lc_new.time, lc_new.counts, freq=freq, trial=k_trial)
            writer.writerows([temp[-1]])
            k_trial += 1
# get_lc_onesource_const(3,'19',num_trials=10000)

def get_LSP_disb_CR(k,CR,num_trials=100):
    k_trial=0
    FP = [];
    period = [];
    while k_trial <num_trials:
        path = '/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep{0}/'.format(k)
        epoch_file=np.loadtxt(path+'CDFS_epoch_ep{0}.txt'.format(k))
        # path = '/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8/'
        # epoch_file = np.loadtxt(path + 'CDFS_epoch.txt')
        tstart=epoch_file[:,0];tstop=epoch_file[:,1];exptime=epoch_file[:,3]
        ev_all = EventList()
        for i in range(len(exptime)):
            dt=100
            cts_rate=0.5*CR  #实际的cts-rate应为这个的2倍
            num_bins=int(exptime[i]/dt)
            sim = simulator.Simulator(N=num_bins, mean=cts_rate, dt=dt)
            w = np.arange(1 / exptime[i], 0.5 / dt, 1 / exptime[i])
            # w = np.fft.rfftfreq(sim.N, d=sim.dt)[1:]
            # spectrum = smoothbknpo(w, [0.01, 2, 1e-2, 1e-3])
            spectrum = bending_po(w, [2.3e-3, 3.4, 0.40, 4.3e-4]) ## 这是RE J1034+396 的参数
            # spectrum =powerlaw + generalized_lorentzian(w, [1 / 1000., 1 / 10000., 0.5*np.max(powerlaw), 2])
            lc = sim.simulate(spectrum)
            lc.counts += cts_rate
            lc.counts[np.where(lc.counts<0)]=0
            ps=Powerspectrum(lc,norm='abs')
            ev = EventList()
            ev.time=sim_evtlist(lc)+tstart[i]
            # print(ev.time)RX J1301.9+2747
            ev_all=ev_all.join(ev)
        lc_new = ev_all.to_lc(dt=dt, tstart=ev_all.time[0]-0.5*dt, tseg=ev_all.time[-1]-ev_all.time[0])
        # T_exp=lc_new.time[-1]-lc_new.time[0]
        T_exp = 11000154.981141508
        dt = 100
        freq = np.arange(1 / T_exp, 0.5 / dt, 1 / (5 * T_exp))
        freq = freq[np.where(freq > 1 / 20000.)]
        # print(len(freq))
        temp=get_LS(lc_new.time, lc_new.counts, freq=freq)
        FP.append(temp[0]);period.append(temp[1])
        k_trial+=1
    result=np.column_stack((FP,period))
    np.savetxt(path+'simulation/'+'trial_out_1hr_{0}_REJ1034+396.txt'.format(cr_str[j]),result,fmt="%10.5f %10.5f")
    return ev_all

def plot_result_scatter(k_num,threshold):
    figlabel=[[0,0],[0,1],[1,0],[1,1]]
    threshold=1-threshold
    figurepath='/Users/baotong/Desktop/aas/AGN_CDFS/figure/'
    cts_rate = ['2e-5']
    # plt.figure(1, (12, 8))
    fig, axes = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(14, 10))
    for i in range(len(k_num)):
        k=k_num[i];
        path = '/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep{0}/simulation/'.format(k)
        qpo_res_file = np.loadtxt(path + 'trial_out_1hr_{0}_REJ1034+396.txt'.format(str(cts_rate[0])))
        fp=qpo_res_file[:,0]
        period=qpo_res_file[:,1]
        #
        noqpo_res_file = np.loadtxt(path + 'trial_out_{0}_REJ1034+396_noQPO.txt'.format(str(cts_rate[0])))
        fp_noqpo=noqpo_res_file[:, 0]
        period_noqpo=noqpo_res_file[:,1]

        res_file_04hr = np.loadtxt(path + 'trial_out_0.4hr_{0}_REJ1034+396.txt'.format(str(cts_rate[0])))
        fp_04hr=res_file_04hr[:,0]
        period_04hr=res_file_04hr[:,1]

        res_file_02hr = np.loadtxt(path + 'trial_out_0.2hr_{0}_REJ1034+396.txt'.format(str(cts_rate[0])))
        fp_02hr=res_file_02hr[:,0]
        period_02hr=res_file_02hr[:,1]

        res_file_01hr = np.loadtxt(path + 'trial_out_0.1hr_{0}_REJ1034+396.txt'.format(str(cts_rate[0])))
        fp_01hr=res_file_01hr[:,0]
        period_01hr=res_file_01hr[:,1]
        ax_temp=axes[figlabel[i][0],figlabel[i][1]]
        ax_temp.fill_between([1. / (2.7e-4 + 2.7e-4 / 16), 1. / (2.7e-4 - 2.7e-4 / 16)], 0, 1, facecolor='lightgrey',
                         alpha=0.5)
        ax_temp.fill_between([1. / (6.68e-4 + 6.68e-4 / 16), 1. / (6.68e-4 - 6.68e-4 / 16)], 0, 1, facecolor='lightgrey',
                         alpha=0.5)
        ax_temp.fill_between([1. / (1.38889e-3 + 1.38889e-3 / 16), 1. / (1.38889e-3 - 1.38889e-3 / 16)], 0, 1,
                             facecolor='lightgrey',
                             alpha=0.5)

        ax_temp.fill_between([1. / (2.777e-3 + 2.777e-3/ 16), 1. / (2.777e-3- 2.777e-3 / 16)], 0, 1, facecolor='lightgrey',
                         alpha=0.5)
        ax_temp.scatter(period_noqpo, fp_noqpo, marker='o',linewidths=1, s=20,color='', edgecolors='black')
        ax_temp.scatter(period,fp,marker='o',s=20,color='green')
        ax_temp.scatter(period_04hr, fp_04hr,marker='x', s=30,color='red')
        ax_temp.scatter(period_02hr, fp_02hr, marker='v', s=30, color='purple')
        ax_temp.scatter(period_01hr, fp_01hr, marker='*', s=30, color='blue')
        ax_temp.legend(['1h','0.4h','0.2h','0.1h','noQPO','QPO (P=1h)','QPO (P=0.4h)','QPO (P=0.2h)','QPO (P=0.1h)',],loc='lower right',fontsize=9)
        ax_temp.text(600,0.001,'Epoch{0}'.format(k),font1)
        if i<3:ax_temp.legend_.remove()
        # ax_temp.set_title('Epoch {0}: LS detection results'.format(k),font1)
        if (i==2 or i==3):ax_temp.set_xlabel('Period (s)',font1)
        if (i==0 or i==2):ax_temp.set_ylabel('FAP',font1)
        ax_temp.tick_params(labelsize=16)
        ax_temp.plot([0,20000],[threshold,threshold],'--')
        ax_temp.loglog()
    # axes[1,1].legend(loc=2, bbox_to_anchor=(1.05, 1.0), borderaxespad=0.)
    plt.subplots_adjust(wspace=0, hspace=0)
    # plt.tick_params(labelsize=16)
    plt.savefig(figurepath+'trial_{0}.eps'.format(cts_rate[0]),bbox_inches='tight',pad_inches=0.0)
    plt.show()
# plot_result_scatter([1,2,3,4],0.99)
# def plot_res_scatter_subplot():

def plot_result_DR(k,threshold):
    ## k 指的是epoch的序号
    ## DR means detection rate
    ## return 的数组中包含该threshold下，预测与实际观测的结果
    threshold=1-threshold
    path = '/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep{0}/simulation/'.format(k)
    CDFS_LS_res = np.loadtxt('/Users/baotong/Desktop/CDFS/fig_LS_ep{0}_ovsamp_5_baluev/LS_result_{0}.txt'.format(k))
    # print(len(CDFS_LS_res))
    cts_rate=['2e-5','5e-5','1e-4','2e-4','3e-4','4e-4','5e-4']
    true_period=[1/6.68e-4,1/2.7e-4,1/2.777e-3]
    period_str=['0.4hr','1hr','0.1hr']
    DR_04hr=[];DR_1hr=[];false_DR_04hr=[];false_DR_1hr=[];DR_01hr=[];false_DR_01hr=[];DR_noQPO_04hr=[];DR_noQPO_1hr=[];DR_noQPO_01hr=[]
    for i in range(len(cts_rate)):
        qpo_res_file_noQPO = np.loadtxt(path + 'trial_out_{0}_REJ1034+396_noQPO.txt'.format(cts_rate[i]))
        fp_noQPO = qpo_res_file_noQPO[:, 0]
        qpo_res_file_04hr = np.loadtxt(path + 'trial_out_{0}_{1}_REJ1034+396.txt'.format(period_str[0],cts_rate[i]))
        fp_04hr = qpo_res_file_04hr[:, 0]
        period_04hr = qpo_res_file_04hr[:, 1]
        qpo_res_file_1hr= np.loadtxt(path + 'trial_out_{0}_{1}_REJ1034+396.txt'.format(period_str[1],cts_rate[i]))
        fp_1hr=qpo_res_file_1hr[:,0]
        period_1hr=qpo_res_file_1hr[:,1]
        qpo_res_file_01hr = np.loadtxt(path + 'trial_out_{0}_{1}_REJ1034+396.txt'.format(period_str[2], cts_rate[i]))
        fp_01hr = qpo_res_file_01hr[:, 0]
        period_01hr = qpo_res_file_01hr[:, 1]
        DR_noQPO_04hr.append(len(np.intersect1d(np.where(fp_noQPO < threshold)[0],
                                          np.where(np.abs(1 / period_04hr - 6.68e-4) < (6.68e-4 / 16))[0])))
        DR_noQPO_1hr.append(len(np.intersect1d(np.where(fp_noQPO < threshold)[0],
                                          np.where(np.abs(1 / period_04hr - 2.7e-4) < (2.7e-4/ 16))[0])))
        DR_noQPO_01hr.append(len(np.intersect1d(np.where(fp_noQPO < threshold)[0],
                                          np.where(np.abs(1 / period_04hr - 2.777e-3) < (2.777e-3 / 16))[0])))
        DR_04hr.append(len(np.intersect1d(np.where(fp_04hr<threshold)[0],np.where(np.abs(1/period_04hr-6.68e-4)<(6.68e-4/16))[0])))
        DR_1hr.append(len(np.intersect1d(np.where(fp_1hr<threshold)[0],np.where(np.abs(1/period_1hr-2.7e-4)<(2.7e-4/16))[0])))
        DR_01hr.append(len(np.intersect1d(np.where(fp_01hr < threshold)[0],
                                          np.where(np.abs(1 / period_01hr - 2.777e-3) < (2.777e-3 / 16))[0])))
        false_DR_04hr.append(len(np.intersect1d(np.where(fp_04hr<threshold)[0],np.where(np.abs(1/period_04hr-6.68e-4)>(6.68e-4/16))[0])))
        false_DR_1hr.append(len(np.intersect1d(np.where(fp_1hr<threshold)[0],np.where(np.abs(1/period_1hr-2.7e-4)>(2.7e-4/16))[0])))
        false_DR_01hr.append(len(np.intersect1d(np.where(fp_01hr < threshold)[0],
                                               np.where(np.abs(1 / period_01hr - 2.777e-3) > (2.777e-3 / 16))[0])))


    DR_1hr=np.array(DR_1hr);DR_04hr=np.array(DR_04hr);false_DR_1hr=np.array(false_DR_1hr);false_DR_04hr=np.array(false_DR_04hr)
    DR_01hr = np.array(DR_01hr);false_DR_01hr=np.array(false_DR_01hr);DR_noQPO_1hr=np.array(DR_noQPO_1hr);DR_noQPO_04hr=np.array(DR_noQPO_04hr);
    DR_noQPO_01hr=np.array(DR_noQPO_01hr);
    # DR_1hr-=DR_noQPO_1hr;DR_04hr-=DR_noQPO_04hr;DR_01hr-=DR_noQPO_01hr

    # CDFS_LS_res = CDFS_LS_res[np.where(CDFS_LS_res[:, 5] > 0)]
    bins_CR=np.array([0,4e-5,7e-5,1.5e-4,2.5e-4,3.5e-4,4.5e-4,5.5e-4])*2
    CRhist = plt.hist(CDFS_LS_res[:, 5], bins=bins_CR, histtype='step')
    # print(CRhist)
    print(len(np.where(CDFS_LS_res[:,5]>4e-4)[0]))
    #计算大于某个流量的源有多少

    plt.close()
    cts_rate=np.array([2e-5,5e-5,1e-4,2e-4,3e-4,4e-4,5e-4])*2
    plt.figure(1)
    plt.title('Detection ({0})'.format(1-threshold))
    # plt.semilogy()
    plt.plot(cts_rate,(DR_noQPO_1hr+DR_noQPO_01hr+DR_noQPO_04hr)/3000)
    plt.plot(cts_rate, DR_01hr / 1000., marker='v', color='blue')
    plt.plot(cts_rate,DR_04hr/1000.,marker='v',color='r')
    plt.plot(cts_rate, DR_1hr/1000., marker='v',color='g')
    plt.plot(cts_rate, false_DR_01hr / 1000., marker='o', linestyle='--', color='blue')
    plt.plot(cts_rate,false_DR_04hr/1000.,marker='o',linestyle='--',color='r')
    plt.plot(cts_rate, false_DR_1hr/1000., marker='o',linestyle='--',color='g')
    plt.legend(['blank FDR','DR P=0.1hr','DR P=0.4hr','DR P=1hr','FDR P=0.1hr','FDR P=0.4hr','FDR P=1hr'])
    # plt.xlabel('Photon flux (cts/s)', font1)
    plt.ylabel('Detection rate', font1)
    plt.tick_params(labelsize=16)
    # plt.subplot(212)
    # # print(DR_04hr)
    # # print(false_DR_04hr)
    # plt.plot(cts_rate, DR_01hr / (1 + DR_01hr + false_DR_01hr), marker='*', color='blue')
    # plt.plot(cts_rate,DR_04hr/(1+DR_04hr+false_DR_04hr),marker='*',color='r')
    # plt.plot(cts_rate, DR_1hr / (1 +DR_1hr+ false_DR_1hr), marker='*', color='g')
    # plt.legend(['P=0.1hr','P=0.4hr','P=1hr'])
    # # plt.xlabel('Photon flux (cts/s)', font1)
    # plt.ylabel('DR/(FR+DR)', font1)
    plt.show()
    # plt.close()

    sim_NUM=np.sum( CRhist[0] * (DR_04hr+false_DR_04hr) / 1000)
    sim_real_NUM = np.sum(CRhist[0] * DR_04hr / 1000)
    # print(DR_01hr)

    src_id=CDFS_LS_res[:,0];conf=CDFS_LS_res[:,1];period=CDFS_LS_res[:,2]
    good_id=[reduce(np.intersect1d, (np.where(conf < threshold), np.where(period>203.),np.where(np.abs(period - 707.) > 5),
                                   np.where(np.abs(period - 2*707.) > 10),np.where(np.abs(period - 999.) > 5),
                                   np.where(np.abs(period - 2*1000.) > 10),np.where(np.abs(period - 200.) > 10),
                                   np.where(period<20000)))]
    CDFS_LS_res=CDFS_LS_res[good_id]
    np.savetxt('/Users/baotong/Desktop/CDFS/fig_LS_ep{0}_ovsamp_5_baluev/LS_good_result_{0}.txt'.format(k),CDFS_LS_res,fmt='%10d %15.5f %15.5f %10d %10d %15.10f')
    det_NUM=len(good_id[0])
    src_id=src_id.astype('int')
    GOODID=src_id[good_id]
    # print(len(GOODID))
    print(sim_real_NUM)
    return [sim_NUM,det_NUM,sim_real_NUM]

plot_result_DR(4,0.9)

def plot_LS_sim_det(k):
    x=[];y1=[];y2=[];y3=[]
    threshold_range=np.linspace(0.9,1.0,50)
    for i in range(len(threshold_range)):
        [sim,det,sim_real]=plot_result_DR(k,threshold_range[i])
        y1.append(sim)
        y2.append(det)
        y3.append(sim_real)
    x=threshold_range;y1=np.array(y1);y2=np.array(y2)
    y1_down=y1.astype(int)
    y1_up=y1_down+1
    y4 = y3 / y1 * y2
    plt.figure(1,(10,8))
    plt.subplot(211)
    plt.plot(x, y1, marker='*', color='r')
    plt.plot(x,y2,marker='*',color='g')
    plt.plot(x, y3, marker='*', color='black')
    plt.legend(['Sim', 'Detect','Sim_trueP'])
    plt.plot(x,y1_down,'--')
    plt.plot(x,y1_up,'--')
    plt.semilogy()
    plt.fill_between(x,y1_down,y1_up,facecolor='grey',alpha=0.5)
    plt.xlabel('Confidence', font1)
    plt.ylabel('Number of source', font1)
    plt.tick_params(labelsize=16)

    plt.subplot(212)
    plt.plot(x,y2/y1)
    plt.savefig('/Users/baotong/Desktop/CDFS/sim_fig/conf_{0}_04hr.eps'.format(k))
    plt.show()
# plot_LS_sim_det(4)