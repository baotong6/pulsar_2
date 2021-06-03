import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import pandas as pd
import sys
import os
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
font_prop = font_manager.FontProperties(size=16)
import warnings
from functools import reduce

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

    x = time
    y = flux
    # dy=np.sqrt(y)
    # plt.scatter(x,y)
    # plt.show()

    # LS = LombScargle(x, y, dy = 1, normalization = 'standard', fit_mean = True,
    #                  center_data = True).power(freq, method = 'cython')
    LS = LombScargle(x, y,normalization = 'psd')
    # LS = LombScargle(x, y, dy, normalization='psd')
    power = LS.power(freq)

    # print('freq_num={0}'.format(len(freq)))
    # FP=LS.false_alarm_probability(power.max(),minimum_frequency = freq[0], maximum_frequency = freq[-1],method='baluev')
    # FP_99 = LS.false_alarm_level(0.0027, minimum_frequency = freq[0], maximum_frequency = freq[-1],method='baluev')
    # FP_90 = LS.false_alarm_level(0.05,  minimum_frequency=freq[0],
    #                              maximum_frequency=freq[-1], method='baluev')
    # FP_68 = LS.false_alarm_level(0.32, minimum_frequency=freq[0],
    #                              maximum_frequency=freq[-1], method='baluev')
    # if FP<0.01:print(dataname)
    # plt.title('FP={0}'.format(FP))
    plt.semilogx()
    plt.plot(freq, power)
    # print(1./freq[np.where(power==np.max(power))])
    # plt.plot([freq[0], freq[-1]], [FP_99, FP_99], '--'
    # plt.plot([freq[0], freq[-1]], [FP_90, FP_90], '--')
    # plt.plot([freq[0], freq[-1]], [FP_68, FP_68], '--')
    # plt.show()
    res=1e2*power
    res=np.round(res,2)
    # res=res[::smooth]
    # np.savetxt(path+'LS_simres_{0}.txt'.format(trial+1),res,fmt='%10.5f')
    period=1. / freq[np.where(power == np.max(power))][0];LSP=np.max(power)
    period=np.round(period,2);LSP=np.round(1e2*LSP,2)
    if LSP<18.36:
        plt.close()
    return [period,LSP]

def sim_evtlist(lc):
    num_evt_bin=np.random.poisson(lc.counts)
    evt_all=[]
    i=0
    while i < len(num_evt_bin):
        evt_all=np.concatenate((evt_all,(np.random.uniform(lc.time[i]-lc.dt/2,lc.time[i]+lc.dt/2,num_evt_bin[i]))))
        i+=1
    return evt_all

def slow_sim_const(k,src_index,num_trials=2):
    def pfold(time,period,binnumber=10):
        turns=time/period
        HIST=np.zeros(binnumber)
        phase=turns-turns.astype('int')
        for i in range(len(phase)):
            HIST[int(phase[i]*binnumber)]+=1

        x = np.linspace(0, 1, binnumber + 1)
        x = x[:-1]

        x2 = np.concatenate((x, x + 1));
        y2 = np.concatenate((HIST, HIST))
        plt.figure(1,(9,6))
        plt.step(x2,y2)

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


    with open(path + '/simulation/{0}_max_LSP_const.txt'.format(src_index), 'a+') as f:
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
            temp = get_LS(lc_new.time, lc_new.counts, freq=freq, trial=k_trial)
            f.writelines((str(temp[0]) + '        ' + str(temp[1]) + '\n'))
            if temp[1]>18.36:
                print('caution')
                plt.savefig(path + '/simulation/LSP_{0}_trial_{1}.eps'.format(src_index,k_trial))
                plt.close()
                np.savetxt(path + '/simulation/{0}_trial_{1}_event.txt'.format(src_index,k_trial),ev_all.time)
                pfold(ev_all.time,temp[0],10)
                plt.savefig(path + '/simulation/pfold_{0}_trial_{1}.eps'.format(src_index, k_trial))
                break

            k_trial += 1
    f.close()
slow_sim_const('3','89',num_trials=4000)