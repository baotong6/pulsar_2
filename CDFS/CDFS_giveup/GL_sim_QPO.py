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
from astropy.stats import poisson_conf_interval
font1 = {'family': 'Normal',
         'weight': 'normal',
         'size': 16, }
def sim_evtlist(lc):
    num_evt_bin=np.random.poisson(lc.counts)
    evt_all=[]
    i=0
    while i < len(num_evt_bin):
        evt_all=np.concatenate((evt_all,(np.random.uniform(lc.time[i]-lc.dt/2,lc.time[i]+lc.dt/2,num_evt_bin[i]))))
        i+=1
    return evt_all

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

    y2_err=np.array(poisson_conf_interval(y2,interval='frequentist-confidence'))
    y2_err[0]=y2-y2_err[0]
    y2_err[1]=y2_err[1]-y2

    plt.errorbar(x2 - 0.5 / binnumber, y2, yerr=y2_err, fmt='.', capsize=1, elinewidth=1, ecolor='red')

    plt.figure(1,(9,6))
    plt.step(x2,y2)
    plt.show()


def make_evt_nobkg(cts_rate,dt,epoch_file,lc_bin=100,fixpds=True,lc=0):
    tstart = epoch_file[:, 0];tstop = epoch_file[:, 1]
    ID=epoch_file[:,2];exptime = epoch_file[:, 3]
    T_tot=tstop[-1]-tstart[0]
    cts_rate=cts_rate/2*dt ##这里cts_rate单位都是per bin
    num_bins = int(T_tot / dt)
    freq=np.arange(1/T_tot,0.5/dt,1/(5*T_tot))
    freq=freq[np.where(freq > 1 / 20000.)]
    if fixpds:
        lc=lc
    else:
        print('not fix!')
        sim = simulator.Simulator(N=num_bins+1, mean=cts_rate, dt=dt)
        w = np.arange(1 / T_tot, 0.5 / dt, 1 / T_tot)
        spectrum = bending_po(w, [2.3e-4, 3.4, 0.40, 4.3e-4]) + generalized_lorentzian(w, [4.01e-4, 4.01e-4 / 16, 0.1, 2])
        lc = sim.simulate(spectrum)
        lc.time=lc.time+tstart[0]
        lc.counts += cts_rate
        lc.counts[np.where(lc.counts < 0)] = 0
        lc.gti=[[lc.time[0],lc.time[-1]]]
    plt.plot(w,spectrum)
    plt.loglog()
    plt.show()
    plt.plot(lc.time,lc.counts)
    plt.show()
    ps = Powerspectrum(lc, norm='abs')
    fig, ax1 = plt.subplots(1, 1, figsize=(9, 6), sharex=True)
    ax1.plot(ps.freq, ps.power, lw=2, color='blue')
    ax1.set_ylabel("Frequency (Hz)", fontproperties=font1)
    ax1.set_ylabel("Power (raw)", fontproperties=font1)
    # ax1.set_yscale('log')
    ax1.tick_params(axis='x', labelsize=16)
    ax1.tick_params(axis='y', labelsize=16)
    ax1.tick_params(which='major', width=1.5, length=7)
    ax1.tick_params(which='minor', width=1.5, length=4)
    for axis in ['top', 'bottom', 'left', 'right']:
        ax1.spines[axis].set_linewidth(1.5)
    plt.loglog()
    plt.show()

    ev_all = EventList();
    ev_all_bkg = EventList()
    for i in range(len(ID)):
        lc_cut = lc.truncate(start=tstart[i], stop=tstop[i], method='time')
        ev = EventList()
        ev.time = sim_evtlist(lc_cut)
        ev_all = ev_all.join(ev)

    lc_src = ev_all.to_lc(dt=lc_bin, tstart=tstart[0], tseg=T_tot)
    ps = Powerspectrum(lc_src, norm='abs')
    fig, ax2 = plt.subplots(1, 1, figsize=(9, 6), sharex=True)
    ax2.plot(ps.freq, ps.power, lw=2, color='blue')
    ax2.set_ylabel("Frequency (Hz)", fontproperties=font1)
    ax2.set_ylabel("Power (raw)", fontproperties=font1)
    # ax1.set_yscale('log')
    ax2.tick_params(axis='x', labelsize=16)
    ax2.tick_params(axis='y', labelsize=16)
    ax2.tick_params(which='major', width=1.5, length=7)
    ax2.tick_params(which='minor', width=1.5, length=4)
    for axis in ['top', 'bottom', 'left', 'right']:
        ax2.spines[axis].set_linewidth(1.5)
    plt.loglog()
    plt.show()

    return [ev_all.time,lc_src]

if __name__ == '__main__':
    path='/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep3/'
    epoch_file = np.loadtxt(path + 'epoch_src_89.txt')
    cts_rate=1e-3;dt=100

    tstart = epoch_file[:, 0];tstop = epoch_file[:, 1]
    ID=epoch_file[:,2];exptime = epoch_file[:, 3]
    T_tot=tstop[-1]-tstart[0]
    cts_rate=cts_rate/2*dt ##这里cts_rate单位都是per bin
    num_bins = int(T_tot / dt)
    sim = simulator.Simulator(N=num_bins+1, mean=cts_rate, dt=dt)
    w = np.arange(1 / T_tot, 0.5 / dt, 1 / T_tot)
    spectrum = bending_po(w, [2.3e-4, 3.4, 0.40, 4.3e-4]) + generalized_lorentzian(w, [4.01e-4, 4.01e-4 / 16, 200, 2])
    lc = sim.simulate(spectrum)
    lc.time=lc.time+tstart[0]
    lc.counts += cts_rate
    lc.counts[np.where(lc.counts < 0)] = 0
    lc.gti=[[lc.time[0],lc.time[-1]]]
    # plt.plot(w,spectrum)
    k_trials=1
    for i in range(k_trials):
        [time,lc_src]=make_evt_nobkg(cts_rate=1e-3, dt=3.2, epoch_file=epoch_file,fixpds=False,lc=lc)
        import lc_analysis as plot_LS
        import GL_algorithm_withepoch as GL
        freq=np.arange(1/2e6,1/250,1/(10*2e6))
        freq=freq[np.where(freq>1/20000.)]
        plt.close()
        plot_LS.get_LS(lc_src.time,lc_src.counts,freq=freq)

        pfold(time,period=1/4.01e-4)

    # w_range=np.arange(1/3000.,1/2000.,1e-7)
    # GL_R=GL.compute_GL(time, epoch_file=epoch_file,m_max=20, w_range=w_range, ni=10, parallel=True)
    # Prob=GL_R[1]
    # wpeak=GL_R[5]
    # mopt=GL_R[2]
    # wconf_lo=GL_R[7][0]
    # wconf_hi=GL_R[7][1]
    # res=np.column_stack(Prob,wpeak,mopt,wconf_hi,wconf_lo)
    # print(res)