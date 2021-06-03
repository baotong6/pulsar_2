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
import useful_functions as func

font1=func.font1
def build_psd(x,p,x_2=1,p_2=1,type='bendp'):
    ##
    if type=='bendp':
        psd=func.bending_po(x,p)
    elif type=='lorentz':
        psd=func.generalized_lorentzian(x,p)
    elif type=='bendp+lorentz':
        psd=func.bending_po(x,p)+func.generalized_lorentzian(x_2,p_2)
    return psd
def make_freq_range(dt,epoch_file):
    (TSTART, TSTOP, OBSID, exptime)=func.read_epoch(epoch_file)
    T_tot=TSTOP[-1]-TSTART[0]
    w = np.arange(1 / T_tot, 0.5 / dt, 1 / T_tot)
    return w
def make_lc_from_psd(psd,cts_rate,dt,epoch_file):
    (TSTART, TSTOP, OBSID, exptime)=func.read_epoch(epoch_file)
    T_tot=TSTOP[-1]-TSTART[0]
    num_bins=int(T_tot/dt)
    sim = simulator.Simulator(N=num_bins + 1, mean=cts_rate, dt=dt)
    lc=sim.simulate(psd)
    lc.time=lc.time+TSTART[0]
    lc.counts += cts_rate
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
    lenbin=1000
    epoch_89='/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep3/epoch_src_89.txt'
    w=make_freq_range(dt=lenbin,epoch_file=epoch_89)
    psd_model=build_psd(x=w,p=[4.3e-4,3.4,0.4,2.3e-3],type='bendp')
    lc=make_lc_from_psd(psd=psd_model,cts_rate=0.0001*lenbin,dt=lenbin,epoch_file=epoch_89)
    print('counts={0}'.format(np.sum(lc.counts)))
    ps=plot_psd(lc)
    lc_evt=Lightcurve(time=lc.time,counts=lc.counts,dt=lc.dt,gti=lc.gti)
    lc_evt.counts=np.random.poisson(lc_evt.counts)
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

def sim_bunch_lc(dt,epoch_file,k_trials,path):

# def make_evt_fromlc(lc):
#     evt=simulator.simluate_times(lc)
#     return evt
if __name__=='__main__':
    test_somefunc()