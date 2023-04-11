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


