'''
Author: baotong && baotong@smail.nju.edu.cn
Date: 2024-03-06 22:22:37
LastEditors: baotong && baotong@smail.nju.edu.cn
LastEditTime: 2024-03-06 22:26:10
FilePath: /pulsar/mod_result/wwztest1.py
Description: 

Copyright (c) 2024 by baotong, All Rights Reserved. 
'''
# # %%
# import warnings
# warnings.filterwarnings('ignore')
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import pandas as pd
import sys
import os
from scipy import interpolate
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from scipy.optimize import curve_fit
import astropy.units as u
import astropy.constants as c
import libwwz
from astropy.timeseries import LombScargle
import stingray as sr
from stingray.events import EventList
from stingray.lightcurve import Lightcurve
from stingray import Lightcurve, Crossspectrum, sampledata,Powerspectrum,AveragedPowerspectrum
from stingray.simulator import simulator, models
def get_hist(t, len_bin):
    ###将输入的time信息，按照len_bin的长度输出为lc
    t_test = t-t[0]
    a = [0 for i in range(int(t_test[-1] / len_bin) + 1)]
    for i in range(len(t_test)):
        a[int(t_test[i] / len_bin)] += 1
    a = np.array(a)
    return a

dt = 2  # seconds
long_dt = 2  # seconds
long_exposure = 160.  # seconds
long_times = np.arange(0, long_exposure, long_dt)  # seconds
# In count rate units here
long_signal = 100 * np.sin(2.*np.pi*long_times/10) + 1000
# Multiply by dt to get count units, then add Poisson noise
long_noisy = np.random.poisson(long_signal*long_dt)
long_lc = Lightcurve(long_times, long_noisy)

# help(libwwz.wwt)


freq=np.arange(1/10.,1./0.5,1/5000)
bin_len=0.1
path='/Users/baotong/Desktop/FRBtime/'
filename='FRB240114A_0306.txt'
time=np.loadtxt(path+filename)*86400
print(time)
T_exp = time[-1] - time[0]
# freq=np.arange(1/T_exp,0.5/bin_len,1/(5*T_exp))            
flux=get_hist(time,bin_len)
x=np.arange(bin_len/2.,(time[-1]-time[0])+bin_len/2.,bin_len)

[Tau, Freq, WWZ, AMP, COEF, NEFF]=libwwz.wwt(timestamps=x,magnitudes=flux,time_divisions=10,
                                             freq_params=[freq[0],freq[-1],freq[-1]-freq[-2],False],decay_constant=1e-3)
     
fig, ax = plt.subplots()
ax.set_title('test title')
libwwz.plot_methods.linear_plotter(ax,Tau,Freq,WWZ)
plt.show()



