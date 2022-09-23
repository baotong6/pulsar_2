# -*- coding: utf-8 -*-
"""
Created on Sun Sep 20 18:13:40 2015
@author: Felix Darvas
@author: baotong
compute the Gregory-Laredo algorithm on arival times
This function computes the likelihood of a set of arrival times originating
from a periodic system rather than constant rate (poisson) background noise
based on
Gregory, P. C. and Thomas. J. Loredo, 1992,
"A New Method For The Detection Of A Periodic Signal Of Unknown Shape And Period"
in
The Astrophysical Journal, Astrophysical J., 398, p.146
inputs:
Tlist    -  list of arrival times, numpy int array
m_max    -  max number of bins typically 12-15, we use 12 as default
w_range  -  frequency range to scan numpy float array of frequency values
           default is  w_lo=20*pi/T at delta w = pi/T to w_hi=pi*N/T
           where N=#arrival times, T=observation time
ni       - number of integration steps, default ni=10
parallel - use parallel execution - default is off
outut:
O_period - Odds ratio for a periodic process vs. constant rate process
p_period - probability of a periodic process 0<=p_period<=1
m_opt    - optimal bin size 1<= m_opt <=m_max
S        - The probability spectrum
w        - The frequency range for S
w_peak   - the expected frequency of the process
w_conf   - 95% confidence interval of w_peak
"""
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp
import functools
import datetime
from astropy.io import fits
import sys
import os
import string
#import correct as correct
import random
from scipy import optimize
from sympy import *
from scipy.optimize import fsolve
from itertools import chain
import hawkeye as hawk
import rednoise

idlist=[217,414,185,366,423,232,263,331,273,317,162,252,283,290,198,312,229]
path = '/Users/baotong/Desktop/period_Tuc/txt_startover/txt_all_obs_p90/'
for srcid in idlist[14:15]:
    print(srcid)
    epoch_file = path + 'epoch_src_' + str(srcid)+'.txt'
    evt_all=rednoise.simulate_srcinfo(srcid)
    src_evt=np.column_stack((evt_all.time,np.zeros(len(evt_all.time))+1000.))
    w_range = 2 * np.pi * np.arange(1 / 50000., 1 / 10000, 1e-7)
    # print(src_evt)
    res=hawk.GL_algorithm_single.write_result(src_evt,epoch_file,w_range,dataname='1',if_filter=False)
    print(res)