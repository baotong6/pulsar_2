#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import csv
import os
import re
import string
#import correct as correct
from datetime import datetime
from scipy.interpolate import lagrange
from scipy import interpolate
from scipy.optimize import curve_fit
from scipy.special import comb, perm
import pandas as pd
from astropy.stats import poisson_conf_interval
from astropy.timeseries import LombScargle
import scipy
import linecache

path='/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep{0}/simulation/'.format('3')
pathsim='/Volumes/pulsar/89_LS_sim_noQPO/'
i=1;num_trials=100000;smooth=100

with open(path+"test.csv","a+") as csvfile:
    T_exp = 11000154.981141508;
    dt = 100
    freq = np.arange(1 / T_exp, 0.5 / dt, 1 / (5 * T_exp))
    freq = freq[np.where(freq > 1 / 20000.)]
    header=freq[::smooth]
    header=header.astype('str')
    writer = csv.writer(csvfile)
    writer.writerow(header)
    while i <num_trials:
        line=np.loadtxt(pathsim+'LS_simres_{0}.txt'.format(i))
        # print(line)
        writer.writerows([line])
        i+=1
