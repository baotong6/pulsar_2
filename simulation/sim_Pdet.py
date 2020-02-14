#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import string
#import correct as correct
from datetime import datetime
from scipy.interpolate import lagrange
from scipy import interpolate
from scipy.optimize import curve_fit
import read_csv as data

path='/Users/baotong/Desktop/period_LW/'

def get_Z2(dataname,freq):
    time = np.loadtxt(path + 'txt_all_obs/' + dataname + '.txt')[:,0]
    N=len(time)
    def turns(t,f):
        ti=t-t[0]
        v=f
        #p_test = 1.0/5.500550055005501e-06
        p=1.0/v
        freq = v # + ti * vdot + ti * ti / 2.0 * vddot
        preq = 1.0 / freq
        turns=v*ti
        INT_turns=np.trunc(turns)
        turns=turns-INT_turns
        turns = 2.0*np.pi*turns #+ vdot * ti * ti / 2.0 + vddot * ti * ti * ti / 6.0
        return turns
    # Z2=[]
    # for fi in freq:
    #     Z2.append((2.0 / N)*(sum(np.cos(turns(time,fi)))**2+sum(np.sin(turns(time,fi))**2)))
    fi=freq
    Z2=(2.0 / N) * (sum(np.cos(turns(time, fi))) ** 2 + sum(np.sin(turns(time, fi)) ** 2))
    cts=len(time)
    return [Z2,cts]

ID=data.ID_LW
period=data.P_LW
backscale=data.net_percent_LW
for i in range(len(ID)):
    dataname=str(int(ID[i]))
    if dataname[-3:] == '001' or dataname[-3:] == '002':
        dataname = dataname[:-3]

    [Z2,cts]=get_Z2(dataname,1.0/period[i])
    bkg_cts=cts*(1-backscale[i])
    A_rms=(Z2/cts)**0.5*(cts/(cts-bkg_cts))
    A0=2**0.5*A_rms
    print(A_rms)
    # print(A0)