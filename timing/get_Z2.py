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
import read_data as data
from scipy.fftpack import fft,ifft
import scipy.signal as ss
import scipy.stats as stats
import random
import read_data as data

#采样频率数
v=44999
##(1-p)**N=0.99
#置信度99%
p99=1-0.99**(1.0/v)
p90=1-0.90**(1.0/v)
x1=[]
y1=[]
x2=[]
y2=[]
for x in np.linspace(1,1000,5000):
    if stats.chi2.cdf(x,2)!=0:
        y1.append(np.abs((1-stats.chi2.cdf(x,20-1))-p99))
        y2.append(np.abs((1-stats.chi2.cdf(x,20-1))-p90))
        #print y1
        x1.append(x)
        x2.append(x)
    else:
        break
y1=np.array(y1)
x1=np.array(x1)
y2=np.array(y2)
x2=np.array(x2)

#print min(y1)
print(x1[np.where(y1==min(y1))])
print ("%e" %(1-stats.chi2.cdf(x1[np.where(y1==min(y1))],2)))
print(p99)

#print min(y2)
print(x1[np.where(y2==min(y2))])
print ("%e" %(1-stats.chi2.cdf(x1[np.where(y2==min(y2))],2)))
print(p90)



# path='/Users/baotong/Desktop/period/'
# dataname="src2338_I_lc.txt"
# time=[];energy=[];obs_ID=[]
# time=data.time
# energy=data.energy
# obs_ID=data.obs_ID
# t=data.t
# E=data.E
# #time=t[0]
# N=len(time)
# dict=data.dict
# epoch=np.loadtxt(path+"SgrA_I_epoch.txt")
# exptime=epoch[:,1]-epoch[:,0]
# obstime=epoch[:,0:2]
# cts_s=len(t[:])/exptime
# print cts_s
# print np.where(exptime==max(exptime))
# #plt.scatter(obstime[:,0],cts_s)
# plt.hist(t[9],bins=int(exptime[9]/4000.0))
# plt.show()
# def get_R2(q,M2=0.5,Porb=10.):
#     R2=0.2478*M2**0.33*Porb**(2/3.)*((q*(1+q))**0.333)/(0.6*q**0.667+np.log(1+q**0.333))
#     return R2
#
# q=np.linspace(0,0.9,1000)
# R2=get_R2(q)
# plt.scatter(q,R2)
# plt.show()