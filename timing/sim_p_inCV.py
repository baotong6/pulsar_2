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
import scipy.signal as signal

#预测轨道倾角对探测到周期信号的影响

def WD_M_R(M):
    #输入白矮星质量，输出白矮星半径
    R=7.8*1e8*((1.44/M)**(2.0/3.0)-(M/1.44)**(2.0/3.0))**0.5
    return R

path='/Users/baotong/Desktop/period_gc/'

sim_info=fits.open(path+'re_model.fit')
M1=sim_info[1].data.field(0)
R1=WD_M_R(M1)/(6.955*1e10)
M2=sim_info[1].data.field(1)
R2=sim_info[1].data.field(2)
Period=sim_info[1].data.field(3)
Sep=sim_info[1].data.field(4)

cosi=(R2+R1)/Sep
prob=(0.5*np.pi-np.arccos(cosi))/(0.5*np.pi)

print(M2)
print(R1)
print(R2)
print(Sep)
print(cosi)
print(np.arccos(cosi)*180/np.pi)
print(prob)
#plt.scatter(Period,M2)
plt.figure(1)
plt.scatter(Period,prob,marker='+',alpha=0.5
            ,color='green')
plt.xlabel('period(h)')
plt.ylabel('probability')
plt.show(1)
plt.figure(2)
plt.xlim(0.8,0.0)
plt.ylim(0.0,0.8)
plt.scatter(M2,R2,marker='+')
plt.xlabel('mass')
plt.ylabel('radius')
#plt.show(2)
