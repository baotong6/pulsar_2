#!/bin/bash
# -*- coding: utf-8 -*-
#plt.rcParams['font.sans-serif'] = ['Microsoft YaHei']
#用于解决matplotlib中文乱码问题
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
from scipy.fftpack import fft,ifft
import scipy.signal as ss
import scipy.stats as stats
import random

#path='/Volumes/halo/chandra_obs_hr_m04_for_baotong/'
path='/Users/baotong/xmm_obs/'
#path='/Users/baotong/Desktop/period/'
#path='/Volumes/halo/redgiant/'
#dataname="Muno4.txt"
#dataname="Muno6.txt"
dataname='src22pn'

def get_Z2(dataname,freq):
    time = np.loadtxt(path + 'txt/' + dataname + '.txt')[:,0]
    #time=   np.loadtxt(path + 'txt/' + dataname + '.txt')[:,0]
    #time = np.loadtxt(path + 'txt/' + dataname + '_s.txt')[:, 0]
    #time = np.loadtxt(path + 'txt_G/' + dataname + '.txt')[:, 0]
    #time=np.loadtxt(path + 'txt/evt_list10956.txt')[:,0]
    N=len(time)
    def turns(t,f):
        ti=t-t[0]
        #print ti
        v=f
        #p_test = 1.0/5.500550055005501e-06
        p=1.0/v
        # pdot=-vdot/(v*v)
        # vddot = 2.0 * pdot * pdot / (p * p * p)
        freq = v # + ti * vdot + ti * ti / 2.0 * vddot
        preq = 1.0 / freq
        turns=v*ti
        INT_turns=np.trunc(turns)
        turns=turns-INT_turns
        turns = 2.0*np.pi*turns #+ vdot * ti * ti / 2.0 + vddot * ti * ti * ti / 6.0
        #print turns
        #turns=list(turns)
        #turns=np.array(turns[0])
        #初始相位
        return turns
    #print time
    Z2=[]
    plt.hist(time-time[0],bins=60,histtype='step')
    plt.show()

    for fi in freq:
        Z2.append((2.0 / N)*(sum(np.cos(turns(time,fi)))**2+sum(np.sin(turns(time,fi))**2)))
    cts=len(time)
    return [Z2,cts]

# border=100
# vary=np.array([i for i in range(0,border)])
# freq=1/10000.+vary*(1/100000.)


def make_period_range(pmin, pmax, expT):
    P = [pmin]
    while P[-1] < pmax:
        dP = 0.05 * P[-1] ** 2 / (expT - P[-1])
        P.append(P[-1] + dP)
    return np.array(P)


freq = 1. / make_period_range(2000, 10000, 2e4)
# border=10000
# vary=np.array([i for i in range(0,border)])
# freq=1.e-4+vary*1.e-6
p99_Z2 = 27.6244
p90_Z2 = 22.9232

#
# evt_list_ID=np.arange(1,49)
# evt_list_ID=np.delete(evt_list_ID,[20,31,36,37,38,39,40,41,45])
# evt_list_name=[]
# for i in range(len(evt_list_ID)):
#     if evt_list_ID[i]<10:
#         evt_list_name.append('src'+'0'+str(evt_list_ID[i]))
#     else:
#         evt_list_name.append('src'+str(evt_list_ID[i]))

# evt_list_name=['src1mos1','src2mos1','src3mos1','src4mos1','src5mos1',
#                'src6mos1','src7mos1','src8mos1','src9mos1','src10mos1',
#                'src11mos1','src12mos2','src13mos2','src14mos1','src15mos1','src15mos1','src16mos2']

# evt_list_name=['src1mos2','src2mos2','src3mos2','src4mos2','src5mos2',
#                'src6mos2','src7mos2','src8mos2','src9mos2','src10mos2',
#                'src11mos2','src12mos2','src13mos2','src14mos2','src15mos2','src15mos2','src16mos2']
# #
# evt_list_name=['src1pn','src2pn','src3pn','src4pn','src5pn',
#                'src6pn','src7pn','src8pn','src9pn','src10pn',
#                'src11pn','src12pn','src14pn','src15pn','src15pn','src16pn']

evt_list_name=['src22pn']
for i in range(len(evt_list_name)):
    dataname=evt_list_name[i]
    Z2=get_Z2(dataname,freq)[0]
    cts=get_Z2(dataname,freq)[1]
    Z2=np.array(Z2)
    plt.figure(1,(7,7))
    plt.title(dataname+',cts={0}'.format(cts))
    plt.semilogx(freq,[p99_Z2 for i in range(len(Z2))],'--',color='black')
    plt.semilogx(freq,[p90_Z2 for i in range(len(Z2))],'--',color='red')
    plt.step(freq,Z2,color='black')
    plt.text(0.0005,p99_Z2+1,"99%")
    plt.text(0.0005,p90_Z2+1,"90%")
    plt.show()
    #plt.savefig(path + 'fig_Z/' + dataname + '.eps')
    #plt.close()
#
# dataname='src22pn'
# Z2=get_Z2(dataname,freq)[0]
# cts=get_Z2(dataname,freq)[1]
# Z2=np.array(Z2)
# plt.figure(1,(7,7))
# plt.title(dataname+',cts={0}'.format(cts))
# #plt.title('10956'+',cts={0}'.format(cts))
# plt.semilogx(freq,[p99_Z2 for i in range(len(Z2))],'--',color='black')
# plt.semilogx(freq,[p90_Z2 for i in range(len(Z2))],'--',color='red')
# plt.step(freq,Z2,color='black')
# plt.text(0.0005,p99_Z2+1,"99%")
# plt.text(0.0005,p90_Z2+1,"90%")
# plt.savefig(path + 'figZ/10-1000/' + dataname + '.eps')
# plt.show()

#plt.close()

# Z2_p=Z2-signal_Z2
# p_possible=1.0/freq[np.where(Z2_p>0)]
# print p_possible
# np.savetxt(path+dataname[3:7]+"/p_possible.txt",p_possible)0
