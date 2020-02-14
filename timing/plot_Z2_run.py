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
import read_data as data
from scipy.fftpack import fft,ifft
import scipy.signal as ss
import scipy.stats as stats
import random
import fits_to_txt

path='/Users/baotong/Desktop/period/'
# dataname="3525.txt"
#dataname="Muno4.txt"
#dataname="Muno6.txt"
epochname='SgrA_I_epoch.txt'

def sim_time(dataname):
    time = data.get_data(dataname)[0]
    energy = data.get_data(dataname)[1]
    obs_ID = data.get_data(dataname)[2]
    t = data.get_data(dataname)[3]
    E = data.get_data(dataname)[4]
    dict = data.get_data(dataname)[5]

    t = np.array(t)

    sim_t = t - t
    for i in range(len(t)):
        sim_t[i] = np.random.random(len(t[i])) * (t[i][-1] - t[i][0]) + t[i][0]
    return sim_t

def get_Z2(dataname,freq):
    time=data.get_data(dataname)[0]
    energy=data.get_data(dataname)[1]
    obs_ID=data.get_data(dataname)[2]
    t=data.get_data(dataname)[3]
    E=data.get_data(dataname)[4]
    dict=data.get_data(dataname)[5]

    sss=0
    for index in t:
        for i in range(len(index)):
            sss+=1
    if sss!=len(time):
        print('error')
        print('sss=%d %sss')
        print('len(time)=%d len(time)')
    #time=t[0]
    # #
    delete_i=[]
    use_ID=[]
    #dict[242,2951][0]

    for i in dict:
        if dict[i][-2]<5000:
            delete_i.append(i)
    for item in delete_i:
        dict.pop(item)
    for i in dict:
        use_ID.append(i)

    epoch0 = np.array([2943,3663])
    epoch1=np.array([2943,3663,3392,3393,3665])
    epoch2=np.array([5950,5951,5952,5953,5954])
    epoch3=np.array([9169,9170,9171,9172])
    epoch4=np.array([10556,11843])
    use_ID=epoch1
    #use_ID=np.concatenate((epoch2,epoch3,epoch4))
    #use_ID=epoch2
    #use_ID=epoch4
    #print use_ID
    #use_ID=[3392,3393]
    # #print dict[10556]
    time=np.concatenate((dict[use_ID[0]][0],dict[use_ID[1]][0]))

    for i in range(2,len(use_ID)):
        time=np.concatenate((time,dict[use_ID[i]][0]))

    N = len(time)

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

    for fi in freq:
        Z2.append((2.0 / N)*(sum(np.cos(turns(time,fi)))**2+sum(np.sin(turns(time,fi))**2)))
    cts=len(time)
    return [Z2,cts]

def get_Z2_sim(dataname):
    sim_t=sim_time(dataname)
    simtime=data.tran_t_time(sim_t)
    Zr_temp = []
# simtime=np.reshape(sim_time(dataname),(1,len(time)))
    for fi in freq:
        Zr_temp.append((2.0 / N)*(sum(np.cos(turns(simtime,fi)))**2+sum(np.sin(turns(simtime,fi))**2)))
    Zr_temp=np.array(Zr_temp)
    return Zr_temp

# Zr=get_Z2_sim(dataname)
# for i in range(10):
#     Zr_temp=get_Z2_sim(dataname)
#     Zr=(Zr+Zr_temp)/2.0
def get_fig_Z2(dataname):
    border=5000
    vary=np.array([i for i in range(0,border)])
    #freq=1/700.+vary*1.e-6
    freq=1.e-5+vary*1.e-6


    # border=10000
    # vary=np.array([i for i in range(0,border)])
    # freq=1.e-4+vary*1.e-6
    p99_Z2 = 23.01325331
    p90_Z2 = 18.31207802

    # p99_Z2=27.6244061
    # p90_Z2=22.92323081

    (Z2,cts)=get_Z2(dataname,freq)
    (Z2_standard,cts_standard)=get_Z2('pwn.txt',freq)
    Z2=np.array(Z2)
    Z2_standard=np.array(Z2_standard)
    plt.subplot(211)
    plt.title(dataname[0:-4]+',cts={0}'.format(cts))
    plt.semilogx(freq,[p99_Z2 for i in range(len(Z2))],'--',color='black')
    plt.semilogx(freq,[p90_Z2 for i in range(len(Z2))],'--',color='red')
    plt.step(freq,Z2,color='black')
    #plt.step(freq,Z2_standard*max(Z2)/max(Z2_standard),color='pink')
    plt.text(0.5,p99_Z2+1,"99%")
    plt.text(0.5,p90_Z2+1,"90%")
    plt.subplot(212)
    plt.semilogx(freq,[p99_Z2 for i in range(len(Z2))],'--',color='black')
    plt.semilogx(freq,[p90_Z2 for i in range(len(Z2))],'--',color='red')
    plt.step(freq,Z2-Z2_standard*max(Z2)/max(Z2_standard))
    #plt.show()
    print(dataname)
    plt.savefig(path+'fig_Z_90_ep1/'+dataname[0:-4]+"_ep1.png")
    plt.close()
#
# item=214
# get_fig_Z2(str(int(item)) + '.txt')

cand_ID=fits_to_txt.cand_ID
for item in cand_ID[499:]:
    get_fig_Z2(str(int(item)) + '.txt')
# Z2_p=Z2-signal_Z2
# p_possible=1.0/freq[np.where(Z2_p>0)]
# print p_possible
# np.savetxt(path+dataname[3:7]+"/p_possible.txt",p_possible)0
