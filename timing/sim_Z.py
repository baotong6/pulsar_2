#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import pylab as pl
import string
from datetime import datetime
from scipy.interpolate import lagrange
from scipy import interpolate
from scipy.fftpack import fft,ifft
import scipy.signal as ss
import random
import pandas as pd
import timing.read_data as data
import scipy.stats as stats
import matplotlib.pyplot as plt
import matplotlib.style as style
from IPython.core.display import HTML

def sim_T():
#包含freq分布

    path = '/Users/baotong/Desktop/period/'
    dataname='pwn.txt'
    time=data.get_data(dataname)[0]
    energy=data.get_data(dataname)[1]
    obs_ID=data.get_data(dataname)[2]
    t=data.get_data(dataname)[3]
    E=data.get_data(dataname)[4]
    dict=data.get_data(dataname)[5]

    Z2=[]
    TURNS=[]
    border = 5000
    vary = np.array([i for i in range(0, border)])
    freq = 1/30000. + vary * 1.e-5
    print(freq)

    # size = 1000
    # sep = np.random.poisson(182, size)
    # t0 = 100 * np.random.random(1)[0]
    # t = [t0]
    # for i in range(1, size):
    #     t.append(t[i - 1] + sep[i - 1])
    t = t[9]

    size = 1000

    sep2 = np.random.poisson(123, size)
    t02 = 100 * np.random.random(1)[0]
    t2 = [t02]
    for i in range(1, size):
        t2.append(t2[i - 1] + sep2[i - 1])

    t2 = np.array(t2)

    #t=np.concatenate((t,t2))

    for i in range(border):
        #t += 10000000.
        w=freq[i]
        turns = w* t
        INT_turns = np.trunc(turns)
        turns = turns - INT_turns
        turns=2*np.pi*turns
        TURNS.append(turns)

        Z2.append((2.0/size)*(sum(np.cos(turns))**2+sum(np.sin(turns))**2))
    print(np.sort(Z2))
    print(max(Z2))

    plt.hist(Z2,bins=20,histtype='step',color='green')

    print(np.where(Z2==max(Z2))[0])

    #plt.hist(TURNS[np.where(Z2==max(Z2))[0][0]], bins=20, histtype='step', color='green')
    # plt.semilogx()
    plt.step(freq,Z2)
    plt.show()

#sim_T()
def dis_Z2_pure_turns():
    size1=100
    size2=5000
    border=5000
    Z2_rand_1= []
    COS_rand_1= []
    SIN_rand_1= []
    TURNS_rand_1= []
    Z2_rand_2= []
    COS_rand_2= []
    SIN_rand_2= []
    TURNS_rand_2= []

    for i in range(border):
        turns_rand_1 = 2 * np.pi * np.random.random(size1)
        temp1_rand_1 = sum(np.cos(turns_rand_1)) * (2.0 / size1) ** 0.5
        temp2_rand_1 = sum(np.sin(turns_rand_1)) * (2.0 / size1) ** 0.5
        TURNS_rand_1.append(turns_rand_1)
        Z2_rand_1.append((temp1_rand_1 ** 2 + temp2_rand_1** 2))

        turns_rand_2 = 2 * np.pi * np.random.random(size2)
        temp1_rand_2 = sum(np.cos(turns_rand_2)) * (2.0 / size2) ** 0.5
        temp2_rand_2 = sum(np.sin(turns_rand_2)) * (2.0 / size2) ** 0.5
        TURNS_rand_2.append(turns_rand_2)
        Z2_rand_2.append((temp1_rand_2 ** 2 + temp2_rand_2** 2))
    weights1=np.ones_like(Z2_rand_1) / float(len(Z2_rand_1))
    weights2 = np.ones_like(Z2_rand_2) / float(len(Z2_rand_2))
    plt.hist(Z2_rand_1, bins=50, weights=weights1, histtype='step', color='blue')
    plt.hist(Z2_rand_2, bins=50, weights=weights2, histtype='step', color='green')
    plt.show()

dis_Z2_pure_turns()

def dis_Z2():
    Z2=[]
    Z2_rand=[]
    COS=[]
    SIN=[]
    COS_rand=[]
    SIN_rand=[]
    border=1000
    TURNS=[]
    TURNS_rand=[]
    for i in range(border):
        size=1000
        sep=np.random.poisson(182,size)
        t0=100*np.random.random(1)[0]
        t=[t0]
        for i in range(1,size):
            t.append(t[i-1]+sep[i-1])

        t=np.array(t)
        #t += 10000000.
        w=1/283.8

        turns = w* t
        INT_turns = np.trunc(turns)
        turns = turns - INT_turns
        turns=2*np.pi*turns

        turns_rand=2*np.pi*np.random.random(size)

        temp1=sum(np.cos(turns))*(2.0/size)**0.5
        temp2=sum(np.sin(turns))*(2.0/size)**0.5

        temp1_rand=sum(np.cos(turns_rand))*(2.0/size)**0.5
        temp2_rand=sum(np.sin(turns_rand))*(2.0/size)**0.5

        #print(sum(np.cos(turns)))
        #print(sum(np.cos(turns_rand)))
        TURNS.append(turns)
        TURNS_rand.append(turns_rand)

        Z2.append((temp1**2+temp2**2))
        Z2_rand.append((temp1_rand**2+temp2_rand**2))
        COS.append(temp1)
        SIN.append(temp2)
        COS_rand.append(temp1_rand)
        SIN_rand.append(temp2_rand)
    print(len(Z2))
    print(len(Z2_rand))
    #plt.hist(Z2,bins=20,histtype='step',color='green')

    weights1 = np.ones_like(Z2_rand) / float(len(Z2_rand))
    plt.plot(np.linspace(0, 40, 1000), stats.chi2.pdf(np.linspace(0, 40, 1000), df=2))
    plt.hist(Z2_rand, bins=50, weights=weights1,histtype='step', color='blue')
    #plt.hist(COS, bins=20, histtype='step', color='red')
    #plt.hist(SIN, bins=20, histtype='step', color='blue')
    #plt.hist(TURNS[1],bins=100,histtype='step', color='purple')
    #plt.hist(TURNS_rand[1],bins=100,histtype='step', color='grey')
    print(sum(np.cos(TURNS[1])))
    print(sum(np.cos(TURNS_rand[1])))
    print(Z2[1])
    print(Z2_rand[1])
    plt.show()
#dis_Z2()

def plot_turns_dis():
    A=[]
    size=1000
    for i in range(1000):
        a=np.random.random(size)
        a=2*np.pi*a
        A.append((2.0/size)**0.5*sum(np.sin(a)))

    A=np.array(A)
    plt.hist(A,100,histtype='step')
    plt.show()


def flatten(a):
    if not isinstance(a, (list, )):
        return [a]
    else:
        b = []
        for item in a:
            b += flatten(item)
    return b

def plot_pwn_sep():
    path = '/Users/baotong/Desktop/period/'
    dataname='pwn.txt'
    time=data.get_data(dataname)[0]
    energy=data.get_data(dataname)[1]
    obs_ID=data.get_data(dataname)[2]
    t=data.get_data(dataname)[3]
    E=data.get_data(dataname)[4]
    dict=data.get_data(dataname)[5]

    t=np.array(t)

    sep=[]
    for i in range(len(t)):
        sep.append((t[i][1:]-t[i][:-1]))


    sep=np.array(sep)
    sep_time=time[1:]-time[:-1]

    a=sep[0]
    for i in range(1,i):
        a=np.concatenate((a,sep[i]))
    flux=np.zeros(len(t))
    for i in range(len(t)):
        flux[i]=(t[i][-1]-t[i][0])/len(t[i])
        flux[i] =np.sort(a)[int(len(a)/2)]


    # plt.plot(np.linspace(0, 20, 100), stats.chi2.cdf(np.linspace(0, 20, 100), df=2))  # 绘制累积概率密度函数
    #
    # weights1 = np.ones_like(Z2) / float(len(Z2))
    # plt.subplot(211)
    # plt.hist(Z2,bins=200,weights=weights1, histtype='step')
    # plt.plot(np.linspace(0, 40, 1000), stats.chi2.pdf(np.linspace(0, 40, 1000), df=2))  # 绘制0到20的卡方分布曲线,给定自由度为2
    # plt.fill_between(np.linspace(0, 40, 1000), stats.chi2.pdf(np.linspace(0, 40, 1000), df=2), alpha=0.15)  # 填充曲线
    # plt.subplot(212)
    # plt.hist(Z2_standard,bins=200,weights=weights1,histtype='step')
    # plt.plot(np.linspace(0, 40, 1000), stats.chi2.pdf(np.linspace(0, 40, 1000), df=2))  # 绘制0到20的卡方分布曲线,给定自由度为2
    # plt.fill_between(np.linspace(0, 40, 1000), stats.chi2.pdf(np.linspace(0, 40, 1000), df=2), alpha=0.15)  # 填充曲线
    # #plt.step(a[1],a[0])
    weights1=np.ones_like(a)/float(len(a))
    #plt.hist(a,bins=100,weights=weights1,histtype='step')

    plt.hist(a[1284:2296], bins=50, histtype='step')
    plt.hist(np.random.poisson(flux[9], len(a[1284:2296])), bins=50, histtype='step')
    print(sum(a[1284:2296]))
    print(sum(np.random.poisson(flux[9],len(a[1284:2296]))))


    # for i in range(len(flux)):
    #     plt.hist(np.random.poisson(flux[i],len(a)),bins=100,histtype='step')

    #plt.plot(np.linspace(0, 40, 1000), stats.chi2.pdf(np.linspace(0, 40, 1000), df=2))  # 绘制0到20的卡方分布曲线,给定自由度为2
    plt.show()
#plot_pwn_sep()

def plot_guess():
    plt.plot(np.linspace(0, 40, 1000), stats.chi2.pdf(np.linspace(0, 40, 1000), df=2))
    print(stats.chi2.pdf(np.linspace(0, 40, 1000),df=2))
    plt.show()

#plot_guess()

#print(np.sort(x))
# pillar = 150
# a = plt.hist(x, bins=pillar,color='g', alpha=0.5)
#plt.plot(a[1][0:pillar], a[0], 'r')
# plt.grid()
# plt.show()
# for i in range(size-1):
#     for j in range(i)