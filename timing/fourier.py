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
import random
#import seaborn

time=data.time
energy=data.energy
obs_ID=data.obs_ID
t=data.t
E=data.E
# t_bin=0.5
# time_col=[]
# sort_time=np.sort(time)
# print sort_time
# temp=time[0]
# while temp<sort_time[-1]:
#     time_col.append(((temp<=sort_time) & (sort_time<temp+t_bin)).sum())
#     temp=temp+t_bin
# print np.sort(time_col)
# x=np.arange(time[0],temp,t_bin)
# print len(x)
# print len(time_col)
# cts_rate=np.array(time_col)/t_bin
# plt.hist(cts_rate,len(x))
# plt.hist(t[0], bins=int((t[0][-1]-t[0][0])/t_bin), normed=0, facecolor='black', edgecolor='black',alpha=1,histtype='step')
# plt.show()

pd=3.76
frq=2*np.pi/pd
a=np.concatenate((t[0],t[1]))
#a=time
b=np.arange(0,len(a),1)
k=len(a)/(a[-1]-a[0])
c=k*(a-a[0])
plt.plot(a,b)
plt.plot(a,c)
y=b-c
plt.plot(a,y)
plt.show()
#z=np.arange(5.e-6,2.e-3,1.e-6)
z=np.arange(0.1,20.,1.e-2)
w=ss.lombscargle(a,y,z)
# plt.subplot(211)
# plt.plot(z,w)
p_value=z-z
epoch=100
for i in range(epoch):
    b1=b;a1=a
    random.shuffle(b1)
    k1= len(a)/(a[-1] - a[0])
    c1= k*(a-a[0])
    y1=b1-c1
    #z1= np.arange(5.e-6, 2.e-3, 1.e-6)
    z1= np.arange(0.1, 20.0, 1.e-2)
    w1= ss.lombscargle(a1, y1, z1)
    plt.subplot(221)
    plt.plot(a, y)
    plt.subplot(222)
    plt.plot(a1, y1)
    plt.subplot(223)
    plt.plot(z, w)
    plt.subplot(224)
    plt.plot(z1, w1)

    plt.show()
    l=w1-w
    l[l>0]=1
    l[l<=0]=0
    p_value=p_value+l
p_value/=epoch
# plt.subplot(212)
# plt.plot(z,p_value)
# plt.show()


# yy=fft(y)                     #快速傅里叶变换
# yreal = yy.real               # 获取实数部分
# yimag = yy.imag               # 获取虚数部分
#
# yf=abs(fft(y))                # 取绝对值
# yf1=abs(fft(y))/len(x)           #归一化处理
# yf2 = yf1[range(int(len(x)/2))]  #由于对称性，只取一半区间
#
# xf = np.arange(len(y))        # 频率
# xf1 = xf
# xf2 = xf[range(int(len(x)/2))]  #取一半区间
#
# print yy
# plt.subplot(221)
# plt.plot(x,y)
# plt.title('Original wave')
#
# plt.subplot(222)
# plt.plot(xf,yf,'r')
# plt.title('FFT of Mixed wave(two sides frequency range)',fontsize=7,color='#7A378B')  #注意这里的颜色可以查询颜色代码表
#
# plt.subplot(223)
# plt.plot(xf1,yf1,'g')
# plt.title('FFT of Mixed wave(normalization)',fontsize=9,color='r')
#
# plt.subplot(224)
# plt.plot(xf2,yf2,'b')
# plt.title('FFT of Mixed wave)',fontsize=10,color='#F08080')
#
#
# plt.show()