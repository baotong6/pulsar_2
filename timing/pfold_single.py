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

#path='/Volumes/halo/chandra_obs_hr_m04_for_baotong/'
#path='/Users/baotong/xmm_obs/'
path='/Volumes/pulsar/xmm_obs_16/'
#path='/Volumes/halo/redgiant/'
# path='/Volumes/pulsar/xmm_obs_16/0112870101/'
# path='/Volumes/halo/chandra_obs_hr_m04_for_baotong/'
path='/Users/baotong/Desktop/period_LW/xmm_CV/0111310101/'
#dataname='src22'
dataname='HT_CAS_pn'
time = np.loadtxt(path + 'txt/' + dataname + '.txt')[:, 0]
#time = np.loadtxt(path + 'txt/' + dataname + '_s.txt')[:, 0]
#time=np.loadtxt(path + 'txt/evt_list10956.txt')[:,0]

N = len(time)

vdot=0
#根据表来查
num_p=1.0
# p_test=49405.45
def trans(t,num_p,p_test,resolution,border):
    ti=t
    #print ti
    vary=[i for i in range(-border,border)]
    vary=np.array(vary)
    #p_test = 1.0/5.500550055005501e-06
    #p_test=4945.45
    p_test=num_p*p_test
    v = 1.0 / (p_test + vary *resolution)
    p=1.0/v
    # pdot=-vdot/(v*v)
    # vddot = 2.0 * pdot * pdot / (p * p * p)
    freq = v # + ti * vdot + ti * ti / 2.0 * vddot
    preq = 1.0 / freq
    turns = v * ti #+ vdot * ti * ti / 2.0 + vddot * ti * ti * ti / 6.0
    shift=0.1
    turns += shift
    #初始相位
    for i in range(len(turns)):
        turns[i]=turns[i] - int(turns[i])
    turns = num_p * turns
    return turns

def getIndexes(y_predict, y_data):
    n = y_data.size
    # SSE为和方差
    SSE = ((y_data - y_predict) ** 2).sum()
    # MSE为均方差
    MSE = SSE / n
    # RMSE为均方根,越接近0，拟合效果越好
    RMSE = np.sqrt(MSE)

    # 求R方，0<=R<=1，越靠近1,拟合效果越好
    u = y_data.mean()
    SST = ((y_data - u) ** 2).sum()
    SSR = SST - SSE
    R_square = SSR / SST
    return SSE, MSE, RMSE, R_square

def func(x,a,b,c):
    return a*np.sin(2*np.pi*x+b)+c

border=100
resolution=0.1
def phasefold(time,num_p,bin,p_test):
    turns_n=[]
    for index in time:
        turns_n.append(trans(index,num_p,p_test,border=100,resolution=0.1))
    turns_n=np.transpose(turns_n)
    #bin=20
    #print len(turns_n[0])
    sort_turn=[]
    for index in turns_n:
        sort_turn.append(np.sort(index))
    #print len(sort_turn)
    loc=np.zeros([2*border,bin])
    for i in range(2*border):
        for index in sort_turn[i]:
            loc[i][int(index / (num_p / bin))] += 1

    x = np.array([(num_p * i / bin + 0.5 / bin) for i in range(bin)])
    return [x,loc]
# print sum(loc[0])
# standard=[sum(loc[0])/float(bin) for i in range(bin)]
# chi=[]
#
# for index in loc:
#     #print index-standard
#     chi.append((index-standard).dot(index-standard))
# print chi.index(max(chi))

# 利用curve_fit作简单的拟合，popt为拟合得到的参数,pcov是参数的协方差矩阵

def error_sin(x,loc):
    error_MSE=[]
    cof=[]
    for i in range(len(loc)):
        ydata=loc[i]
        popt, pcov = curve_fit(func, x, ydata,bounds=(0, [1000., 2*np.pi, 3000.]))
        y_predict = func(x, *popt)
        indexes_1 = getIndexes(y_predict, ydata)
        error_MSE.append(indexes_1[1])
        cof.append(popt)
    return [error_MSE,cof]
p_test=6363.1008
x=phasefold(time,num_p=1.0,bin=20,p_test=p_test)[0]
loc=phasefold(time,num_p=1.0,bin=20,p_test=p_test)[1]
error_MSE = error_sin(x, loc)[0]
cof = error_sin(x, loc)[1]
print('done')
# print len(loc)
plt.figure(1)
plt.step(x, loc[error_MSE.index(min(error_MSE))])
plt.plot(x, func(x, *cof[error_MSE.index(min(error_MSE))]))
plt.step(x, loc[border])
plt.plot(x, func(x, *cof[border]))
plt.legend(["min","min_sin","original","original_sin"])
#print p_test+(error_MSE.index(min(error_MSE))-border)*resolution/num_p
print(error_MSE.index(min(error_MSE)))
plt.figure(2)
x1=x+1
x2=np.concatenate((x,x1))
y2=np.concatenate((loc[border],loc[border]))
plt.step(x2,y2)
plt.plot(x2,func(x2,*cof[border]))
per=p_test
plt.text(0.5,max(loc[border]),"p="+str(per))
plt.show()

