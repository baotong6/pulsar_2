#!/bin/bash
# -*- coding: utf-8 -*-
#pulsation=0.0331

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import string
import correct as correct
#import lc as lc
from datetime import datetime
path='/Users/baotong/desktop/pulsar/xpnavdata/'
TB_new=np.load("right_time.npy")
t2=np.load("obs_time.npy")
#print TB-t2
#print t2[-1]-t2[0]
epoch=321321600+0.007781
#2018年03月08到2008年1月1日之间的秒数
p_test=1/29.6328827000
#2018年03月10日的数据
vary=-4.
v=1/(p_test+vary*1.e-8)
vdot=-369161.97*1e-15
#根据表来查
border=100
#print t2[-1]-t2[0]
#print t2
TB=[]
for index in TB_new:
    TB.append(index[0][0])
t2=np.array(t2)
TB=np.array(TB)
#print t2-TB

def trans(t):
    ti=t-epoch
    #print ti
    vary=[i for i in range(-border,border)]
    vary=np.array(vary)
    p_test = 1 / 29.6471534323
    v = 1 / (p_test + vary * 1.e-8)
    p=1.0/v
    pdot=-vdot/(v*v)
    vddot = 2.0 * pdot * pdot / (p * p * p)
    freq = v + ti * vdot + ti * ti / 2.0 * vddot
    preq = 1.0 / freq
    turns = v * ti + vdot * ti * ti / 2.0 + vddot * ti * ti * ti / 6.0
    shift=0
    #turns=list(turns)
    #turns=np.array(turns[0])
    #初始相位
    for i in range(len(turns)):
        turns[i]=turns[i]-shift - int(turns[i]-shift)
    return turns
#print TB
#print trans(TB[0])
turns_n=[]
for index in TB:
    turns_n.append(trans(index))
turns_n=np.transpose(turns_n)
#turns_n=[border,length of time]
bin=100
#print len(turns_n[0])
sort_turn=[]
for index in turns_n:
    sort_turn.append(np.sort(index))
#sort_turn=[border,length of time]
loc=np.empty([border,bin])
#print loc
#print sort_turn[0]
for i in range(border):
    for index in sort_turn[i]:
        loc[i][int(index / (1.0 / bin))] += 1
#print sum(loc[0])
standard=[sum(loc[0])/float(bin) for i in range(bin)]
chi=[]
for index in loc:
    #print index-standard
    chi.append((index-standard).dot(index-standard))
print chi.index(max(chi))

loc=loc*bin/((t2[-1]-t2[0]))
#loc=loc*bin
x=[(1.0*i/bin+0.5/bin) for i in range(bin)]
#print len(x)
#print(len((x)))
#print(len(loc))
print('done')
#print loc
#plt.plot(x,loc[chi.index(max(chi))])
#plt.plot(x,loc[border/2-1])
plt.scatter(x,loc[chi.index(max(chi))],color='blue',marker='+')
plt.scatter(x,loc[border/2-1],color='green',marker='+')

plt.legend(["max","original"])
#原始数据
plt.show()



