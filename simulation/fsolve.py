# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 19:21:40 2019

@author: H.S.Wang
"""
# delta = 0.1
# w = 1.0/5000
# t = 2500
#
#
# def f(t1):
#     return t1-delta/w*np.cos(w*t1)
#
#
# result = fsolve(f,2800)

import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp
import functools
import datetime
from astropy.io import fits
import sys
import os
import string
#import correct as correct
import random
from scipy import optimize
from sympy import *
from scipy.optimize import fsolve
from itertools import chain
def get_O1m(Om1w,w_range):
  O1m = np.zeros(12)

  w = w_range
  w_hi = np.max(w_range)
  w_lo = np.min(w_range)
  pw = 1 / w / np.log(w_hi / w_lo)
  for i in range(12):
    O1m[i] = np.trapz(pw * np.random.random(len(w))*1e-2, w)
  O_period = np.sum(O1m[1:])
  return O_period
a=[]
b=[]
c=[]
for i in range(1000):
  a.append(get_O1m(np.random.random(12)*1e-2,w_range=np.arange(1./50000.,1./10000.,1e-9)))
  b.append(get_O1m(np.random.random(12)*1e-2,w_range=np.arange(1./10000.,1./3000.,1e-8)))
  c.append(get_O1m(np.random.random(12)*1e-2,w_range=np.arange(1./3000.,1./300.,1e-7)))

print(np.mean(a))
print(np.mean(b))
print(np.mean(c))
#
# data_file='data.txt'
# p_test=156.9
# bin=15
# shift=0.1
#
# time=np.loadtxt(data_file)[:,1]
# def trans(t,p_test,shift):
#   ti =t
#   v = 1.0 /p_test
#   turns = v * ti
#   turns += shift
#   # 初始相位
#   for i in range(len(turns)):
#     turns[i] = turns[i] - int(turns[i])
#   return turns
#
# turns=trans(time,p_test,shift)
# loc=np.zeros(bin)
# for index in turns:
#   loc[int(index*bin)] += 1
#
# x = np.array([(i / bin + 0.5 / bin) for i in range(bin)])
#
# x2=np.concatenate((x,x+1))
# y2=np.concatenate((loc,loc))
#
# plt.figure(1,(8,8))
#
# plt.xlabel('phase')
# plt.ylabel('counts/bin')
# plt.step(x2,y2,color='red')
# # plt.savefig('phase' + str(bin) + '.png')
# plt.show()