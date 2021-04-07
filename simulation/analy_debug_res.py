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

sim_N=100
path='/Users/baotong/Desktop/period_LW/simulation/test/'
P1=50000.;P2=5000.;P3=500.
w1 = 2 * np.pi * np.arange(1. / 55000, 1. / 45000., 1e-9)
w2 = 2 * np.pi * np.arange(1. / 5500, 1. / 4500., 1e-8)
w3 = 2 * np.pi * np.arange(1. / 550, 1. / 450., 1e-7)

def get_info(P):
    mean=[]
    out=0
    for k in range(sim_N):
        temp=np.loadtxt(path+'res_{0}/res_{1}.txt'.format(int(P),k))
        out+=temp
    return out

w_s1=get_info(P1)
w_s2=get_info(P2)
w_s3=get_info(P3)
out_w1=np.zeros(len(w_s1[0]))
out_w2=np.zeros(len(w_s2[0]))
out_w3=np.zeros(len(w_s3[0]))
# w_s1=[np.sum(w_s1[0]) for i in range(len(w_s1))]
# w_s2=[np.sum(w_s2) for i in range(len(w_s2))]
# w_s3=[np.sum(w_s3) for i in range(len(w_s3))]
for i in range(len(w_s1[0])):
    out_w1[i]=np.sum(w_s1[:,i])
for i in range(len(w_s2[0])):
    out_w2[i]=np.sum(w_s2[:,i])
for i in range(len(w_s3[0])):
    out_w3[i]=np.sum(w_s3[:,i])

plt.semilogx()
# plt.semilogy()
# plt.ylim(1e-4,1e4)
plt.plot(w1,out_w1)
plt.plot(w2,out_w2)
plt.plot(w3,out_w3)
plt.legend(['P=50000s','P=5000s','P=500s'])
plt.show()
