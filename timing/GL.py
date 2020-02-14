##Gregory-Loredo method  written by baotong 2019-09-29
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
import datetime
from scipy.interpolate import lagrange
from scipy import interpolate
import read_data as data
from scipy.fftpack import fft,ifft
import scipy.signal as ss
import scipy.stats as stats
import random
import scipy.stats as stats
import matplotlib.pyplot as plt
import matplotlib.style as style
from IPython.core.display import HTML
import read_data as data

path='/Users/baotong/Desktop/period/'
dataname='1266.txt'
def compute_bin(Tlist,m,w,fi):
    time = data.get_data(dataname)[0]
    energy = data.get_data(dataname)[1]
    obs_ID = data.get_data(dataname)[2]
    t = data.get_data(dataname)[3]
    E = data.get_data(dataname)[4]
    dict = data.get_data(dataname)[5]

time=data.get_data(dataname)[0]
energy=data.get_data(dataname)[1]
obs_ID=data.get_data(dataname)[2]
t=data.get_data(dataname)[3]
E=data.get_data(dataname)[4]
dict=data.get_data(dataname)[5]

y=compute_bin(time,30,2*np.pi/5409.6,0.0)
x=np.arange(1,31)
print(x)
plt.step(x,y/(sum(y)/30.))
plt.plot([0,31],[1.0,1.0],'--')
plt.show()
