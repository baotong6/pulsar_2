# import numpy as np
# import matplotlib.pyplot as plt
# from astropy.io import fits
# import sys
# import os
# import string
# import datetime
# from scipy.interpolate import lagrange
# from scipy import interpolate
# from scipy.fftpack import fft,ifft
# import scipy.signal as ss
# import scipy.stats as stats
# import pandas as pd
# path='/Users/baotong/Desktop/period/txt/'
# time=np.loadtxt(path+'pwn.txt')[:,0]
# sep=time[1:-1]-time[0:-2]
# print(sep)
#
# plt.hist(sep,bins=np.linspace(0,500,1000))
# plt.show()


import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from pylab import *
matplotlib.rcParams['axes.unicode_minus'] =False
matplotlib.rcParams['font.sans-serif'] = ['SimHei']
v = np.linspace(-1,10,1000)
y=np.zeros(len(v))
for i in range(len(v)):
    if v[i]<=0:
        y[i]=0
    else:
        y[i]=np.e**(-v[i])
print(y)
plt.plot(v,y)
plt.show()