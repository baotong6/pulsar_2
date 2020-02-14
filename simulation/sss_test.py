import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import string
import datetime
from scipy.interpolate import lagrange
from scipy import interpolate
from scipy.fftpack import fft,ifft
import scipy.signal as ss
import scipy.stats as stats
import pandas as pd
path='/Users/baotong/Desktop/period/txt/'
time=np.loadtxt(path+'pwn.txt')[:,0]
sep=time[1:-1]-time[0:-2]
print(sep)

plt.hist(sep,bins=np.linspace(0,500,1000))
plt.show()
