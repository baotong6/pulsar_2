import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import string
from datetime import datetime
from scipy.interpolate import lagrange
from scipy import interpolate
from scipy.optimize import curve_fit
import pandas as pd
import linecache

path='/Users/baotong/Desktop/HST/phot/'
phot=np.loadtxt(path+'F373N.phot')
ra=phot[:,2]
dec=phot[:,3]
SN=phot[:,5]
with open(path+'region_all.reg','a+') as f:
    k=0
    while k<len(ra):
        if SN[k]>5:
            line='circle({0},{1},4)'.format(ra[k],dec[k])+'\n'
            f.writelines(line)
        k+=1
f.close()