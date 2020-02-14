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
from scipy.fftpack import fft,ifft
import scipy.signal as ss
import random
from astropy.wcs import WCS
path='/Users/baotong/Desktop/period_gc/'
p90_list='mod_expweighted_90_mean_i4.psfmap'
srcname_list=path+'combineobs_info_box.txt'

phy_x=np.loadtxt(srcname_list)[:,2]
phy_y=np.loadtxt(srcname_list)[:,3]

phy_x=np.rint(phy_x)
phy_y=np.rint(phy_y)
phy_x_int=phy_x.astype(np.int)
phy_y_int=phy_y.astype(np.int)

src_x=phy_x_int-2896
src_y=phy_y_int-2896
#将physical坐标转换为img坐标
hdul_p90=fits.open(path+p90_list)

p90_data=hdul_p90[0].data
p90_data=p90_data.T
src_radius=p90_data[src_x,src_y]
src_radius*=2.03252
#src_radius*=2
#pixel size=0.492arcsec
#the default units of p90 fits is arcsec, 1/0.492=2.03252
print(np.sort(src_radius))
for i in range(len(phy_x)):
    with open(path+'region/{0}.reg'.format(i+1),'w+') as f1:
        reg='circle('+str(phy_x[i])+','+str(phy_y[i])+','+str(src_radius[i])+')'
        f1.writelines(reg)
    with open(path+'region/all.reg','a+') as f2:
        f2.writelines(reg+'\n')

