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
#path='/Users/baotong/Desktop/period/'
path='/Users/baotong/Desktop/period_LW/'
path='/Volumes/pulsar/LimWin_damage/merge_data/timing/'
#evt_list='LW_merged_evt.fits'
#p90_list='SgrA_grt_2000_8000_psf90.fits'
#p90_list='SgrA_2000_8000_psf90.fits'
p90_list='expweighted_90_mean_e1.psfmap'
#srcname_list='LimWin_p50_2000_8000_src.fits'
fitsname='LW_e1.fits'
w=WCS(path+fitsname)
# lon, lat = w.all_world2pix(266.3980091,-29.02585626,1)

source_info = np.loadtxt(path+'catalog_LW.txt')
ra = source_info[:,1]
dec = source_info[:,2]

#print(ra)
#print(dec)

src_x,src_y=w.all_world2pix(ra,dec,1)
print(src_x,src_y)
src_x=np.rint(src_x)
src_y=np.rint(src_y)
src_x_int=src_x.astype(np.int)
src_y_int=src_y.astype(np.int)
phy_x=src_x+2896
phy_y=src_y+2896
# phy_x=src_x+3584
# phy_y=src_y+3584
#在LW中

#将img坐标转换为physical坐标

hdul_p90=fits.open(path+p90_list)

p90_data=hdul_p90[0].data
p90_data=p90_data.T
src_radius=p90_data[src_x_int,src_y_int]
src_radius*=2.03252
#src_radius*=2
#pixel size=0.492arcsec
#the default units of p90 fits is arcsec, 1/0.492=2.03252
print(np.sort(src_radius))
for i in range(len(phy_x)):
    with open(path+'region_90/{0}.reg'.format(i+1),'w+') as f1:
        reg='circle('+str(phy_x[i])+','+str(phy_y[i])+','+str(src_radius[i])+')'
        f1.writelines(reg)
    with open(path+'region_90/all.reg','a+') as f2:
        f2.writelines(reg+'\n')
