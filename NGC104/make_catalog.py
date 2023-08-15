#!/bin/bash
# -*- coding: utf-8 -*-
# written by Tong
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
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.timeseries import LombScargle
import stingray
from stingray import Lightcurve, Powerspectrum, AveragedPowerspectrum
from stingray.events import EventList
import hawkeye as hawk

pos_all={'Tuc':[6.0236250, -72.0812833, 3.17 * 60, 3.17 / 8.8 * 60],
         'terzan5':[267.0202083, -24.7790556, 0.72 * 60, 0.72 / 3.4 * 60],
         'M28':[276.1363750, -24.8702972, 1.97 * 60, 1.97 / 8.2 * 60],
         'omg':[201.69700, -47.47947, 5 * 60, 5 / 2.1 * 60],
         'NGC6397':[265.17539, -53.67433, 2.9 * 60, 2.9 / 58 * 60],
         'NGC6752':[287.71713, -59.98455, 1.91* 60, 1.91 / 11.24 * 60],
         'NGC6266':[255.303333,-30.113722,0.92* 60,0.92/4.2*60],
         'M30':[325.092167,-23.179861,1.03* 60,1.03/17.2*60],
         'NGC6121':[245.89675,-26.52575,4.33*60, 1.17*60],
         'NGC6304':[258.63438, -29.46203,1.42*60,0.21*60],
         'NGC6656':[279.09975, -23.90475,3.36*60, 1.34*60]}

# gcname=['Tuc','terzan5','omg','M28','NGC6397','NGC6752','NGC6266','M30','NGC6121','NGC6304','NGC6656']
gcname=['NGC6121','NGC6304','NGC6656']
srcnum=[592,489,300,502,376,244,146,84]
# catname=['xray_properties-592.fits','cheng2019_terzan.fit','cheng2020_omg.fit',
#          'cheng2020_M28.fit','ngc6397_catalog.fits','ngc6752_catalog.fits',
#          'NGC6266_p50_i5_src_1_2_4_8.fits','M30_p50_i5_src_1_2_4_8.fits',
#          'NGC6121_p50_i5_src_1_2_4_8.fits','NGC6304_p50_i5_src_1_2_4_8.fits','NGC6656_p50_i5_src_1_2_4_8.fits']
catname=['NGC6121_p50_i5_src_1_2_4_8.fits','NGC6304_p50_i5_src_1_2_4_8.fits','NGC6656_p50_i5_src_1_2_4_8.fits']
def extract_info(path,catname,keyword):
    catalog_info = fits.open(path + catname)
    ra=catalog_info[1].data[keyword[0]]
    dec=catalog_info[1].data[keyword[1]]
    print(keyword[2])
    [ra_center,dec_center,rhl,rc]=pos_all[keyword[2]]
    c1 = SkyCoord(ra * u.deg, dec * u.deg, frame='fk5')
    c2 = SkyCoord(ra_center * u.deg, dec_center * u.deg, frame='fk5')
    dist_all = c1.separation(c2)
    dist_all=dist_all.arcsec

    seq=np.arange(1,len(ra)+1,1)
    path_txt=path+'txt_all_obs_p90/'
    counts_all=[];exptime_all=[];VI_all=[]
    for i in range(len(seq)):
        src_evt=np.loadtxt(path_txt+str(seq[i])+'.txt')
        src_epoch=np.loadtxt(path_txt+'epoch_src_'+str(seq[i])+'.txt')
        if src_epoch.ndim == 1: src_epoch = np.array([src_epoch])
        CR = hawk.plot_longT_V(src_evt=src_evt, bkg_file=None, epoch_info=src_epoch,show=False)
        VI_all.append(np.max(CR)/np.min(CR))
        counts_all.append(len(src_evt));exptime_all.append(np.sum(src_epoch[:,3]))
    src_info=np.column_stack((seq,ra,dec,counts_all,exptime_all,VI_all,dist_all))
    np.savetxt(path+'src_info.txt',src_info,fmt='%10d %10.5f %10.5f %10d %15d %10.5f %10.5f')

if __name__=='__main__':
    for i in range(len(gcname)):
        print(i)
        # if i==1 or i==2 or i==3:keyword = ['RAJ2000', 'DEJ2000',gcname[i]]
        # elif i==0 or i==4 or i==5:keyword=['RAdeg','DEdeg',gcname[i]]
        # else:keyword=['RA','DEC',gcname[i]]
        keyword = ['RA', 'DEC', gcname[i]]
        extract_info(path='/Users/baotong/Desktop/period_' + gcname[i] + '/', catname=catname[i],
                     keyword=keyword)
    # plot_somefig()
