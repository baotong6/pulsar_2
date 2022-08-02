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
         'M30':[325.092167,-23.179861,1.03* 60,1.03/17.2*60]}

gcname=['Tuc','terzan5','omg','M28','NGC6397','NGC6752','NGC6266','M30']
srcnum=[537,489,300,502,376,244,146,84]
catname=['cheng2019_Tuc.fit','cheng2019_terzan.fit','cheng2020_omg.fit',
         'cheng2020_M28.fit','ngc6397_catalog.fits','ngc6752_catalog.fits',
         'NGC6266_p50_i5_src_1_2_4_8.fits','M30_p50_i5_src_1_2_4_8.fits']

def choose_ecf(seq,ra,dec,ecf=[50,75,90],rad_ecf=None):
    close_dist=[];ecf_all=[]
    for i in range(len(seq)):
        c1=SkyCoord(ra[i] * u.deg, dec[i] * u.deg, frame='fk5')
        c2=SkyCoord(ra * u.deg, dec * u.deg, frame='fk5')
        dist_all=c1.separation(c2)
        dist_all = dist_all.arcsec
        close_dist.append(np.sort(dist_all)[1])
        j=0
        while j < len(ecf):
            if close_dist[i]<rad_ecf[j][i]:
                ecf_all.append(ecf[j])
                continue
            else:
                j+=1

    return (ecf_all,close_dist)

# def read_psf():
#     path='/Volumes/pulsar/47Tuc/merge_data/timing/'
#

def extract_info(path,catname,keyword,epoch_name,bkg=0):
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
    counts_all=[];exptime_all=[];VI_all=[];bkg_counts=[]
    for i in range(len(seq)):
        src_evt=np.loadtxt(path_txt+str(seq[i])+'.txt')
        bkg_evt=np.loadtxt(path_txt+str(seq[i])+'_bkg.txt')
        src_epoch=np.loadtxt(path_txt+'epoch_src_'+str(seq[i])+'.txt')
        if src_epoch.ndim == 1: src_epoch = np.array([src_epoch])
        CR = hawk.plot_longT_V(src_evt=src_evt, bkg_file=None, epoch_info=src_epoch,show=False)
        VI_all.append(np.max(CR)/np.min(CR))
        counts_all.append(len(src_evt));exptime_all.append(np.sum(src_epoch[:,3]));bkg_counts.append(len(bkg_evt))
    if bkg:
        epoch_info = np.loadtxt(path + f'{epoch_name}.txt')
        path_bkgscale=path+'region_startover/exp/'
        path_alltxt=path+'txt_startover/'
        obsidlist=epoch_info[:,2]
        bkg_cts=np.zeros(len(seq));bkg_scale=np.zeros(len(seq))
        for obsid in obsidlist:
            obsid=int(obsid)
            bkg_pix=np.loadtxt(path_bkgscale+f'{obsid}_bkg_psf90.txt')[:,3]
            src_pix=np.loadtxt(path_bkgscale+f'{obsid}_psf90.txt')[:,3]
            for i in range(len(seq)):
                srcall = np.loadtxt(path_txt + f'{seq[i]}.txt')
                bkgall = np.loadtxt(path_txt + f'{seq[i]}_bkg.txt')
                if src_pix[i]==0 or len(bkgall)==0 or len(srcall)==0:continue
                if bkgall.ndim==1: bkgall=np.array([bkgall])
                if srcall.ndim==1: srcall=np.array([srcall])
                else:
                    bkg_scale[i] = src_pix[i]/bkg_pix[i]
                    src_cts_temp=len(np.where(srcall[:,2]==obsid)[0])
                    bkg_cts_temp=len(np.where(bkgall[:,2]==obsid)[0])
                    bkg_cts[i]+=bkg_cts_temp*bkg_scale[i]

        src_info = np.column_stack((seq, ra, dec, counts_all,bkg_cts, 1/bkg_scale,exptime_all, VI_all, dist_all))
        np.savetxt(path + 'src_info.txt', src_info, fmt='%10d %10.5f %10.5f %10d %10.2f %10.2f %15d %10.5f %10.5f')
    else:
        src_info=np.column_stack((seq,ra,dec,counts_all,exptime_all,VI_all,dist_all))
        np.savetxt(path + 'src_info.txt', src_info, fmt='%10d %10.5f %10.5f %10d %15d %10.5f %10.5f')

if __name__ == '__main__':
    path='/Users/baotong/Desktop/period_Tuc/'
    extract_info(path, catname=catname[0], keyword=['RAJ2000','DEJ2000',gcname[0]], epoch_name='47Tuc_epoch', bkg=1)

