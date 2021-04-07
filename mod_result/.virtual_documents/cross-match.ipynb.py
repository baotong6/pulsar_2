import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import pylab as pl
import string
from datetime import datetime
from scipy.interpolate import lagrange
from scipy import interpolate
from scipy.fftpack import fft,ifft
import scipy.signal as ss
import random
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle
import astropy.units as u
import astropy.coordinates.matching as match_catalog


def getlen(a,b):
    length=((a[0]-b[0])**2+(a[1]-b[1])**2)**0.5
    return length


def unit_convert(ra,dec):
    test = SkyCoord(ra,dec, unit = (u.hourangle,u.deg), frame='icrs')
    return [test.ra.degree,test.dec.degree]


# unit_convert(['17:45:43.036 -28:59:49.60','17:45:43.036 -28:59:49.60'])


def get_input():
    path='/Users/baotong/Desktop/period/'
    GCCR_file='GCCR_tab.txt'
    Xray_file='zhu18_3.fits'
    info_GCCR=[]
    ## GCCR的信息
    with open(path+GCCR_file,'r') as file_to_read:
        while True:
            lines = file_to_read.readline() # 整行读取数据
            info_GCCR.append(lines)
            if not lines:
                break
                pass
            
    info_GCCR=info_GCCR[0:-1]  ##去掉末尾空行
    label1=[];ra1=[];dec1=[];
    for i in range(len(info_GCCR)):  
        label_i,ra_i,dec_i=[str(i) for i in info_GCCR[i][0:-1].split(';')]   ##去掉末尾换行符
        label1.append(label_i);ra1.append(ra_i);dec1.append(dec_i)
    [ra1,dec1]=unit_convert(ra1,dec1)
#     ra1=np.array(ra1).astype(float)
#     dec1=np.array(dec1).astype(float)
    
    ## Zhu src的信息
    srclist=fits.open(path+Xray_file)
    label2=srclist[1].data['Seq']
    ra2 = srclist[1].data['RAJ2000']
    dec2= srclist[1].data['DEJ2000']
    
    ra_dec1=[ra1,dec1]
    ra_dec2=[ra2,dec2]
    return [ra_dec1,ra_dec2,label1,label2]
    
    


# get_input()


def compare_counterpart(ra_dec1,ra_dec2,seq1,seq2):

    offset=1
    ##arcsec

    match=[[] for i in range(len(seq1))]
    for i in range(len(seq1)):
        for j in range(len(seq2)):
            dis=getlen([ra_dec2[0][j],ra_dec2[1][j]],[ra_dec1[0][i],ra_dec1[1][i]])*3600
            if dis<offset:
                match[i].append([seq1[i],seq2[j],dis])
    print(np.sort(match))


input=get_input()
Radio_catalog = SkyCoord(ra=input[0][0]*u.degree, dec=input[0][1]*u.degree)
Xray_catalog = SkyCoord(ra=input[1][0]*u.degree, dec=input[1][1]*u.degree)
print(Xray_catalog)


id=match_catalog.match_coordinates_sky(Radio_catalog, Xray_catalog, nthneighbor=1)[0]+1
offset=match_catalog.match_coordinates_sky(Radio_catalog, Xray_catalog, nthneighbor=1)[1]
offset_arcsec=[]
for i in range(len(offset)):
    offset_arcsec.append(Angle(offset[i]).degree*3600)
# temp=SkyCoord(offset[0], unit = (u.deg), frame='icrs').to_string('decimal')


xray_counterpart=np.column_stack((id.astype(int),offset_arcsec))
np.savetxt('/Users/baotong/Desktop/period/GCCR_xrayc.txt',xray_counterpart,fmt='get_ipython().run_line_magic("d", " %10.4f')")



