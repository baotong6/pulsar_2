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
import pandas as pd
#path='/Volumes/halo/chandra_obs_hr_m04_for_baotong/'
path='/Volumes/halo/redgiant/'
#catalog=pd.read_csv(path+'zhu3.csv')
#regname='2338.reg'


def get_txt(evt_list):
    reg_file=[]
    #hdul_evt= fits.open(path+evt_list+'/'+'timing/'+'merged_evt.fits')
    hdul_evt = fits.open(path + evt_list + '/' + 'timing/' + '15761_reproj_evt.fits')
    x=hdul_evt[1].data.field(10)
    y=hdul_evt[1].data.field(11)
    energy=hdul_evt[1].data.field(14)
    time=hdul_evt[1].data.field(0)
    #obs_ID=hdul_evt[0].data.field(11)


    def read_region(regname):
        #with open(path+regname+'_s'+'.reg', 'r') as file_to_read:
        with open(path + regname + '.reg', 'r') as file_to_read:
            while True:
                lines = file_to_read.readline() # 整行读取数据
                reg_file.append(lines)
                if not lines:
                    break
                    pass
        region=reg_file[-2][7:-2]
        #region = reg_file[-2][11:-2]
        reg_x,reg_y,reg_r=[float(i) for i in region.split(',')]
        return [reg_x,reg_y,reg_r]
    reg=read_region(evt_list)

    def where_region(x,y,reg):
        r=np.array((x-reg[0],y-reg[1]))
        len_r=np.sqrt(r[0]**2+r[1]**2)
        temp=len_r-reg[2]
        return np.where(temp<=0)

    src_index=where_region(x,y,reg)
    src_x=x[src_index]
    src_y=y[src_index]
    src_t=time[src_index]
    src_E=energy[src_index]
    #src_ID=obs_ID[src_index]

    def delete_photon_ID(time,energy):
        i=0
        while i < len(energy):
            if energy[i]>8000 or energy[i]<500:
                energy=np.delete(energy,i)
                time=np.delete(time,i)
                i=i-1
            i=i+1
        return [time,energy]

    [src_t,src_E]=delete_photon_ID(src_t,src_E)

    src_txt=np.column_stack((src_t,src_E))
    src_txt = src_txt[src_txt[:,0].argsort()]
    #print src_txt
    #np.savetxt(path + 'txt/' + evt_list + '.txt', src_txt, fmt="%.7f  %5.3f ")
    #np.savetxt(path+'txt/'+evt_list+'_s'+'.txt',src_txt,fmt="%.7f  %5.3f ")
    np.savetxt(path+'txt/'+'evt_list15761'+'.txt',src_txt,fmt="%.7f  %5.3f ")

evt_list_ID=np.arange(1,49)
evt_list_ID=np.delete(evt_list_ID,[20,31,36,37,38,39,40,41,45])
#evt_list_ID=[33,34,35,36,43,44,45,47,48]
#evt_list_ID=np.arange(20,49)
evt_list_name=[]
for i in range(len(evt_list_ID)):
    if evt_list_ID[i]<10:
        evt_list_name.append('src'+'0'+str(evt_list_ID[i]))
    else:
        evt_list_name.append('src'+str(evt_list_ID[i]))

get_txt('src14')
# for i in range(len(evt_list_name)):
#     get_txt(evt_list_name[i])

#get_txt('src29')