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

path='/Users/baotong/Desktop/period/'
evt_list='SgrA_G_evt.fits'


# path='/Users/baotong/Desktop/period_gc/'
# evt_list='GC_merged_evt_FID.fits'

# path = '/Users/baotong/Desktop/period_LW/'
# evt_list='LW_merged_evt_FID.fits'

#catalog=pd.read_csv(path+'zhu3.csv')
#regname='pwn.reg'

def get_txt(regname):
    reg_file=[]
    hdul_evt= fits.open(path+evt_list)
    # x=hdul_evt[1].data.field(10)
    # y=hdul_evt[1].data.field(11)
    # energy=hdul_evt[1].data.field(14)
    # time=hdul_evt[1].data.field(0)
    # obs_ID=hdul_evt[1].data.field(19)

    x = hdul_evt[1].data.field(0)
    y = hdul_evt[1].data.field(1)
    energy = hdul_evt[1].data.field(2)
    time = hdul_evt[1].data.field(3)
    obs_ID = hdul_evt[1].data.field(11)

    def read_region(regname):
        with open(path+'region_G/'+regname, 'r') as file_to_read:
            while True:
                lines = file_to_read.readline() # 整行读取数据
                reg_file.append(lines)
                if not lines:
                    break
                    pass
        region=reg_file[-2][7:-2]
        reg_x,reg_y,reg_r=[float(i) for i in region.split(',')]
        return [reg_x,reg_y,reg_r]
    reg=read_region(regname)


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
    src_ID=obs_ID[src_index]

    def delete_photon_ID(time,energy,ID):
        i=0
        while i < len(energy):
            if energy[i]>8000 or energy[i]<2000:
                energy=np.delete(energy,i)
                time=np.delete(time,i)
                ID=np.delete(ID,i)
                i=i-1
            i=i+1
        return [time,energy,ID]

    [src_t,src_E,src_ID]=delete_photon_ID(src_t,src_E,src_ID)

    src_txt=np.column_stack((src_t,src_E,src_ID))
    src_txt = src_txt[src_txt[:,0].argsort()]
    #print src_txt
    np.savetxt(path+'txt_G_2_8k/'+regname[0:-4]+'.txt',src_txt,fmt="%.7f  %5.3f  %d")

cand_ID=[i+1 for i in range(518)]
cand_ID=np.loadtxt('G_src.txt')
cand_ID=cand_ID.astype(int)
# cand_ID=catalog['seq']
# for i in range(2401):
# # item='Cannonball'
#     get_txt(str(i)+'.reg')

# np.savetxt('cand_id.txt',cand_ID,format('%5d'))
# a=np.loadtxt('cand_id.txt')
# a=a.astype(int)
# print(a)
# print cand_ID
for item in cand_ID:
    get_txt(str(item)+'.reg')