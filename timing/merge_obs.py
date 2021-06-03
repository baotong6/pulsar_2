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
from astropy.table import Table, Column

#path='/Volumes/halo/GC/merge_data/xdata/'
#obsID=['945','14897','17236','17239','17237','18852','17240','17238','20118','17241','20807','20808']

# path='/Volumes/pulsar/LimWin/merge_data/xdata/'
# obsID=['6362','5934','6365','9505','9855','9502','9500','9501','9854','9503','9892','9893','9504']
# label=['ACIS-I','ACIS-S','LW']

path='/Volumes/pulsar/ECDFS/merge_data/xdata/'
# obsID=['3798','10059','13225','13252','13705','13706' ,
#         '14339','14475','14476','14477','14478','14479',
#         '14625','15615','15750','16638','17779','18881']
# label=['terzan5']
obsID=['5019','5020']
label=['ECDFS']
def get_epoch_file(obsID,label):
    TSTART = []
    TSTOP = []

    for item in obsID:
        #hdul=fits.open(path+'GC_'+item+'_reproj_evt.fits')
        hdul = fits.open(path + 'all_bcc_'+item+'_reproj_evt.fits')
        TSTART.append(hdul[1].header['TSTART'])
        TSTOP.append(hdul[1].header['TSTOP'])

    TSTART = np.array(TSTART)
    TSTOP = np.array(TSTOP)
    obsID=np.array(obsID)
    obsID=obsID.astype('int')
    epoch_info=np.column_stack((TSTART,TSTOP,obsID,TSTOP-TSTART))
    epoch_info = epoch_info[epoch_info[:, 0].argsort()]
    print(epoch_info)
    np.savetxt(path+label[0]+'_epoch.txt',epoch_info,fmt='%15.2f %15.2f %10d %20.2f')
get_epoch_file(obsID,label)

def get_index_in_list(list,a):
    #返回list中最接近a的index
    #其中list必须为np.array格式
    list_u=list-a
    list_u=np.abs(list_u)
    NUM=np.where(list_u==np.min(list_u))[0]
    return NUM

def mod_evt(label):
    print('run')
    epoch_info=np.loadtxt(path+label[2]+'_epoch.txt')
    TSTART=epoch_info[:,0]
    TSTOP=epoch_info[:,1]
    obsID=epoch_info[:,2]
    exptime=epoch_info[:,3]
    file_o='LW_merged_evt.fits'
    hdul=fits.open(path+file_o)
    hd1=hdul[1]
    time=hdul[1].data.field(0)
    FID=[]
    for i in range(len(obsID)):
        if i >=1:
            if get_index_in_list(time,TSTART[i])[0]<=get_index_in_list(time,TSTOP[i-1])[-1]:
                ##防止出现几个光子时间相同的情况而重复
                FID.append(np.zeros(get_index_in_list(time,TSTOP[i])[-1]-get_index_in_list(time,TSTOP[i-1])[-1])+obsID[i])
            else:
                FID.append(np.zeros(get_index_in_list(time, TSTOP[i])[-1] - get_index_in_list(time, TSTART[i])[0] + 1) + obsID[i])
        else:
            FID.append(
                np.zeros(get_index_in_list(time, TSTOP[i])[-1] - get_index_in_list(time, TSTART[i])[0] + 1) + obsID[i])




    temp=0
    for i in range(len(FID)):
        temp+=len(FID[i])
    if temp!=len(time):
        print(temp)
        print(len(time))
        print('error')
        return 'error'

    FID_1d=np.concatenate((FID[0],FID[1]))
    for i in range(2,len(FID)):
        FID_1d=np.concatenate((FID_1d,FID[i]))
    FID_1d=FID_1d.astype('int')
    #得到前一个表
    t1 = Table.read(path+file_o)
    col1 = Table.Column(name='FID', data=FID_1d)
    t1.add_column(col1)
    #得到第二个表
    print(t1)
    t1.write(path+'LW_merged_evt_FID.fits',format='fits',overwrite='True')

#mod_evt(label)