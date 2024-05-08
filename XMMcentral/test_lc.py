'''
Author: baotong && baotong@smail.nju.edu.cn
Date: 2024-02-27 09:18:44
LastEditors: baotong && baotong@smail.nju.edu.cn
LastEditTime: 2024-02-27 10:17:23
FilePath: /pulsar/XMMcentral/test_lc.py
Description: 

Copyright (c) 2024 by baotong, All Rights Reserved. 
'''
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.stats import poisson_conf_interval
from scipy.stats import norm
from stingray.events import EventList
from astropy.timeseries import LombScargle
from stingray.lightcurve import Lightcurve
import os
font1 = {'family': 'Normal',
         'weight': 'normal',
         'size': 18, }

def readconfirmed():
    path = '/Users/baotong/Desktop/XMMcentral/all_lc/confirmed'
    file_names = os.listdir(path);list=[]
    for file_name in file_names:
        if file_name.endswith(".pdf"):
            list.append(file_name[:-7])
    return list
def readlc():
    path = '/Users/baotong/Desktop/XMMcentral/all_lc/'
    # file_names = os.listdir(path);list=[]
    # for file_name in file_names:
    #     if file_name.endswith(".lc"):
    #         list.append(file_name)
    list=readconfirmed()
    for i in range(len(list)):
        dataname=list[i]+'.lc'
        a = fits.open(path + dataname)
        rate = a[1].data['RATE']
        time = a[1].data['TIME']
        GTI=a[2].data
        GTI=np.array(GTI)
        epoch=[]
        for j in range(len(GTI)):
            tstart=GTI[j][0]
            tstop=GTI[j][1]
            tlast=tstop-tstart
            epoch.append([tstart,tstop,tlast])
        epoch=np.array(epoch)
        gaptime=np.sum(epoch[:,0][1:]-epoch[:,1][0:-1])
        print(gaptime)
        print(dataname)
        T = time[-1] - time[0]
        lc = Lightcurve(time=time, counts=rate * (time[1] - time[0]))
        # print(lc.dt)

if __name__=='__main__':
    readlc()

