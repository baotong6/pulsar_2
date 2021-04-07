#!/bin/bash
# -*- coding: utf-8 -*-
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
from scipy.special import comb, perm
import pandas as pd
from astropy.timeseries import LombScargle
import scipy
import ephem

def DSOAZ(ra,dec,mydate):
    zuodijiang = ephem.Observer()
    zuodijiang.lat = '32.06'  #纬度
    zuodijiang.long = '118.77' #经度
    zuodijiang.elevation = 0  # 海拔
    zuodijiang.pressure = 0  # 压强
    # zuodijiang.date = now()
    zuodijiang.date = mydate  # 时间根据需要修改
    p = ephem.FixedBody(mydate)
    p._ra = ra
    p._dec = dec
    p.compute(zuodijiang)
    sun=ephem.Sun()
    sun.compute(zuodijiang)
    return (p.az,p.alt,sun.alt)

def make_obs_objlist(path,date,maglimits):
    #path中存放输出列表，CV星表文件'RK14.fit'
    #date为观测日期
    #maglimits为观测极限星等
    #该程序时间采样精度为1h

    BJtims = ["20:00:00","21:00:00","22:00:00","23:00:00","24:00:00","1:00:00","2:00:00","3:00:00"]
    EDTtims= ["12:00:00","13:00:00","14:00:00","15:00:00","16:00:00","17:00:00","18:00:00","19:00:00"] ##BJtime-8
    RK14=fits.open(path+'RK14.fit')
    RAJ2000=RK14[1].data['RAJ2000'];DECJ2000=RK14[1].data['DEJ2000'];
    NAME=RK14[1].data['Name'];mag1=RK14[1].data['mag1']
    NAME_list=[]
    os.chdir(path)
    os.system('rm CV_OBS_{0}.txt'.format(date.replace('/','')))
    f=open('CV_OBS_{0}.txt'.format(date.replace('/','')),'a+')
    headerline='NAME'.ljust(20,' ')+'RAJ2000'.ljust(20,' ')+'DECJ2000'.ljust(20,' ')+\
               'TSTART'.ljust(15,' ')+'TSTOP'.ljust(15,' ')+'Alt_max'.ljust(15,' ')+'mag1'.ljust(10,' ')+'\n'
    f.writelines(headerline)
    f.writelines('============================================================================================================='+'\n')
    for j in range(len(RAJ2000)):
        if mag1[j]>maglimits or mag1[j]==0:
            continue
        obstime=[];ALT_MAX=0
        for i in range(len(EDTtims)):
            ra = RAJ2000[j]
            dec = DECJ2000[j]
            (az,alt,sun_alt) = DSOAZ(ra, dec, mydate='2021/05/01'+' '+EDTtims[i])
            az = ephem.degrees(az)/2/np.pi*360
            alt = ephem.degrees(alt) / 2 / np.pi * 360
            sun_alt=ephem.degrees(sun_alt) / 2 / np.pi * 360
            if (alt>30) & (sun_alt<-18):
                obstime.append(BJtims[i])
                if alt>ALT_MAX:ALT_MAX=alt
        if len(obstime)>0:
            f.writelines(NAME[j].ljust(20,' ')+
                         str(RAJ2000[j].ljust(20,' '))+
                         str(DECJ2000[j].ljust(20, ' '))+
                         str(obstime[0]).ljust(15,' ')+
                         str(obstime[-1]).ljust(15,' ')+
                         str(np.round(ALT_MAX,2)).ljust(15,' ')+
                         str(mag1[j])+'\n')
    f.close()
path='/Users/baotong/Desktop/doc_file/CV_plan/'
make_obs_objlist(path,'2021/01/01',maglimits=14.0)
