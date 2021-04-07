#!/bin/bash
# -*- coding: utf-8 -*-
#pulsation=0.0331
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import string
import math
from jplephem.spk import SPK
from astropy.time import Time
from astropy.coordinates import FK5
import sofa_test as sofa_test

#constant
#np.set_printoptions(suppress=True)
def bary(time,orbit):
    c=299792.458
    #km
    G=6.67259*10**(-11)
    Msun=1.9891*10**(30)
    ra=(5+34./60+31.97232/3600)*15
    dec=(22+0./60.+52.0690/3600)*15
    distance=6300*9460730472580.8
    #(取6500光年,最后单位为km)

    #jd = np.array([2457061.5, 2457062.5, 2457063.5, 2457064.5])

    def jpl(number,jd):
        kernel = SPK.open('/Users/baotong/desktop/pulsar/de430.bsp')
        position, velocity = kernel[0, number].compute_and_differentiate(jd)

        #position=np.array(position)
        #velocity=np.array(velocity)
        velocity=velocity/86400.0
        return[position,velocity]

    #print jpl(3,jd)[0]

    def angle(x,y):
        Lx=x.dot(x)
        Ly=y.dot(y)
        x=x/(Lx**0.5)
        #print x
        y=y/(Ly**0.5)
        #print y
        cos_angle=x.dot(y)
        angle=np.arccos(cos_angle)
        return angle

    # with open(path+name[-1], 'r') as file_to_read:
    #     while True:
    #         lines = file_to_read.readline() # 整行读取数据
    #         if not lines:
    #             break
    #             pass
    #         p_tmp, E_tmp = [float(i) for i in lines.split()] # 将整行数据分割处理，如果分割符是空格，括号里就不用传入参数，如果是逗号， 则传入‘，'字符。
    #         time.append(p_tmp)  # 添加新读取的数据
    #         energy.append(E_tmp)
    #     pass
    #     time = np.array(time)
    #     energy = np.array(energy)

    #tb = tobs + (clock) - (dispersion) + (geometric) + (Einstein) - (Shapiro)
    # (dispersion) - are dispersion corrections of the form D/fb^2
    # where fb is the barycentric frequency.  D=0 for X-rays.

    def TDB_UTC(t):
    # t must be format "jd"
        time=Time(t,format='jd',scale='utc')
        time_utc=time
        time_tt=time.tt
        time_tdb=time.tdb
        tdb_tt=time_tt.delta_tdb_tt
        tt_utc=86400*(time_tt.value-time.value)
        #print tdb_tt
        #print tt_utc
        dt=tdb_tt+tt_utc

        return dt
    tobs=time
    clock=TDB_UTC(time)
    dispersion=0
    #robs is the vector from the SSB to the observatory
    #rearth is the vector from the SSB to the earth geocenter
    #rsat is the vector from the earth geocenter to the spacecraft
    rsat = sofa_test.r_CIS(orbit,tobs)
    #print rsat
    rearth=jpl(3,tobs)[0]
    vearth=jpl(3,tobs)[1]
    rsun=jpl(10,tobs)[0]
    n=np.array([np.cos(ra)*np.cos(dec),np.sin(ra)*np.cos(dec),np.sin(dec)])
    rpulsar=distance*n
    #print rearth
    #print rsun
    #print rpulsar

    th=angle(rearth-rsun,rpulsar-rsun)
    #print np.cos(th)
    robs = rearth + rsat
    # n is a unit vector pointing to the astrophysical source
    # The value of n depends on the position of the source on the sky (ie,
    # RA and DEC), is assumed to be known already, or is solved for in using
    # an iterative process.
    geometric = (robs.dot(n))/c
    Einstein_geo=1.48*10**(-8)
    Einstein=Einstein_geo+(rsat.dot(vearth))/(c**2)
    #print Einstein
    Shapiro = - (2*G*Msun/((c*1000)**3))*np.log(1+np.cos(th))
    #print Shapiro
    #print geometric
    #print Einstein


    tb = (tobs-2454466.5)*86400 + clock - dispersion + geometric + Einstein - Shapiro
    #print tb-(tobs-2454466.5)*86400
    # print robs
    # print geometric
    # print clock
    # print Einstein
    # print Shapiro
    return tb

x=np.array([5.353913581681e+006,-1.043885119003e+006,4.196499926732e+006])
T_test=2.829595769999998e8
# y=index/86400.0+2454466.5
# print y
y=bary(T_test/86400.0+2454466.5,x)
# print y-index
#z=(2457061.5*86400)

