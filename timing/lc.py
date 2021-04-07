#!/bin/bash
# -*- coding: utf-8 -*-
#pulsation=0.0331
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import string
# import correct as correct
from datetime import datetime
from scipy.interpolate import lagrange
from scipy import interpolate
path='/Users/baotong/desktop/pulsar/pulsardata/B0531+21/'

name=[];time=[];energy=[];name_o=[];orbit=[0,0,0];v=[0,0,0]
TIME=[];ORB=[];VE=[]
with open(path+'evt/'+'data.txt') as f:
    for line in f.readlines():
        name.append(line.strip())

    name = np.array(name)  # 将数据从list类型转换为array类型。
with open(path+'obt/'+'data.txt') as f:
    for line in f.readlines():
        name_o.append(line.strip())

    name_o = np.array(name_o)  # 将数据从list类型转换为array类型。
#print(name)

# def readtxt(name):
#     time=[];energy=[]
#     for i in range(len(name)):
#         with open(path+name[i], 'r') as file_to_read:
#             while True:
#                 lines = file_to_read.readline() # 整行读取数据
#                 if not lines:
#                     break
#                     pass
#                 p_tmp, E_tmp = [float(i) for i in lines.split()] # 将整行数据分割处理，如果分割符是空格，括号里就不用传入参数，如果是逗号， 则传入‘，'字符。
#                 time.append(p_tmp)  # 添加新读取的数据
#                 energy.append(E_tmp)
#             pass
#     time = np.array(time)
#     energy = np.array(energy)
#     return [time,energy]
# time=readtxt(name)[0];energy=readtxt(name)[1]

with open(path+'evt/'+name[-1], 'r') as file_to_read:
    while True:
        lines = file_to_read.readline() # 整行读取数据
        if not lines:
            break
            pass
        p_tmp, E_tmp = [float(i) for i in lines.split()] # 将整行数据分割处理，如果分割符是空格，括号里就不用传入参数，如果是逗号， 则传入‘，'字符。
        time.append(p_tmp)  # 添加新读取的数据
        energy.append(E_tmp)
    pass
    time = np.array(time)
    energy = np.array(energy)

with open(path+'obt/'+name_o[-1], 'r') as file_to_read:
    while True:
        lines = file_to_read.readline() # 整行读取数据
        #print lines
        if not lines:
            break
            pass
        p_tmp,orbit[0],orbit[1],orbit[2],v[0],v[1],v[2] = [float(i) for i in lines.split()] # 将整行数据分割处理，如果分割符是空格，括号里就不用传入参数，如果是逗号， 则传入‘，'字符。
        TIME.append(p_tmp)
        #print orbit
        ORB.append(orbit)  # 添加新读取的数据
        VE.append(v)
    pass
    data = np.genfromtxt(path+'obt/'+name_o[-1])

    TIME = data[:,0]
    ORB=np.transpose([data[:,1],data[:,2],data[:,3]])
    VE=np.transpose([data[:,4],data[:,5],data[:,6]])
    #print ORB[:,0]
    #print TIME

#print ORB
#print TIME

cut=[]
for i in range(len(time)-1):
    if time[i+1]-time[i]>100:
        cut.append(i+1)
t1=time[0:cut[0]]
t2=time[cut[0]+1:cut[1]]
#t3=time[cut[1]+1:]
energy_2=energy[cut[0]+1:cut[1]]
#t2=time
#energy_2=energy
#print len(t2),len(energy_2)
print(t2[1]-t2[0])

def delete_photon(time,energy):
    i=0
    while i < len(energy):
        if energy[i]>5000 or energy[i]<500:
            energy=np.delete(energy,i)
            time=np.delete(time,i)
            i=i-1
        i=i+1

    return [time,energy]
[t2,energy_2]=delete_photon(t2,energy_2)
#plt.scatter(t2,energy_2)
#plt.show()

#print len(delete_photon(t2,energy_2)[1]
# ZOOM=np.append(t2,TIME)
# print len(ZOOM)
# location=[]
# ZOOM=np.sort(ZOOM)
# for i in range(len(TIME)):
#     #print np.argwhere(ZOOM==index)
#     lll=int(np.argwhere(ZOOM==TIME[i]))
#     location.append([i,lll])
#
# location=np.array(location)
# print location[:,1]
# check=[]
# print TIME
#
# def find_orbit(x,orbit_time):
#      a=[]
#      for i in range(len(orbit_time)):
#          a.append(np.fabs(x-orbit_time[i]))
#      #print a
#      return a.index(min(a))
#  # print t2
#  # x=correct.bary(t2[0]/86400.0+2454466.5,ORB[find_orbit(t2[0],TIME)])
#  # print x
# #print ORB[find_orbit(t2,TIME)]
# #print ORB


def find_orbit(t):
    x=TIME
    x=x-282614400
    #print x
    y0=ORB[:,0]
    y1=ORB[:,1]
    y2=ORB[:,2]
    # x=x[0:10]
    #y0=y0[0:10]
    #y1=y1[0:10]
    #y2=y2[0:10]
    kind='cubic'
    f0=interpolate.interp1d(x, y0, kind)
    f1=interpolate.interp1d(x, y1, kind)
    f2=interpolate.interp1d(x,y2, kind)

    return np.array([float(f0(t)),float(f1(t)),float(f2(t))])
#print find_orbit(t2[0]-282614400)

def TB():
    TB=[]
    for index in t2:
        TBindex=correct.bary(index/86400.0+2454466.5,find_orbit(index-282614400))
        #print index-TBindex
        TB.append(TBindex)

        #print index


    return TB

TB=TB()
np.save("right_time.npy",TB)
np.save("obs_time.npy",t2)
print('done')
