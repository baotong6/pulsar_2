#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import string
from scipy.interpolate import lagrange
from scipy import interpolate

path='/Users/baotong/desktop/pulsar/pulsardata/B0531+21+fi/'
name=[];name_o=[];
with open(path+'evt/'+'data.txt','r+') as f:
    for line in f.readlines():
        name.append(line.strip())

    name = np.array(name)  # 将数据从list类型转换为array类型。
with open(path+'obt/'+'data.txt','r+') as f:
    for line in f.readlines():
        name_o.append(line.strip())

    name_o = np.array(name_o)  # 将数据从list类型转换为array类型。

def find_orbit(t):
    TIME_inter=np.append(TIME,TIME[-1]+1)
    ORB_inter=np.row_stack((ORB,ORB[-1]))
    x = TIME_inter
    #print ORB
    #print ORB_inter
    # if t>TIME[-1]:
    #     return ORB[-1]
    # if t<TIME[0]:
    #     return ORB[0]
    y0=ORB_inter[:,0]
    y1=ORB_inter[:,1]
    y2=ORB_inter[:,2]
    # x=x[0:10]
    #y0=y0[0:10]
    #y1=y1[0:10]
    #y2=y2[0:10]
    kind='cubic'
    f0=interpolate.interp1d(x, y0, kind)
    f1=interpolate.interp1d(x, y1, kind)
    f2=interpolate.interp1d(x,y2, kind)

    return np.array([f0(t),f1(t),f2(t)])
def get_fi():
    time = [];
    energy = [];
    TIME = [];
    ORB = [];
    VE = [];

    def find_orbit(t):
        TIME_inter = np.append(TIME, TIME[-1] + 1)
        TIME_inter=np.append(TIME[0]-1,TIME_inter)
        ORB_inter = np.row_stack((ORB, ORB[-1]))
        ORB_inter=np.row_stack((ORB[0],ORB_inter))
        x = TIME_inter
        # print ORB
        # print ORB_inter
        # if t>TIME[-1]:
        #     return ORB[-1]
        # if t<TIME[0]:
        #     return ORB[0]
        y0 = ORB_inter[:, 0]
        y1 = ORB_inter[:, 1]
        y2 = ORB_inter[:, 2]
        kind = 'cubic'
        f0 = interpolate.interp1d(x, y0, kind)
        f1 = interpolate.interp1d(x, y1, kind)
        f2 = interpolate.interp1d(x, y2, kind)

        return np.array([f0(t), f1(t), f2(t)])
    for i in range(len(name)):
        evtdata = np.genfromtxt(path + 'evt/' + name[i])
        time = evtdata[:, 0]
        energy = evtdata[:, 1]
        orbitdata = np.genfromtxt(path + 'obt/' + name_o[i])
        TIME = orbitdata[:, 0]
        ORB = np.transpose([orbitdata[:, 1], orbitdata[:, 2], orbitdata[:, 3]])
        VE = np.transpose([orbitdata[:, 4], orbitdata[:, 5], orbitdata[:, 6]])

        #print time
        #print TIME

        #print find_orbit(np.array([321544171,321544172]))
        #sate_r=np.sqrt((orb_event.dot(np.transpose(orb_event))).diagonal())
        orb_event=np.transpose(find_orbit(time))
        fi=[]
        for index in orb_event:
            fi.append(np.arcsin(index[2]/(index.dot(index)**0.5)))
        fi=np.array(fi)
        fi=fi/np.pi*180
        #print len(fi)

        evtdata_new=np.column_stack((evtdata,fi))
        evtdata_new=np.mat(evtdata_new)
        print(evtdata_new.shape[1])
        np.savetxt(path+'evt/new/'+name[i],evtdata_new,fmt="%.7f  %10.3f  %6.4f")

get_fi()


