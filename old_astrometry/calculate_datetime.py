#!/bin/bash
# -*- coding: utf-8 -*-
#pulsation=0.0331
import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
import datetime as datetime




path='/Users/baotong/desktop/pulsar/pulsardata/B0531+21_chandra/'
name=[];name_o=[];
with open(path+'evt/'+'data.txt','r+') as f:
    for line in f.readlines():
        name.append(line.strip())

    name = np.array(name)  # 将数据从list类型转换为array类型。

for i in range(len(name)):
    evtdata = np.genfromtxt(path + 'evt/' + name[i])
    time = evtdata[:, 0]
    print(time)
    energy = evtdata[:, 1]
    t_jd = 2454466.5 + time / 86400.0

    def clock_tt(t):

        # t must be format "jd"
        time = Time(t, format='jd', scale='utc')
        time2=Time(2400000.5+50814.0,format='jd',scale='tt')
        time_tt = time.tt
        return (time_tt.value-time2.value)*86400

    t_new=clock_tt(t_jd)

    evtdata_new=np.column_stack((t_new,energy))
    np.savetxt(path + 'evt/new/' + name[i], evtdata_new,fmt="%.7f  %.2f")
