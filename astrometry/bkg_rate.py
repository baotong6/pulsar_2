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
from scipy.optimize import curve_fit
import pandas as pd
import itertools

path='/Users/baotong/Desktop/pulsar/pulsardata/B0531+21/'
os.chdir(path)
file_obt=[]
with open ('./obt/data.txt','r') as f_obt:
    while True:
        lines = f_obt.readline()  # 整行读取数据
        file_obt.append(lines)
        if not lines:
            break
            pass
file_evt=['B0531+21'+file_obt[i][5:-1] for i in range(len(file_obt))]
file_obt=[file_obt[i][0:-1] for i in range(len(file_obt))]

bin=100
grid_back=np.zeros((bin,2*bin))
grid_cts=np.zeros((bin,2*bin))
grid_time=np.zeros((bin,2*bin))

def get_fi(x,y):
    if y>=0 and x>=0:
        fi=(np.arctan(y/x))/np.pi*180
    elif y>=0 and x<0:
        fi=180+np.arctan(y/x)/np.pi*180
    elif y<=0 and x>0:
        fi = 360+(np.arctan(y / x)) / np.pi * 180
    elif x<0 and y<0:
        fi = 180+(np.arctan(y /x)) / np.pi * 180
    return fi

def get_bkg_rate_inday(evtname,obtname):
    bkg_rate=[]
    os.chdir(path)
    if os.path.exists('./evt/{0}'.format(evtname)) and os.path.exists('./obt/{0}'.format(obtname)):
        time=np.loadtxt('./evt/{0}'.format(evtname))[:,0]
        energy=np.loadtxt('./evt/{0}'.format(evtname))[:,1]
        time_obt=np.loadtxt('./obt/{0}'.format(obtname))[:,0]
        x=np.loadtxt('./obt/{0}'.format(obtname))[:,1]
        y=np.loadtxt('./obt/{0}'.format(obtname))[:,2]
        z=np.loadtxt('./obt/{0}'.format(obtname))[:,3]
        r=(x**2+y**2+z**2)**0.5
        theta=np.arccos(z/r)/np.pi*180
        fi=np.zeros(len(y))


        for i in range(len(y)):
            fi[i]=get_fi(x[i],y[i])

        for i in range(len(time_obt)-1):
            index=np.where((time>time_obt[i])&(time<time_obt[i+1]))[0]
            cts=len(index)
            theta_temp = (theta[i] + theta[i + 1]) / 2
            fi_temp = (fi[i] + fi[i + 1]) / 2
            if cts!=0 :
                energy_bkg=energy[index]
                bkg_cts=len(np.where(energy_bkg>8000)[0])
                grid_back[int(theta_temp*100 / 180), int(fi_temp*200 / 360)] +=bkg_cts
                grid_cts[int(theta_temp*100 / 180), int(fi_temp*200 / 360)] +=cts
                grid_time[int(theta_temp * 100 / 180), int(fi_temp * 200 / 360)] += time_obt[i + 1] - time_obt[i]
            elif time_obt[i+1]-time_obt[i]<100:
                grid_time[int(theta_temp*100 / 180), int(fi_temp*200 / 360)]+=time_obt[i+1]-time_obt[i]


        return [time,energy,time_obt,r/1000.,theta,fi]
# time_all=0;cts_all=0
# for i in range(len(file_obt)):
#     get_bkg_rate_inday(file_evt[i],file_obt[i])
#     print(i)

# np.save('grid_cts.npy',grid_cts)
# np.savetxt('grid_cts.txt',grid_cts)
# np.save('grid_back.npy',grid_back)
# np.savetxt('grid_back.txt',grid_back)
# np.save('grid_time.npy',grid_time)
# np.savetxt('grid_time.txt',grid_time)
#
grid_cts=np.load('grid_cts.npy')
grid_back=np.load('grid_back.npy')
grid_time=np.load('grid_time.npy')

grid=grid_back/grid_time

grid_cts_rate=list(itertools.chain.from_iterable(grid))
i=0
while i <len(grid_cts_rate):
    if np.isnan(grid_cts_rate[i]):
        grid_cts_rate=np.delete(grid_cts_rate,i)
    else:
        i+=1

# np.savetxt('grid.txt',grid)
for i in range(len(grid)):
    for j in range(len(grid[i])):
        if grid[i][j]<0.01:
            grid[i][j]=0.01

plt.contourf(np.linspace(0,360,200),np.linspace(90,-90,100),np.log10(grid))
plt.colorbar()
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.savefig('grid.eps')

# plt.hist(grid_cts_rate,bins=100000,cumulative=True,normed=True,histtype = 'step')
# plt.xlim(0,10)
# plt.ylim(0,1)
# plt.xlabel('Counts rate (>8 keV)')
# plt.ylabel('Probability density')
# plt.savefig('cdf.eps')
plt.show()
