# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import string
import correct_test as correct_test
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

evtdata=np.genfromtxt(path+'evt/'+name[-1])
time=evtdata[:,0]
energy=evtdata[:,1]
orbitdata = np.genfromtxt(path+'obt/'+name_o[-1])
TIME = orbitdata[:,0]
ORB=np.transpose([orbitdata[:,1],orbitdata[:,2],orbitdata[:,3]])
VE=np.transpose([orbitdata[:,4],orbitdata[:,5],orbitdata[:,6]])

# def readtxt(name):
#     time=[];energy=[]
#     for i in range(len(name)):
#         with open(path+'evt/'+name[i], 'r') as file_to_read:
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

#print ORB
cut=[]
for i in range(len(time)-1):
    if time[i+1]-time[i]>100:
        cut.append(i+1)
#print len(cut)
T=[0 for i in range(len(cut)+1)]
E=[0 for i in range(len(cut)+1)]
T[0]=time[0:cut[0]]
E[0]=energy[0:cut[0]]
for i in range(1,len(cut)):
    T[i]=time[cut[i-1]+1:cut[i]]
    E[i]=energy[cut[i-1]+1:cut[i]]
T[-1]=time[cut[-1]+1:]
E[-1]=energy[cut[-1]+1:]
def delete_photon(time,energy):
    i=0
    while i < len(energy):
        if energy[i]>5000 or energy[i]<500:
            energy=np.delete(energy,i)
            time=np.delete(time,i)
            i=i-1
        i=i+1

    return [time,energy]

v=29.6328827000
vdot=-369161.97*1e-15
v=v+vdot*222598


evt_time=T[-1]
evt_energy=E[-1]
[evt_time,evt_energy]=delete_photon(T[-1],E[-1])



turns=[]
for item in evt_time:
    turns.append(v*item-int(v*item))
#print turns
bin=500
loc=[0 for i in range(bin)]
turns=sorted(turns)

for item in turns:
    loc[int(item/(1.0/bin))]+=1

x=np.arange(0,1,1.0/bin)
plt.plot(x,loc)
plt.show()

