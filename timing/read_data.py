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
import scipy.signal as signal

from collections import OrderedDict
path='/Users/baotong/Desktop/period/'
#dataname="zhu170.txt"
#dataname='Magnetar_S_lc.txt'
def get_data(dataname):
    time=[];energy=[];obs_ID=[]
    #with open(path+'txt_90/'+dataname, 'r') as file_to_read:
    with open(path + 'txt_all_obs_I/' + dataname, 'r') as file_to_read:
    #with open(path + 'txt_merge_2_8k/' + dataname, 'r') as file_to_read:
    #with open(path + 'txt_2_8k/' + dataname, 'r') as file_to_read:
        while True:
            lines = file_to_read.readline() # 整行读取数据
            if not lines:
                break
                pass
            p_tmp, E_tmp, obsID_temp = [float(i) for i in lines.split()] # 将整行数据分割处理，如果分割符是空格，括号里就不用传入参数，如果是逗号， 则传入‘，'字符。
            #print obsID_temp
            #print dataname
            time.append(p_tmp)  # 添加新读取的数据
            energy.append(E_tmp)
            obs_ID.append(obsID_temp)
        pass
        time = np.array(time)
        energy = np.array(energy)
        obs_ID=np.array(obs_ID)
    # print time
    # print energy
    #print obs_ID
    dict_obsinfo = OrderedDict()
    #epoch=np.loadtxt(path+'txt_90/'+'SgrA_I_epoch.txt')
    epoch = np.loadtxt(path + 'txt_all_obs_I/' + 'epoch_src_'+dataname)
    #epoch = np.loadtxt(path + 'txt_merge/' + 'SgrA_IG_epoch.txt')
    #epoch = np.loadtxt(path + 'txt_merge_2_8k/' + 'SgrA_IG_epoch.txt')
    #epoch = np.loadtxt(path + 'txt_2_8k/' + 'SgrA_I_epoch.txt')
    for i in range(len(epoch)):
        dict_obsinfo[int(epoch[i][2])]=[epoch[i][3],np.array([epoch[i][0],epoch[i][1]])]
    #dict_obsinfo[epoch[2]]=
    # exptime=epoch[:,1]-epoch[:,0]
    #np.savetxt(path+'txt/'+'exptime_G.txt',exptime,fmt='%.2f')
    obstime=epoch[:,0:2]

    ID=[]
    cut=[]
    for i in range(len(time)-1):
        if obs_ID[i+1]!=obs_ID[i]:
            cut.append(i+1)
    t=[0 for j in range(len(cut)+1)]
    E=[0 for j in range(len(cut)+1)]
    t[0] = time[0:cut[0]]
    E[0] = energy[0:cut[0]]

    dict=OrderedDict()
    dict[int(obs_ID[0])]=[t[0],E[0],dict_obsinfo[int(obs_ID[0])][0],dict_obsinfo[int(obs_ID[0])][1]]

    ID.append(int(obs_ID[0]))
    for j in range(1,len(cut)):
        t[j]=time[cut[j-1]:cut[j]]
        E[j]=energy[cut[j-1]:cut[j]]

        dict[int(obs_ID[cut[j-1]])]=[t[j],E[j],dict_obsinfo[int(obs_ID[cut[j-1]])][0],dict_obsinfo[int(obs_ID[cut[j-1]])][1]]
        #print int(obs_ID[cut[j-1]]),dict[int(obs_ID[cut[j-1]])][-2]
        ID.append(int(obs_ID[cut[j-1]]))
    t[-1]=time[cut[-1]:]
    E[-1]=energy[cut[-1]:]
    dict[int(obs_ID[-1])]=[t[-1],E[-1],dict_obsinfo[int(obs_ID[-1])][0],dict_obsinfo[int(obs_ID[-1])][1]]
    ID.append(int(obs_ID[-1]))

    #print ID
    #print 'run1'


    sss=0
    for index in t:
        for i in range(len(index)):
            sss+=1
    if sss!=len(time):
        print('error')
    # print len(time)
    return [time,energy,obs_ID,t,E,dict]

def tran_t_time(t):
    time=np.array(t[0])
    for i in range(1,len(t)):
        time=np.concatenate((time,t[i]))

    return time
#dict=get_data('Muno4.txt')[5]
#print dict[10556]
#print get_data('Muno1.txt')[2]
#print dict[242][-1]
#print exptime
# print sss
# print len(time)
# print sum(exptime)
#plt.scatter(plot_obstime,plot_obstime-plot_obstime+1)
#plt.show()
# print exptime
# print 1/exptime
#plt.scatter(time,[1 for i in range(len(time))])
# t_test=time
#a=plt.hist(t_test[0],int((t_test[-1]-t_test[0])/12.6))

