# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import string

path='/Users/baotong/Desktop/period_LW/simulation/simulation_LW_554/'
#cts_range=[1,2,3,4,5]
cts_range=[0.5,1.0,2.0,3.0,4.0,5.0]
# cts_range=[2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9]
amp_range=[0.5,0.6,0.7,0.8,0.9]
# amp_range=[0.5]
#sim_N=10 #单次模拟的源数目
sim_N=100
threshold=0.9
period_real = 554.
sim_id_range=np.linspace(1,sim_N,sim_N)
sim_id_range=sim_id_range.astype(int)

detect_rate=[]
mean=[]
var=[]
std=[]
for cts_rate in cts_range:
    res_detect = []
    res_mean = []
    res_var = []
    res_std = []
    for amp in amp_range:
        detect=0
        period_get=[]
        for i in sim_id_range:
            temp_info=np.loadtxt(path+'result_{0}_{1}/result_sim_{2}.txt'.format(str(cts_rate),str(amp),str(i)))

            if temp_info[2]>threshold and 0.2*period_real<temp_info[4]<1.2*period_real:
                detect+=1
                period_get.append(temp_info[4])
        period_get_array=np.array(period_get)

        if len(period_get_array)==0:
            period_mean=0
            period_var=0
            period_std=0
        else:
            period_mean=np.mean(period_get_array)  ##均值
            period_var=np.var(period_get_array)    ##方差
            period_std=np.std(period_get_array,ddof=1)  ##标准差

        res_detect.append(detect/sim_N)
        res_mean.append(period_mean)
        res_var.append(period_var)
        res_std.append(period_std)

    detect_rate.append(res_detect)
    mean.append(res_mean)
    var.append(res_var)
    std.append(res_std)
print(detect_rate)
print(mean)

def make_plot():
    path_out='/Users/baotong/Desktop/period_LW/simulation/fig_554/'
    os.chdir(path_out)
    plt.figure(1)
    # plt.title('detection')
    for i in range(len(cts_range)):
        plt.plot(amp_range,detect_rate[i],marker='v')

    #plt.legend(['cr=1','cr=2','cr=3','cr=4','cr=5'])
    plt.xlabel('Amplitude')
    plt.ylabel('Detection-rate')
    plt.legend(title='P=554s',labels=['cts=50','cts=100', 'cts=200', 'cts=300', 'cts=400', 'cts=500'],fontsize='small',loc=[0.81,0.08])
    # plt.text(0.7,0.4,'P=5540s')
    #plt.savefig('detection.eps')

    plt.savefig('Detection_{0}.eps'.format(str(int(period_real))))
    plt.close(1)

    plt.figure(2)
    #plt.semilogy()
    plt.title('period_mean')
    plt.xlabel('amplitude')
    plt.ylabel('p-mean')
    plt.ylim(period_real*(1-0.001), period_real*(1+0.001))
    plt.plot(amp_range, np.zeros(len(amp_range)) + period_real, '--')
    for i in range(len(cts_range)):
        plt.plot(amp_range, mean[i], marker = 'v')
    #plt.ylim(5400,5700)
    #plt.legend(['true_p', 'cr=1', 'cr=2', 'cr=3', 'cr=4', 'cr=5'])
    plt.legend(['true_p','cts=50','cts=100', 'cts=200', 'cts=300', 'cts=400', 'cts=500'],loc='best')
    plt.savefig('period_mean_{0}.eps'.format(str(period_real)))
    plt.close(2)

    plt.figure(3)
    plt.title('period_var')
    plt.xlabel('amplitude')
    plt.ylabel('p-var')
    plt.semilogy()
    for i in range(len(cts_range)):
        plt.plot(amp_range,var[i],marker='v')

    #plt.legend(['cr=1','cr=2','cr=3','cr=4','cr=5'])
    plt.legend(['cts=50','cts=100', 'cts=200', 'cts=300', 'cts=400', 'cts=500'],loc='best')
    plt.savefig('period_var_{0}.eps'.format(str(period_real)))
    plt.close(3)

    plt.figure(4)
    plt.title('period_std')
    plt.xlabel('amplitude')
    plt.ylabel('p-std')
    plt.semilogy()
    for i in range(len(cts_range)):
        plt.plot(amp_range,std[i],marker='v')

    #plt.legend(['cr=1','cr=2','cr=3','cr=4','cr=5'])
    plt.legend(['cts=50','cts=100', 'cts=200', 'cts=300', 'cts=400', 'cts=500'],loc='best')
    plt.savefig('period_std_{0}.eps'.format(str(period_real)))
    plt.close(4)
make_plot()