# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import string

path='/Users/baotong/Desktop/period/simulation/simulation_NSC_eclipse_554/'
#cts_range=[1,2,3,4,5]
cts_range=[10.0,11.0,12.0,13.0,14.0,15.0]
width_range=[0.2,0.3,0.4,0.5]
ec_amp=0.9
#sim_N=10 #单次模拟的源数目
sim_N=50
sim_id_range=np.linspace(1,sim_N,sim_N)
threshold=0.9
sim_id_range=sim_id_range.astype(int)

detect_rate=[]
mean=[]
var=[]
std=[]
cts=[]
period_real = 554.
for width in width_range:
    res_detect = []
    res_mean = []
    res_var = []
    res_std = []
    res_cts=[]
    for cts_rate in cts_range:
        detect=0
        period_get=[]
        cts_get=[]
        for i in sim_id_range:
            temp_info=np.loadtxt(path+'result_{0}_{1}_{2}/result_sim_{3}.txt'.format(str(cts_rate),str(ec_amp),str(width),str(i)))
            cts_get.append(temp_info[-1])
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
        res_cts.append(np.mean(np.array(cts_get)))
        res_detect.append(detect/sim_N)
        res_mean.append(period_mean)
        res_var.append(period_var)
        res_std.append(period_std)

    detect_rate.append(res_detect)
    mean.append(res_mean)
    var.append(res_var)
    std.append(res_std)
    cts.append(res_cts)
print(detect_rate)
print(cts)

def make_plot():
    # plt.figure(1)
    # plt.title('detection')
    # for i in range(len(cts)):
    #     plt.plot(cts[i],detect_rate[i],marker='v')
    #
    # #plt.legend(['cr=1','cr=2','cr=3','cr=4','cr=5'])
    # plt.xlabel('cts')
    # plt.ylabel('detection-rate')
    # plt.legend(['width=0.2', 'width=0.3', 'width=0.4', 'width=0.5'],loc='best')
    # #plt.savefig('detection.eps')
    # plt.savefig('./fig/'+'detection_{0}.eps'.format(str(period_real)))
    # plt.show()
    # plt.close(1)

    plt.figure(2)
    #plt.semilogy()
    plt.title('period_mean')
    plt.xlabel('cts')
    plt.ylabel('p-mean')
    #plt.ylim(period_real*(1-0.04), period_real*(1+0.04))
    plt.plot([0,10000], np.zeros(2) + period_real, '--')
    for i in range(len(cts)):
        plt.plot(cts[i], mean[i], marker = 'v')
    plt.ylim(period_real*0.99,period_real*1.01)
    #plt.legend(['true_p', 'cr=1', 'cr=2', 'cr=3', 'cr=4', 'cr=5'])
    plt.legend(['true_p','width=0.2', 'width=0.3', 'width=0.4', 'width=0.5'], loc = 'best')
    plt.savefig('./fig/'+'period_mean_{0}.eps'.format(str(period_real)))
    plt.show()
    plt.close(2)

    # plt.figure(3)
    # plt.title('period_var')
    # plt.xlabel('amplitude')
    # plt.ylabel('p-var')
    # plt.semilogy()
    # for i in range(len(cts_range)):
    #     plt.plot(amp_range,var[i],marker='v')
    #
    # #plt.legend(['cr=1','cr=2','cr=3','cr=4','cr=5'])
    # plt.legend(['cr=2.1', 'cr=2.2', 'cr=2.3', 'cr=2.4', 'cr=2.5','cr=2.6','cr=2.7','cr=2.8','cr=2.9'],loc='best')
    # plt.savefig('./fig/'+'period_var_{0}.eps'.format(str(period_real)))
    # plt.close(3)
    #
    # plt.figure(4)
    # plt.title('period_std')
    # plt.xlabel('amplitude')
    # plt.ylabel('p-std')
    # plt.semilogy()
    # for i in range(len(cts_range)):
    #     plt.plot(amp_range,std[i],marker='v')
    #
    # #plt.legend(['cr=1','cr=2','cr=3','cr=4','cr=5'])
    # plt.legend(['cr=2.1', 'cr=2.2', 'cr=2.3', 'cr=2.4', 'cr=2.5','cr=2.6','cr=2.7','cr=2.8','cr=2.9'],loc='best')
    # plt.savefig('./fig/'+'period_std_{0}.eps'.format(str(period_real)))
    # plt.close(4)
make_plot()