# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import string

path_1='/Users/baotong/Desktop/period_LW/simulation/simulation_LW_eclipse_5258/'
path_2='/Users/baotong/Desktop/period_LW/simulation/simulation_LW_eclipse_15258/'
path_3='/Users/baotong/Desktop/period_LW/simulation/simulation_LW_eclipse_45258/'
cts_range=[1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,11.0,15.0,20.0]
# cts_range=[10.0,11.0,12.0,13.0,14.0,15.0]
width=0.9
amp_all=[0.9]
#sim_N=10 #单次模拟的源数目
sim_N=100
sim_id_range=np.linspace(1,sim_N,sim_N)
threshold=0.9
sim_id_range=sim_id_range.astype(int)
def get_info(path,period_real):
    detect_rate=[]
    mean=[]
    var=[]
    std=[]
    cts=[]
    for ec_amp in amp_all:
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
                temp_info=np.loadtxt(path+'result_{0}_{1}/result_sim_{3}.txt'.format(str(cts_rate),str(ec_amp),str(width),str(i)))
                cts_get.append(temp_info[-1])
                if temp_info[2]>threshold and 0.99*period_real<temp_info[4]<1.01*period_real \
                    or 1.98*period_real<temp_info[4] <2.01*period_real\
                    or 2.97*period_real<temp_info[4] <3.03*period_real :
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

        return [detect_rate,mean,var,std,cts]

def make_plot():
    plt.figure(1)
    #plt.title('detection')
    detect_rate=get_info(path_1,5258.)[0]
    cts=get_info(path_1,5258.)[-1]
    for i in range(len(cts)):
        plt.plot(cts[i],detect_rate[i],marker='v',color='green')

    detect_rate = get_info(path_2,15258.)[0]
    cts = get_info(path_2,15258.)[-1]
    for i in range(len(cts)):
        plt.plot(cts[i],detect_rate[i],marker='v',color='red')

    detect_rate = get_info(path_3, 45258.)[0]
    cts = get_info(path_3, 45258.)[-1]
    for i in range(len(cts)):
        plt.plot(cts[i], detect_rate[i], marker = 'v', color = 'magenta')

    #plt.legend(['cr=1','cr=2','cr=3','cr=4','cr=5'])
    plt.xlabel('Counts')
    plt.ylabel('Detection rate')
    plt.legend(['P=5258s','P=15258s','P=45258s'])
    print(cts)
    #plt.legend(['width=0.2', 'width=0.3', 'width=0.4', 'width=0.5'],loc='best')
    #plt.savefig('/Users/baotong/Desktop/aas/pCV_in_LW/figure/sim_LW/eclipse.eps')
    #plt.savefig('./fig/'+'detection_{0}.eps'.format(str(period_real)))
    plt.show()
   # plt.close(1)
    #
    # plt.figure(2)
    # #plt.semilogy()
    # plt.title('period_mean')
    # plt.xlabel('cts')
    # plt.ylabel('p-mean')
    # #plt.ylim(period_real*(1-0.04), period_real*(1+0.04))
    # plt.plot([0,10000], np.zeros(2) + period_real, '--')
    # for i in range(len(cts)):
    #     plt.plot(cts[i], mean[i], marker = 'v')
    # plt.ylim(period_real*0.99,period_real*1.01)
    # #plt.legend(['true_p', 'cr=1', 'cr=2', 'cr=3', 'cr=4', 'cr=5'])
    # plt.legend(['true_p','width=0.2', 'width=0.3', 'width=0.4', 'width=0.5'], loc = 'best')
    # plt.savefig('./fig/'+'period_mean_{0}.eps'.format(str(period_real)))
    # plt.show()
    # plt.close(2)

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