# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import string

#path_p1='/Users/baotong/Desktop/period_LW/simulation/simulation_LW_554/'
path_p1='/Users/baotong/Desktop/period/simulation/simulation_NSC_I_554/'
path_p2='/Users/baotong/Desktop/period/simulation/simulation_NSC_I_5540/'
#path_p2='/Users/baotong/Desktop/period_LW/simulation/simulation_LW_5540/'
#path_p3='/Users/baotong/Desktop/period_LW/simulation/simulation_LW_45540/'
path_p3='/Users/baotong/Desktop/period/simulation/simulation_NSC_I_45540/'
path_out='/Users/baotong/Desktop/aas/V63/figure/sim_NSC_I/'
#cts_range=[1,2,3,4,5]
#cts_range=[1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,15.0]
# cts_range=[2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9]
#cts_range=[0.5,1.0,2.0,3.0,4.0,5.0]
#amp_range=[0.5,0.6,0.7,0.8,0.9]
cts_range=[1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,15.0]
amp_range=[0.2,0.3,0.4,0.5,0.6,0.7,0.8]
# amp_range=[0.5]
#sim_N=10 #单次模拟的源数目
sim_N=100
threshold=0.9
def get_info(path,period):
    period_real=period
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

                if temp_info[2]>threshold and 0.9*period_real<temp_info[4]<1.1*period_real:
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

    return [detect_rate,mean,var,std]
# print(detect_rate)
# print(mean)

def make_plot():
    period_real=554.
    detect_rate=get_info(path_p1,period_real)[0]
    print(detect_rate)
    os.chdir(path_out)
    plt.figure(1)
    # plt.title('detection')
    for i in range(len(cts_range)):
        plt.plot(amp_range,detect_rate[i],marker='v')

    #plt.legend(['cr=1','cr=2','cr=3','cr=4','cr=5'])
    #plt.ylim(ymax=1.01)
    plt.xlabel('Amplitude')
    plt.ylabel('Detection-rate')
    plt.legend(title='P=554s',labels=['CR=1','CR=2', 'CR=3', 'CR=4', 'CR=5', 'CR=6','CR=7','CR=8','CR=9','CR=15'],fontsize='small')
    # plt.text(0.7,0.4,'P=5540s')
    #plt.savefig('detection.eps')
    yticks=np.linspace(0,1,11)
    plt.yticks(yticks)

    plt.savefig('Detection_{0}.eps'.format(str(int(period_real))))
    plt.show()
make_plot()

def plot_dither_40():
    path='/Users/baotong/Desktop/period_LW/simulation/40_dither/'
    detect=0
    threshold=0.9
    period_real=5074.

    for i in range(1,101):
        temp_info = np.loadtxt(path + 'result_sim_{0}.txt'.format(i))
        if temp_info[2] > threshold and 0.2 * period_real < temp_info[4] < 1.2 * period_real:
            detect+=1

    print(detect)



#plot_dither_40()