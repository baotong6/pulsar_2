# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import string

def get_info_sin(path,period,cts_range,amp_range,sim_N=100,threshold=0.99):
    period_real=period
    sim_id_range=np.arange(1,sim_N+1,1)
    detect_rate=[];mean=[];var=[];std=[]
    for cts_rate in cts_range:
        res_detect = [];res_mean = [];res_var = [];res_std = []
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

def get_LSinfo_sin(path,period,cts_range,amp_range,sim_N=100,threshold=0.9973):
    bin_len=100;ifvary_value=10;varydelta=0
    period_real=period
    sim_id_range=np.arange(1,sim_N+1,1)
    DE=np.zeros((len(cts_range),len(amp_range)))
    for i in range(len(cts_range)):
        cts_rate=cts_range[i]
        for j in range(len(amp_range)):
            amp=amp_range[j]
            path_file=path+'result_{0}_{1}/'.format(cts_rate,amp)
            # LS_info=np.loadtxt(path_file+'sim_LS_bin{0}_VI{1}_N{2}.txt'.format(bin_len,ifvary_value,sim_N))
            LS_info = np.loadtxt(path_file + 'sim_LS_bin{0}_VD{1}_N{2}.txt'.format(bin_len, varydelta, sim_N))
            FP_det=LS_info[:,1];period_det=LS_info[:,2]
            DE[i,j]=len(np.intersect1d(np.where(FP_det<1-threshold),np.where(np.abs(period_det-period_real)<0.1*period_real)))/sim_N
    return DE

def make_plot(detection_rate,cts_range,amp_range,figurepath):
    plt.figure(1)
    # plt.title('detection')
    for i in range(len(cts_range)):
        plt.plot(amp_range,detection_rate[i],marker='v')
    plt.xlabel('Amplitude')
    plt.ylabel('Detection rate')
    plt.legend(title='P={0}s'.format(period_real),labels=['CR=0.5','CR=0.7', 'CR=1.0', 'CR=1.5', 'CR=2.0'],fontsize='small')
    yticks=np.linspace(0,1,11)
    plt.yticks(yticks)
    plt.savefig(figurepath+'Detection_{0}_GL.eps'.format(str(int(period_real))))
    plt.show()

# def estimate_num(DE,)
if __name__ == '__main__':
    path = '/Users/baotong/eSASS/data/raw_data/47_Tuc/simulation/'
    # cts_range = [0.5, 0.7, 1.0,1.5,2.0]
    # amp_range = [0.2, 0.3, 0.4, 0.5,0.6,0.7,0.8]
    cts_range = [0.5, 0.7, 1.0]
    amp_range = [0.2, 0.3, 0.4, 0.5]
    period_real=3600.
    detection_rate1 = get_info_sin(path+'const_CR_GL/', period_real, cts_range, amp_range)[0]
    # detection_rate2 = get_info_sin(path + 'vary_delta_GL/', period_real, cts_range, amp_range)[0]
    # detection_rate3 = get_LSinfo_sin(path + 'vary_delta_LS/', period=3600., cts_range=cts_range, amp_range=amp_range, sim_N=100,
    #                      threshold=0.9973)
    # detection_rate4 = get_LSinfo_sin(path + 'const_delta_LS/period_1h/', period=period_real, cts_range=cts_range, amp_range=amp_range, sim_N=100,
    #                      threshold=0.9973)
    # detection_rate4 = get_LSinfo_sin(path + 'const_delta_25k_LS/period_1h/', period=period_real, cts_range=cts_range, amp_range=amp_range, sim_N=100,
    #                      threshold=0.9973)
    # make_plot(detection_rate1, cts_range, amp_range)
    # make_plot(detection_rate2, cts_range, amp_range)
    make_plot(detection_rate1, cts_range, amp_range,figurepath=path+'fig_sim/')
    # DE1=get_info_sin(path+'const_CR/', period=3600, cts_range=cts_range, amp_range=amp_range, sim_N=100, threshold=0.9973)[0]
    # DE2=get_LSinfo_sin(path+'sim_LS/',period=3600.,cts_range=cts_range,amp_range=amp_range,sim_N=100,threshold=0.9973)
    # DE3 = get_LSinfo_sin(path + 'vary_delta_LS/', period=3600., cts_range=cts_range, amp_range=amp_range, sim_N=100,
    #                      threshold=0.9973)
    # print(DE3)
    # print(DE2)