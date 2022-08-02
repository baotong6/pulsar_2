# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import string
from scipy.interpolate import Rbf,interp2d
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
        plt.plot(amp_range,detection_rate[i],marker='v',label=f'Counts={cts_range[i]*250}')
    plt.xlabel('Amplitude')
    plt.ylabel('Detection rate')
    ##====specific plot====###
    plt.legend(title='P=7200s')
    # x1=100/(np.array(cts_range)*0.01*25000)
    # for i in range(len(x1)):
    #     plt.plot([x1[i],x1[i]],[0,1])
    yticks=np.linspace(0,1,11)
    plt.yticks(yticks)
    # plt.savefig(figurepath+'Detection_{0}_GL.eps'.format(str(int(period_real))))
    plt.show()

def plot_Backlevel():
    ecf=0.75
    path='/Users/baotong/eSASS/data/raw_data/47_Tuc/'
    srcid=np.arange(1,889,1)
    obsIDlist = [700011,700163,700013,700014,700173,700174,700175]
    expT=[25808.86,25268.79,25208.83,25208.82,8809.16,8389.71,8330.68]
    A_all=[];C_all=[]
    for k in range(len(obsIDlist)):
        obsid=obsIDlist[k]
        srcinfo=np.loadtxt(path + 'txt/txt_psf{0}_{1}/src_info.txt'.format(int(100*ecf), obsid))
        C=srcinfo[:,1];B=srcinfo[:,2]
        filter_index=np.loadtxt(path+'txt/inter_srcID.txt')
        filter_index=np.array(filter_index)-1
        filter_index=filter_index.astype('int')
        C=C[filter_index];B=B[filter_index]
        # A=(C-B)/(C+B)
        A=(C-B)/C
        # print(C[195])
        A[np.isnan(A)]=0
        A_all.append(A);C_all.append(C)
    A_all=np.array(A_all);C_all=np.array(C_all)
    A_all_mean=[np.sum(A_all[:,i])/(np.count_nonzero(A_all[:,i])) for i in range(len(A_all[0]))]
    C_all_mean=np.array([np.sum(C_all[:,i])/(np.count_nonzero(C_all[:,i])) for i in range(len(C_all[0]))])
    A_all_mean=np.nan_to_num(A_all_mean)
    plt.hist(A_all_mean,bins=20,histtype='step')
    plt.show()
    plt.hist(C_all_mean, bins=np.logspace(1,5,20), histtype='step')
    print(len(np.where(C_all_mean>100)[0]))
    plt.semilogx()
    plt.show()
    plt.scatter(A_all_mean,C_all_mean)
    plt.show()

def estimate_num(detection_rate,cts_range,amp_range):
    ecf=0.75
    path='/Users/baotong/eSASS/data/raw_data/47_Tuc/'
    srcid=np.arange(1,889,1)
    obsIDlist = [700011,700163,700013,700014,700173,700174,700175]
    expT=[25808.86,25268.79,25208.83,25208.82,8809.16,8389.71,8330.68]
    A_all=[];C_all=[]
    for k in range(len(obsIDlist)):
        obsid=obsIDlist[k]
        srcinfo=np.loadtxt(path + 'txt/txt_psf{0}_{1}/src_info.txt'.format(int(100*ecf), obsid))
        C=srcinfo[:,1];B=srcinfo[:,2]
        filter_index=np.loadtxt(path+'txt/inter_srcID.txt')
        filter_index=np.array(filter_index)-1
        filter_index=filter_index.astype('int')
        C=C[filter_index];B=B[filter_index]
        A=(C-B)/C
        # print(C[195])
        A[np.isnan(A)]=0
        A_all.append(A);C_all.append(C)
    A_all=np.array(A_all);C_all=np.array(C_all)
    A_all_mean=[np.sum(A_all[:,i])/(np.count_nonzero(A_all[:,i])) for i in range(len(A_all[0]))]
    C_all_mean=np.array([np.sum(C_all[:,i])/(np.count_nonzero(C_all[:,i])) for i in range(len(C_all[0]))])
    A_all_mean=np.nan_to_num(A_all_mean)
    # func=interp2d(amp_range,cts_range,detection_rate,kind='cubic')
    amp_range,cts_range=np.meshgrid(amp_range,cts_range,indexing='xy')
    func = Rbf(amp_range, cts_range, detection_rate, function='cubic')
    est_num=0
    print(func)
    net_CTS=[];est_DR_all=[]
    for i in range(len(A_all_mean)):
        est_DR=func(A_all_mean[i],C_all_mean[i]/25000*100)
        if C_all_mean[i]<100 :est_DR=0
        if est_DR<0:est_DR=0
        if est_DR>1:est_DR=1
        if est_DR>0.99:print(C_all_mean[i],A_all_mean[i])
        net_CTS.append(C_all_mean[i]*A_all_mean[i]);est_DR_all.append(est_DR)
        est_num+=est_DR
    # print(est_DR)
    plt.scatter(net_CTS,est_DR_all)
    plt.show()
    print(est_num)
if __name__ == '__main__':
    path = '/Users/baotong/eSASS/data/raw_data/47_Tuc/simulation/'
    # cts_range = [0.5, 0.7, 1.0,1.5,2.0]
    # amp_range = [0.2, 0.3, 0.4, 0.5,0.6,0.7,0.8]
    cts_range = np.array([1.0,2.0,3.0,4.0])
    amp_range = np.array([0.15,0.2,0.25,0.3,0.4,0.5])
    period_real=7200.
    detection_rate1 = get_info_sin(path+'GL_25k/period_2h/', period_real, cts_range, amp_range,threshold=0.9)[0]
    # detection_rate2 = get_info_sin(path + 'vary_delta_GL/', period_real, cts_range, amp_range)[0]
    # detection_rate3 = get_LSinfo_sin(path + 'vary_delta_LS/', period=3600., cts_range=cts_range, amp_range=amp_range, sim_N=100,
    #                      threshold=0.9973)
    # detection_rate4 = get_LSinfo_sin(path + 'const_delta_LS/period_1h/', period=period_real, cts_range=cts_range, amp_range=amp_range, sim_N=100,
    #                      threshold=0.9973)
    # detection_rate4 = get_LSinfo_sin(path + 'const_delta_25k_LS/period_1h/', period=period_real, cts_range=cts_range, amp_range=amp_range, sim_N=100,
    #                      threshold=0.9973)
    # make_plot(detection_rate1, cts_range, amp_range,figurepath=path+'fig_sim/')
    # make_plot(detection_rate2, cts_range, amp_range)
    # make_plot(detection_rate1, cts_range, amp_range,figurepath=path+'fig_sim/')
    # DE1=get_info_sin(path+'const_CR/', period=3600, cts_range=cts_range, amp_range=amp_range, sim_N=100, threshold=0.9973)[0]

    # cts_range = [0.5,0.7,1.0,1.5,2.0]
    # amp_range = np.array([0.7,0.8])
    # DE2=get_LSinfo_sin(path+'const_delta_25k_LS/period_2h/',period=7200.,cts_range=cts_range,amp_range=amp_range,sim_N=100,threshold=0.9973)

    # make_plot(DE2, cts_range, amp_range,None)
    # DE3 = get_LSinfo_sin(path + 'vary_delta_LS/', period=3600., cts_range=cts_range, amp_range=amp_range, sim_N=100,
    #                      threshold=0.9973)
    # print(DE3)
    # detection_rate1=np.array(detection_rate1)
    # estimate_num(detection_rate1, cts_range, amp_range)
    # print(detection_rate1)
    plot_Backlevel()