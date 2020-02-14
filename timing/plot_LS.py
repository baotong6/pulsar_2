#!/bin/bash
# -*- coding: utf-8 -*-
#plt.rcParams['font.sans-serif'] = ['Microsoft YaHei']
#用于解决matplotlib中文乱码问题
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
import read_data as data
#import read_data_gc as data
#import read_data_lw as data
from scipy.fftpack import fft,ifft
import scipy.signal as ss
import scipy.stats as stats
import random
from astropy.timeseries import LombScargle

#path='/Users/baotong/Desktop/period_gc/'
#import fits_to_txt
# standard_name='pwn.txt'
#dataname="Muno2.txt"
def get_fig_LS(dataname):
    # def sim_time(dataname):
    #     time = data.get_data(dataname)[0]
    #     energy = data.get_data(dataname)[1]
    #     obs_ID = data.get_data(dataname)[2]
    #     t = data.get_data(dataname)[3]
    #     E = data.get_data(dataname)[4]
    #     dict = data.get_data(dataname)[5]
    #
    #     t=np.array(t)
    #
    #     sim_t=t-t
    #     for i in range(len(t)):
    #
    #         sim_t[i]=np.random.random(len(t[i]))*(t[i][-1]-t[i][0])+t[i][0]
    #
    #     return sim_t

    def get_LS(freq,len_bin,dataname):
        time=data.get_data(dataname)[0]
        energy=data.get_data(dataname)[1]
        obs_ID=data.get_data(dataname)[2]
        t=data.get_data(dataname)[3]
        E=data.get_data(dataname)[4]
        dict=data.get_data(dataname)[5]
        N=len(time)

        ### 选用某些符合需求的obs ###
        delete_i=[]
        use_ID=[]
        use_ID_o= []
        for i in dict:
            if dict[i][-2] < 5000:
                delete_i.append(i)
        for item in delete_i:
            dict.pop(item)
        for i in dict:
            use_ID_o.append(i)

        # epoch0=np.array([2943,3663])
        #epoch1=np.array([2943,3663,3392,3393,3665])
        # epoch2=np.array([5950,5951,5952,5953,5954])
        # epoch3=np.array([9169,9170,9171,9172])
        #
        # epoch1 = np.array([13847, 14427, 13848, 13849, 13846, 14438, 13845])
        # epoch2 = np.array([14461, 13853, 13841, 14465, 14466, 13842, 13839, 13840, 14432, 13838, 13852, 14439])
        # epoch3 = np.array([14462, 14463, 13851, 15568, 13843, 15570, 14468])

        #epoch_gc=np.array([945,14897,17236,17239,17237,18852,17240,17238,20118,17241,20807,20808])
        #epoch_gc_1 = np.array([17236, 17239, 17237, 18852, 17240, 17238])

        #use_ID_o=epoch1

        use_ID = use_ID_o
        not_in_id = []
        for i in range(len(use_ID_o)):
            if use_ID_o[i] in dict:
                continue
            else:
                not_in_id.append(i)
        use_ID = np.delete(use_ID_o, not_in_id)
        if len(use_ID) < 2:
            return [['none']]
        else:
            time = np.concatenate((dict[use_ID[0]][0], dict[use_ID[1]][0]))
            for i in range(2, len(use_ID)):
                time = np.concatenate((time, dict[use_ID[i]][0]))


        #use_ID=[14897,17236,17239,17237,18852,17240,17238,20118,17241,20807,20808]

        # use_ID=np.concatenate((epoch1,epoch2,epoch3))
        time=np.concatenate((dict[use_ID[0]][0],dict[use_ID[1]][0]))
        for i in range(2,len(use_ID)):
            time=np.concatenate((time,dict[use_ID[i]][0]))

        ##确立第一次观测的开始时间为start time
        ##并将所有的时间调至这个零点
        obs_d=[dict[i][-1] for i in dict]
        start_time=obs_d[0][0]
        time -= start_time
        for i in range(len(obs_d)):
            obs_d[i]=obs_d[i]-start_time

        def get_hist(t,len_bin):
            ###将输入的time信息，按照len_bin的长度输出为lc
            t_test=t
            a=[0 for i in range(int(obs_d[-1][1]/len_bin)+1)]
            for i in range(len(t_test)):
                a[int(t_test[i]/len_bin)]+=1
            a=np.array(a)
            return a
        def get_window(t,len_bin,obs_d):
            ###给定观测信息，输出窗函数，即在观测时间内每个bin的值为1，非观测时间为0
            t_test=t
            b=[0 for i in range(int(obs_d[-1][1]/len_bin)+1)]
            for i in range(len(obs_d)):
                temp=int(obs_d[i][0])/len_bin
                while temp<int(obs_d[i][1]/len_bin):
                    b[int(temp)]+=1
                    temp+=1
            return b
        cts = len(time)
        y=get_hist(time,len_bin)
        window=get_window(time,len_bin,obs_d)
        x=np.arange(obs_d[0][0],obs_d[-1][1],len_bin)
        # plt.scatter(x,window)
        # plt.show()
        print('run1')
        LS=LombScargle(x,y,dy=1,normalization='standard',fit_mean=True, center_data=True).power(freq,method='cython')
        LS_W=LombScargle(x,window,normalization='standard',fit_mean=False, center_data=False).power(freq,method='cython')
        FP_99=LombScargle(x,y,dy=1,normalization='standard',fit_mean=True, center_data=True).false_alarm_level(0.01,
                                                                                                               minimum_frequency=freq[0],maximum_frequency=freq[-1])
        FP_90=LombScargle(x,y,dy=1,normalization='standard',fit_mean=True, center_data=True).false_alarm_level(0.1,
                                                                                                               minimum_frequency=freq[0],maximum_frequency=freq[-1])
        FP_68=LombScargle(x,y,dy=1,normalization='standard',fit_mean=True, center_data=True).false_alarm_level(0.32,
                                                                                                               minimum_frequency=freq[0],maximum_frequency=freq[-1])

        def get_LS_sim(dataname,len_bin):
            sim_t = sim_time(dataname)
            simtime = data.tran_t_time(sim_t)
            simtime = np.sort(simtime)
            y = get_hist(simtime, len_bin)
            y = np.array(y)
            x=np.arange(simtime[0],simtime[-1],len_bin)
            LS_sim=LombScargle(x,y,dy=1,normalization='psd',fit_mean=True,center_data=True).power(freq,method='cython')

            return LS_sim
        LS_sim=[]
        cts_window=sum(window)
        return [LS,LS_W,LS_sim,cts,cts_window,[FP_99,FP_90,FP_68]]

    freq_1=[1/50000.,1e-8,3000]
    #freq_1=[1e-3,1e-5,1000]
    freq=freq_1[0]+freq_1[1]*np.arange(freq_1[2])
    len_bin=50

    ### uneven sampled freq ###
    exptime=131610914.53
    p_unsamp = []
    freq_unsamp = []
    p_unsamp.append(1e-6 * exptime)
    while p_unsamp[-1] < exptime:
        if p_unsamp[-1] < 0.3 * exptime:
            p_unsamp.append(p_unsamp[-1] + p_unsamp[-1] ** 2 / (2 * exptime * 2))
        elif 0.3 * exptime < p_unsamp[-1] < 0.5 * exptime:
            p_unsamp.append(p_unsamp[-1] + p_unsamp[-1] ** 2 / (3 * exptime * 2))
        else:
            p_unsamp.append(p_unsamp[-1] + p_unsamp[-1] ** 2 / (4 * exptime * 2))
    p_unsamp = np.array(p_unsamp)
    freq_unsamp = 1.0 / p_unsamp
    ### uneven sampled freq ###

    def make_period_range(pmin, pmax, expT):
        P = [pmin]
        while P[-1] < pmax:
            dP = 0.1 * P[-1] ** 2 / (expT - P[-1])
            P.append(P[-1] + dP)
        return np.array(P)

    freq=1./make_period_range(3000,6000,exptime)


    res_normal=get_LS(freq,len_bin=64,dataname=dataname)
    LS=res_normal[0]
    LS_W=res_normal[1]
    cts=res_normal[3]
    cts_W=res_normal[4]
    FP_all=res_normal[5]

    # res_standard=get_LS(freq,len_bin=100,dataname=standard_name)
    # if res_standard[0][0]=='none':
    #     print('none')
    #     return 'none'
    # LS_standard=res_standard[0]
    # cts_standard=res_standard[3]
    #
    #
    # LS_standard=LS_standard*cts**2/cts_standard**2
    LS_W=LS_W*cts**2/cts_W**2

    def get_conf90(Np,Ns):
        return(-np.log(1-0.99**(1.0/(Ns*Np))))

    conf90 = get_conf90(len(p_unsamp), 50)
    print('running')
    plt.figure(1,(9,6))

    #plt.subplot(221)
    plt.title(dataname[0:-4]+',cts={0}'.format(cts))
    plt.semilogx()
    plt.step(freq,LS)
    print(max(LS))
    print()
    plt.plot([freq[0], freq[-1]], [FP_all[0], FP_all[0]], '--')
    plt.plot([freq[0], freq[-1]], [FP_all[1], FP_all[1]], '--')
    plt.plot([freq[0], freq[-1]], [FP_all[2], FP_all[2]], '--')
    #plt.savefig('/Users/baotong/Desktop/li_pr/result_final/fig_NSC_I/{0}_LS.eps'.format(dataname[0:-4]))

    # plt.subplot(222)
    # plt.title('window_function')
    # plt.semilogx()
    # plt.step(freq,LS_W)
    #
    # plt.subplot(223)
    # plt.title(standard_name[0:-4]+',cts={0}'.format(cts_standard))
    # plt.semilogx()
    # plt.step(freq,LS_standard)
    #
    # plt.subplot(224)
    # plt.title('minus')
    # plt.semilogx()
    # plt.step(freq,LS-LS_standard)
    # plt.plot([freq[0], freq[-1]], [FP_all[0], FP_all[0]], '--')
    # plt.savefig(path+'fig_ep2_200_30k/'+dataname[0:-4]+'.jpeg')
    plt.show()
    # plt.close()


def light_curve(dataname):
    time = data.get_data(dataname)[0]
    energy = data.get_data(dataname)[1]
    obs_ID = data.get_data(dataname)[2]
    t = data.get_data(dataname)[3]
    E = data.get_data(dataname)[4]
    dict = data.get_data(dataname)[5]

    cts_rate_all=[]
    exp_time_all=[]
    cts_all=[]
    x_lc=[0]
    y_lc=[]
    text_obsID=[]
    for i in dict:
        cts_rate_all.append(len(dict[i][0])/dict[i][2])
        cts_all.append(len(dict[i][0]))
        exp_time_all.append(dict[i][2])
        y_lc.append(len(dict[i][0])/dict[i][2])
        y_lc.append(len(dict[i][0])/dict[i][2])
        x_lc.append(x_lc[-1]+dict[i][2])
        text_obsID.append(i)
    for i in range(1,2*len(x_lc)-1,2):
        x_lc.insert(i,x_lc[i])
    x_lc=x_lc[:-1]
    print(len(y_lc))

    N = len(time)
    delete_i = []
    use_ID = []
    for i in dict:
        if dict[i][-2]<5000:
            delete_i.append(i)
    for item in delete_i:
        dict.pop(item)
    for i in dict:
        use_ID.append(i)
    print(len(use_ID))
    #use_ID=[4683,4684]#,4683,4684,5950]
    epoch1=np.array([2943,3663,3392,3393,3665])
    epoch2=np.array([5950,5951,5952,5953,5954])
    epoch3=np.array([9169,9170,9171,9172])
    #use_ID=epoch1
    #use_ID=np.concatenate((epoch1,epoch2,epoch3))
    not_in_id = []
    for i in range(len(use_ID)):
        if use_ID[i] in dict:
            continue
        else:
            not_in_id.append(i)
    use_ID = np.delete(use_ID, not_in_id)

    if len(use_ID) < 1:
        return [0, 0]
    elif len(use_ID) == 1:
        time = dict[use_ID[0]][0]
    else:
        time = np.concatenate((dict[use_ID[0]][0], dict[use_ID[1]][0]))
        for i in range(2, len(use_ID)):
            time = np.concatenate((time, dict[use_ID[i]][0]))
    cts_rate_use = []
    cts_use=[]
    exp_time_use = []
    x_lc_use = [0]
    y_lc_use = []
    use_ID=np.array(use_ID)
    for i in range(len(use_ID)):
        cts_rate_use.append(len(dict[use_ID[i]][0]) / dict[use_ID[i]][2])
        cts_use.append(len(dict[use_ID[i]][0]))
        exp_time_use.append(dict[use_ID[i]][2])
        y_lc_use.append(len(dict[use_ID[i]][0]) / dict[use_ID[i]][2])
        y_lc_use.append(len(dict[use_ID[i]][0]) / dict[use_ID[i]][2])
        x_lc_use.append(x_lc_use[-1] + dict[use_ID[i]][2])
    for i in range(1, 2 * len(x_lc_use) - 1, 2):
        x_lc_use.insert(i, x_lc_use[i])
    x_lc_use = x_lc_use[:-1]

    print(len(y_lc))
    plt.figure(1,(12,9))
    plt.subplot(211)
    plt.semilogy()
    plt.title('lc for source {0}'.format(dataname[0:-4]))
    plt.ylabel('cts rate')
    plt.xlabel('obs time')
    plt.step(x_lc,y_lc)
    # plt.subplot(222)
    # plt.title('used obs')
    # plt.step(x_lc_use,y_lc_use)

    plt.subplot(212)
    #plt.title('lc for source {0}'.format(dataname[0:-4]))
    plt.ylim(1e-5,1e-2)
    plt.semilogy()
    mean_ctsr=np.mean(cts_rate_all)
    for i in range(len(cts_rate_all)):
        plt.annotate(text_obsID[i],xy=(i+1,cts_rate_all[i]),xytext=(i+1,cts_rate_all[i]),weight='light')
        plt.ylabel('cts rate')
        plt.xlabel('obs ID')
        plt.scatter(i+1,cts_rate_all[i],color='red')
        plt.errorbar(i + 1, cts_rate_all[i], yerr=cts_all[i]**0.5/ exp_time_all[i])

    # plt.subplot(224)
    # plt.title('used obs')
    # for i in range(len(cts_rate_use)):
    #     plt.annotate(str(use_ID[i]), xy=(i+1, cts_rate_use[i]), xytext=(i + 1, cts_rate_use[i]))
    #     plt.scatter(i+1,cts_rate_use[i],color='red')
    #     plt.errorbar(i+1,cts_rate_use[i],yerr=cts_use[i]**0.5/exp_time_use[i])

    y_lc=np.array(y_lc)
    print(np.where(y_lc==max(y_lc))[0]/2)
    print(text_obsID[int(np.where(y_lc==max(y_lc))[0][0]/2)])
    plt.savefig('/Users/baotong/Desktop/li_pr/result_final/fig_LW/{0}_lc.eps'.format(dataname[0:-4]))
    plt.close()
    #plt.show()


# cand_ID=fits_to_txt.cand_ID
# for item in cand_ID[0:500]:
#     get_fig_LS(str(item)+'.txt',standard_name)

# cand_ID=np.arange(1,518)
# for item in cand_ID:
#     get_fig_LS(str(item) + '.txt', standard_name)

import pandas as pd
path_table='/Users/baotong/Desktop/period/table/'
result_NSC=pd.read_excel(path_table+'final_all_del.csv','result_NSC')
result_LW=pd.read_excel(path_table+'final_all_del.csv','result_LW')
result_ND=pd.read_excel(path_table+'final_all_del.csv','result_ND')

ID_NSC=result_NSC['seq']
ID_LW=result_LW['seq']
ID_ND=result_ND['seq']
# item='1677'
for item in ID_LW:
#get_fig_LS(str(item)+'.txt')
    light_curve(str(item)+'.txt')