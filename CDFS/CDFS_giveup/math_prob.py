#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import re
import string
#import correct as correct
from datetime import datetime
from scipy.interpolate import lagrange
from scipy import interpolate
from scipy.optimize import curve_fit
from scipy.special import comb, perm
import pandas as pd
from astropy.stats import poisson_conf_interval
from astropy.timeseries import LombScargle
import scipy
import stingray as sr
from stingray.events import EventList
from stingray.lightcurve import Lightcurve
from stingray import Lightcurve, Crossspectrum, sampledata,Powerspectrum,AveragedPowerspectrum
from stingray.simulator import simulator, models
import linecache
figurepath='/Users/baotong/Desktop/aas/AGN_CDFS/figure/'
font1 = {'family': 'Normal',
         'weight': 'normal',
         'size': 16, }
#3.13日：NB,两个概率都算对了#
#3.23日：这俩概率好像没啥用#
def Pi(i,N,D,f):
    # 给定样本N，探测效率D，误报率f
    # 返回值为一个长度为i+1的数组a
    # a[j]为，样本N中探测到i个信号，但其中j个为假信号的可能性
    a=np.zeros(i+1);j=0
    while j<(i+1):
        a[j]=(1-D-f)**(N-i)*comb(N,i)*(comb(i,i-j)*D**(i-j)*f**j)
        j+=1
    return a

def Qi(i,N,Nq,D,f1,f2):
    # 给定样本N，其中intrinsic有QPO的样本为Nq，对应的探测效率为D，误报率为f1；
    # 空白样本的探测效率为f2
    # 返回值应为一个长度为i+1的数组bi
    # bi[j]为，样本N中探测到i个信号，但其中j个为假信号的可能性
    # i-j应当小于等于Nq;i-k应当小于等于Nq
    bi=np.zeros(i+1);j=np.max([0,i-Nq])
    if (i-Nq)>0:
        bi[0:i-Nq]=0
    ##这个步骤并没有用，只是强调一下这些数值为0
    while j<(i+1):
        k=np.max([0,i-Nq]);
        while k<np.min([j+1,N-Nq+1]):
            Qijk=Pi(i-k,Nq,D,f1)[j-k]*comb(N-Nq,k)*f2**k*(1-f2)**(N-Nq-k)
            bi[j]+=Qijk
            k+=1
        j+=1
    return bi

def get_mostP_Nq(i,N,D,f1,f2):
    # 已知样本N,探测到了i个信号
    # 返回最有可能的Nq值
    # 其中Pq为N=(Nq,N-Nq)样本下，探测到i个信号的可能性
    Pq=np.zeros(N+1)
    Nq=0
    while Nq<(N+1):
        Pq[Nq]=np.sum(Qi(i,N,Nq,D,f1,f2))
        Nq+=1
    plt.scatter(np.linspace(0,N,N+1),Pq)
    plt.show()
# get_mostP_Nq(2,43,0.00897,0.0047,6e-4)

def get_sim_Prob(k,threshold):
    path = '/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep{0}/simulation/'.format(k)
    cts_rate_str = ['2e-5', '5e-5', '1e-4', '2e-4', '3e-4', '4e-4','5e-4']
    CR=[4e-5,1e-4,2e-4,4e-4,6e-4,8e-4,1e-3]

    true_period = [1/ 2.777e-3,1/1.38889e-3,1 / 6.68e-4, 1 / 2.7e-4]
    period_str = ['0.1hr','0.2hr','0.4hr', '1hr']
    DR=np.zeros((len(cts_rate_str),len(period_str)))
    fDR=np.zeros((len(cts_rate_str),len(period_str)))
    fDR_noQPO=np.zeros(len(cts_rate_str))
    for i in range(len(cts_rate_str)):
        noQPOfile=np.loadtxt(path+'trial_out_{0}_REJ1034+396_noQPO.txt'.format(cts_rate_str[i]))
        fp_noQPO=noQPOfile[:,0];period_noQPO=noQPOfile[:,1]
        fDR_noQPO[i] = len(np.where((fp_noQPO < threshold))[0])
        # fDR_noQPO[i]=len(np.where((fp_noQPO<threshold)&(period_noQPO>210.))[0])
        for j in range(len(period_str)):
            file=np.loadtxt(path+'trial_out_{0}_{1}_REJ1034+396.txt'.format(period_str[j],cts_rate_str[i]))
            fp = file[:,0]
            period_det=file[:,1]
            DR[i][j]=len(np.intersect1d(np.where(fp < threshold)[0],np.where(np.abs(1 / period_det - 1/true_period[j]) < (1/(true_period[j]*10)))[0]))
            fDR[i][j]=len(np.intersect1d(np.where(fp < threshold)[0],np.where(np.abs(1 / period_det - 1/true_period[j]) > (1/ (true_period[j]*10)))[0],np.where(period_det<10000)))
    return [DR,fDR,fDR_noQPO]

def plot_sim_DR(k_num,threshold):
    figlabel = [[0, 0], [0, 1], [1, 0], [1, 1]]
    fig, axes = plt.subplots(2, 2,figsize=(15,10))
    for i in range(len(k_num)):
        k = k_num[i];
        CDFS_LS_res = np.loadtxt('/Users/baotong/Desktop/CDFS/fig_LS_ep{0}_ovsamp_5_baluev/LS_result_{0}.txt'.format(k))
        bins_CR = np.array([0, 4e-5, 7e-5, 1.5e-4, 2.5e-4, 3.5e-4, 4.5e-4, 5.5e-4]) * 2
        plt.figure(2)
        CRhist = plt.hist(CDFS_LS_res[:, 5], bins=bins_CR, histtype='step')
        # print(CRhist)
        print(len(np.where(CDFS_LS_res[:, 5] > 4e-4)[0]))
        # 计算大于某个流量的源有多少
        plt.close(2)

        ax_temp = axes[figlabel[i][0], figlabel[i][1]]
        [DR, fDR, fDR_noQPO]=get_sim_Prob(k,threshold)


        # DR/=1000;fDR/=1000
        x=np.array([4e-5,1e-4,2e-4,4e-4,6e-4,8e-4,1e-3])*1e5
        y1=DR[:,0];y2=DR[:,1];y3=DR[:,2];y4=DR[:,3]
        z1=fDR[:,0];z2=fDR[:,1];z3=fDR[:,2];z4=fDR[:,3]

        m=fDR_noQPO

        sim_NUM_04h = np.sum(CRhist[0] * y3/ 1000)
        print(CRhist[0] * y3/ 1000)
        print('sim_NUM_04h={0}'.format(sim_NUM_04h))
        sim_NUM_false=np.sum(CRhist[0] * m/ 1000)
        print('sim_NUM_false={0}'.format(sim_NUM_false))
        # sim_real_NUM = np.sum(CRhist[0] * DR_04hr / 1000)

        # ax_temp.set_yscale('log')
        ax_temp.plot(x,y1/10,marker='v', linestyle='-', color='blue')
        ax_temp.plot(x,y2/10,marker='v', linestyle='-', color='red')
        ax_temp.plot(x,y3/10,marker='v', linestyle='-', color='green')
        ax_temp.plot(x,y4/10,marker='v', linestyle='-', color='brown')

        ax_temp.plot(x,m/10,marker='o',linestyle='--',color='grey')
        # ax_temp.plot([0,100],[0,0],'--',color='yellow')
        # ax_temp.plot(x,z1,marker='o', linestyle='--', color='blue')
        # ax_temp.plot(x,z2,marker='o', linestyle='--', color='red')
        # ax_temp.plot(x,z3,marker='o', linestyle='--', color='green')
        # ax_temp.plot(x,z4,marker='o', linestyle='--', color='brown')
        #
        ax_temp.legend(['DR P=0.1h','DR P=0.2h','DR P=0.4h','DR P=1h','fDR (no QPO input)'])
        # ax_temp.legend(['DR P=0.1h','DR P=0.2h','DR P=0.4h','DR P=1h','FDR P=0.1h','FDR P=0.2h','FDR P=0.4h','FDR P=1h'])
        ax_temp.text(2,10,'Epoch{0}'.format(k),font1)
        if i<3:ax_temp.legend_.remove()
        # ax_temp.set_title('Epoch {0}: LS detection results'.format(k),font1)
        if (i==2 or i==3):ax_temp.set_xlabel('Photon flux ($10^{-5}$ counts $s^{-1}$ )',font1)
        if (i==0 or i==2):ax_temp.set_ylabel('Detection rate(%)',font1)
        ax_temp.tick_params(labelsize=16)
    # plt.subplots_adjust(wspace=0, hspace=0)
    plt.savefig(figurepath+'DR_thres_{0}.eps'.format(1-threshold),bbox_inches='tight',pad_inches=0.0)

    plt.show()

plot_sim_DR([1,2,3,4],1-0.99)

def read_LS_info(k,threshold):
    path='/Users/baotong/Desktop/CDFS/fig_LS_ep{0}_ovsamp_5_baluev/'.format(k)
    LS_res=np.loadtxt(path+'LS_result_{0}.txt'.format(k))
    LS_res=LS_res[np.where((LS_res[:,5]>2e-4)&(LS_res[:,5]<1e-3))]
    src_id = LS_res[:, 0];conf = LS_res[:, 1];period = LS_res[:, 2];
    src_cts=LS_res[:,3];bkg_cts=LS_res[:,4];cr = LS_res[:, 5]

    CR=np.array([4e-5,1e-4,2e-4,4e-4,6e-4,8e-4,1e-3])
    [DR,fDR,fDR_noQPO]=get_sim_Prob(k,threshold)
    DR_weight_P=DR.mean(axis=1);fDR_weight_P=fDR.mean(axis=1);
    #只能假设光子数足够多的时候，探测效率为100%
    inter_DR = interpolate.interp1d(np.concatenate(([0],CR,[1e-2])), np.concatenate(([0],DR_weight_P,[300])), kind='cubic')
    inter_fDR = interpolate.interp1d(np.concatenate(([0],CR,[1e-2])), np.concatenate(([0],fDR_weight_P,[5])), kind='cubic')
    inter_fDR_noQPO = interpolate.interp1d(np.concatenate(([0],CR)), np.concatenate(([0],fDR_noQPO)), kind='cubic')
    plt.plot(CR,DR_weight_P/1000,color='purple')
    plt.scatter(cr,inter_DR(cr)/1000,color='black')
    plt.plot(CR,fDR_weight_P/1000,color='red')
    plt.scatter(cr,inter_fDR(cr)/1000,color='green')
    plt.plot(CR,fDR_noQPO/1000,color='grey')

    return [np.mean(inter_DR(cr))/1000,np.mean(inter_fDR(cr))/1000,np.mean(fDR_noQPO)/1000]
#
# RES=read_LS_info(3,1-0.991222)
# print(RES)
# a=Pi(2,43,RES[0],RES[1])
# print(a)
# print(1-a[-1]/np.sum(a))
# print(a[0]/np.sum(a))
# print(get_mostP_Nq(2,43,RES[0],RES[1],RES[2]))
def plot_scatter_bootstrap(k,srcid):
    path='/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep{0}/simulation/'.format(k)
    res=np.loadtxt(path+'trial_out_src{0}_REJ1034+396_noQPO.txt'.format(srcid))
    withQPO_res=np.loadtxt(path+'trial_out_src{0}_REJ1034+396.txt'.format(srcid))

    fig, axes = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(15, 10))
    ax_temp = axes[0]
    ax_temp.scatter(withQPO_res[:,1],withQPO_res[:,0])
    ax_temp.plot([0,20000],[1-0.9912,1-0.9912],'--')
    ax_temp.set_xscale('log')
    ax_temp.set_yscale('log')
    print('DR={0}'.format(len(np.where((withQPO_res[:,0]<1-0.99)&(np.abs(1/withQPO_res[:,1]-1/950.73162)<1/(10*950.73162)))[0])))
    print('fDR={0}'.format(len(np.where((withQPO_res[:,0]<1-0.99)&(np.abs(1/withQPO_res[:,1]-1/950.73162)>1/(10*950.73162)))[0])))
    ax_temp = axes[1]
    ax_temp.plot([0, 20000], [1 - 0.9912, 1 - 0.9912], '--')
    ax_temp.scatter(res[:,1],res[:,0])
    ax_temp.set_xscale('log')
    ax_temp.set_yscale('log')
    print('blankFD={0}'.format(len(np.where((res[:,0]<1-0.99)&(np.abs(1/res[:,1]-1/950.73162)<1/(10*950.73162)))[0])))
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.show()
# plot_scatter_bootstrap(3,'89')
# read_LS_power.read_LSP('3','89',num_trials=10000)

def plot_LS_bootstrapFAP(k,srcid,threshold):
    path='/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep{0}/simulation/'.format(k)
    time=np.loadtxt('/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep{0}/{1}.txt'.format(k,srcid))[:,0]
    bkg_time=np.loadtxt('/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep{0}/{1}_bkg.txt'.format(k,srcid))[:,0]
    T_exp = time[-1] - time[0];
    T_exp=11000154.981141508
    dt=100
    freq = np.arange(1 / T_exp, 0.5 / dt, 1 / (5 * T_exp))
    freq = freq[np.where(freq > 1 / 20000.)]
    # print(freq[55100])
    peakP_index=55100
    ##即在这个频率取样下第55100个是周期
    peakfreq=[1/950.73161926-1/(T_exp),1/950.73161926+1/(T_exp)]
    # peak_sim=np.zeros(len(freq)-1)
    # x=[];y=[]
    # for i in range(len(peak_sim)):
    #     peak_sim[i]=np.mean(peakP[np.where((period>1/freq[i+1])&(period<1/freq[i]))])
    #     if np.isnan(peak_sim[i])==False:x.append(freq[i]);y.append(peak_sim[i])

    simP=np.loadtxt('/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep{0}/simulation/{1}_LS_sim_noQPO/LSP_{2}.txt'.format(k,srcid,threshold))
    x=simP[:,0];y=simP[:,1]*1e-2
    x=freq   ## 这是因为savetxt的时候保存精度不够，所以用freq
    def get_smoothxy(x,y,smooth):
        i=0;smoothx=np.zeros(len(freq));smoothy=np.zeros(len(y))
        while i <len(freq):
            smoothy[i:i+smooth]=np.max(y[i:i+smooth])
            smoothx[i:i+smooth]=freq[i]
            i+=smooth
        return (smoothx,smoothy)
    (smoothx, smoothy)=get_smoothxy(x,y,smooth=1)
    plt.figure(1,(10,7))
    # plt.step(smoothx,smoothy,color='grey',linewidth=0.2)
    # plt.text(5e-5,0.19,'99.95%',color='r',weight=2,fontsize=14)
    # plt.plot([0,20000],[0.00878,0.00878],'--')

    def get_max_trial(LSP,threshold=0.9973,num_trials=10000):
        path='/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep{0}/simulation/{1}_LS_sim_noQPO/'.format(k,srcid)
        maxP=np.loadtxt(path+'maxP_trial_const.txt')
        idex = np.lexsort([maxP[:,0]])
        sorted_data = maxP[idex, :]
        sort_maxP=sorted_data[:,0]
        freq_sortP=sorted_data[:,1]
        bins=np.logspace(np.log10(100),np.log10(20000),15)
        freq_out_overLSP=freq_sortP[np.where(sort_maxP>LSP)]
        print(len(freq_out_overLSP))
        freq_out=freq_sortP[int(num_trials*threshold):]
        plt.title('const Distribution of false detection: XID={0}'.format(srcid))
        plt.hist(1/freq_out,bins=bins,histtype='step',linewidth=5,color='green')
        plt.semilogx()
        plt.xlabel('Period')
        plt.ylabel('Num')
        plt.show()

        return sort_maxP[int(num_trials*threshold)]*0.01

    maxP=get_max_trial(LSP=12.43)
    print(maxP)
    plt.plot([freq[0],freq[-1]],[maxP,maxP],'--',color='orange')
    plt.text(freq[0],maxP*1.01,'99.73%')


    def get_hist_withbkg(t, t_bkg, len_bin):
        ###将输入的time信息，按照len_bin的长度输出为lc
        t_test = t - t[0];
        t_bkg_test = t_bkg - t[0];
        dt = len_bin
        t_bkg_test = np.delete(t_bkg_test, t_bkg_test < 0)
        ev = EventList();
        ev_bkg = EventList()
        ev.time = t_test;
        ev_bkg.time = t_bkg_test
        lc_new = ev.to_lc(dt=dt, tstart=ev.time[0] - 0.5 * dt, tseg=ev.time[-1] - ev.time[0])
        lc_bkg = ev_bkg.to_lc(dt=dt, tstart=ev.time[0] - 0.5 * dt, tseg=ev.time[-1] - ev.time[0])
        lc_out = lc_new
        lc_out.counts = lc_new.counts - (1 / 12.) * lc_bkg.counts
        return lc_out

    def get_LS(time, flux, freq, dataname, k):
        x = time
        y = flux
        # dy=np.sqrt(y)
        # plt.scatter(x,y)
        # plt.show()

        # LS = LombScargle(x, y, dy = 1, normalization = 'standard', fit_mean = True,
        #                  center_data = True).power(freq, method = 'cython')
        LS = LombScargle(x, y, normalization='psd')
        # LS = LombScargle(x, y, dy, normalization='psd')
        power = LS.power(freq)
        FP=1
        # FP = LS.false_alarm_probability(power.max(), minimum_frequency=freq[0], maximum_frequency=freq[-1],
        #                                 method='baluev')
        # FP_99 = LS.false_alarm_level(0.0027, minimum_frequency=freq[0], maximum_frequency=freq[-1], method='baluev')
        # FP_90 = LS.false_alarm_level(0.05, minimum_frequency=freq[0],
        #                              maximum_frequency=freq[-1], method='baluev')
        # FP_68 = LS.false_alarm_level(0.32, minimum_frequency=freq[0],
        #                              maximum_frequency=freq[-1], method='baluev')

        # if FP < 0.01: print(dataname)
        # plt.title('{0}, FP={1}%'.format(dataname, np.round(100*FP,3)))
        # plt.semilogx()
        plt.title(' XID={0}'.format(dataname),font1)
        plt.plot(freq, power)
        # plt.loglog()
        plt.xlabel('Frequency (Hz)',font1)
        plt.ylabel('LS Periodogram',font1)
        plt.tick_params(labelsize=16)
        plt.semilogx()
        # print(1. / freq[np.where(power == np.max(power))])
        # if FP < 0.01: print(1. / freq[np.where(power == np.max(power))])
        print(np.max(power))
        out_period = 1. / freq[np.where(power == np.max(power))]
        # plt.plot([freq[0], freq[-1]], [FP_99, FP_99], '--')
        # plt.plot([freq[0], freq[-1]], [FP_90, FP_90], '--')
        # plt.plot([freq[0], freq[-1]], [FP_68, FP_68], '--')
        # plt.savefig('/Users/baotong/Desktop/CDFS/fig_LS_ep{0}_ovsamp_5_baluev/{1}.eps'.format(k,dataname))
        # plt.close()
        return [FP, out_period]

    # lc = get_hist_withbkg(time, bkg_time, dt)
    # x = lc.time;
    # flux = lc.counts
    # (FP, out_period) = get_LS(x, flux, freq, str(srcid), k)
    # plt.savefig(figurepath+'LS_{0}_simFAP_{1}.eps'.format(srcid,threshold),bbox_inches='tight',pad_inches=0.0)
    # plt.show()
# plot_LS_bootstrapFAP('3','89',threshold=0.9999)