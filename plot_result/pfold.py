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
from scipy.optimize import curve_fit
import pandas as pd
from astropy.stats import poisson_conf_interval
import scipy
# import mpmath
#import read_csv as data

font1 = {'family': 'Normal',
         'weight': 'normal',
         'size': 18, }

###----------------read table---------------###
path_table = '/Users/baotong/Desktop/period/table/'
result_NSC = pd.read_excel(path_table + 'final_all_del.csv', 'result_NSC')
result_LW=pd.read_excel(path_table+'final_all_del.csv','result_LW')
result_ND=pd.read_excel(path_table+'final_all_del.csv','result_ND')
result_NSC_IG=pd.read_excel(path_table + 'final_all_del.csv', 'result_NSC_IG')

ID_NSC_IG=result_NSC_IG['seq']
P_NSC_IG=result_NSC_IG['P']
ra_NSC_IG=result_NSC_IG['ra']
dec_NSC_IG=result_NSC_IG['dec']
net_percent_NSC_IG=result_NSC_IG['net_percent']
###----------------read table---------------###

###----------------read table---------------###
path_table = '/Users/baotong/Desktop/period_Tuc/'
result_Tuc = pd.read_excel(path_table + 'result_all.csv', '47Tuc')

ID_Tuc=result_Tuc['seq']
P_Tuc=result_Tuc['P']
# ra_NSC_IG=result_NSC_IG['ra']
# dec_NSC_IG=result_NSC_IG['dec']
net_percent_Tuc=result_Tuc['net_percent']
###----------------read table---------------###

def filter_energy(time,energy,band):
    T=time
    E=energy
    i=0
    if len(time)!=len(energy):
        print('error')
        return None
    else:
        while i <len(E):
            if E[i]<band[0] or E[i]>band[1]:
                E=np.delete(E,i)
                T=np.delete(T,i)
            else:
                i+=1

    return T

def get_T_in_mbins(epoch_info,w,m,fi):
    T=2*np.pi/w
    T_in_perbin = np.zeros(m)
    # 每个bin的总积分时间
    tbin = T/m
    # 每个bin的时间长度
    t_start=epoch_info[:,0];t_end = epoch_info[:, 1]
    N_bin_t_start=t_start/tbin+m*fi/(2*np.pi)
    N_bin_t_end=t_end/tbin+m*fi/(2*np.pi)
    intN_bin_t_start=np.floor(N_bin_t_start)+1
    intN_bin_t_end=np.floor(N_bin_t_end)
    intN_bin_t_start=intN_bin_t_start.astype(int)
    intN_bin_t_end=intN_bin_t_end.astype(int)
    for i in range(len(N_bin_t_start)):
        if intN_bin_t_end[i]>=intN_bin_t_start[i]:
            T_in_perbin+=int((intN_bin_t_end[i]-intN_bin_t_start[i])/m)*tbin
            #print(intN_bin_t_start[i]-1)
            T_in_perbin[np.mod(intN_bin_t_start[i],m)-1]+=(intN_bin_t_start[i]-N_bin_t_start[i])*tbin
            T_in_perbin[np.mod(intN_bin_t_end[i],m)]+=(N_bin_t_end[i]-intN_bin_t_end[i])*tbin
            rest=np.mod(intN_bin_t_end[i]-intN_bin_t_start[i],m)
            for k in range(rest):
                T_in_perbin[int(np.mod((intN_bin_t_start[i] + k), m))] += tbin
            #print(rest)
        else:
            T_in_perbin[np.mod(intN_bin_t_start[i],m)-1]+=(N_bin_t_end[i]-N_bin_t_start[i])*tbin
    return T_in_perbin
def plot_longT_V(data_file,bkg_file,epoch_file):
    epoch_info = np.loadtxt(epoch_file)
    if epoch_info.ndim == 1:
        epoch_info=np.array([epoch_info])
    t_start = epoch_info[:, 0]
    t_end = epoch_info[:, 1]
    obsID = epoch_info[:, 2]
    expT = epoch_info[:, 3]
    time = np.loadtxt(data_file)
    cts=[];bkg_cts=[]
    if not bkg_file:
        for i in range(len(obsID)):
            cts.append(len(np.where(time[:, 2] == obsID[i])[0]))
        cts = np.array(cts)
        CR = cts / expT
        CR_ERR = np.sqrt(CR * expT) / expT

    else:
        time_bkg = np.loadtxt(bkg_file)
        for i in range(len(obsID)):
            cts.append(len(np.where(time[:,2]==obsID[i])[0]))
            bkg_cts.append(len(np.where(time_bkg[:, 2] == obsID[i])[0]))
        cts=np.array(cts)
        bkg_cts=np.array(bkg_cts)
        CR=(cts-bkg_cts/12)/expT
        CR_ERR=np.sqrt(CR*expT)/expT

    # print(obsID[np.where(CR>0.0006)])
    plt.semilogy()
    plt.errorbar(t_start,CR,CR_ERR,fmt='o',capsize=3, elinewidth=1, ecolor='red')
    plt.show()

def filter_obs(src_evt,useid):
    src_evt_use = src_evt[np.where(src_evt[:-1] == useid[0])[0]]
    i=1
    while i < len(useid):
        id=useid[i]
        src_evt_use_temp=src_evt[np.where(src_evt[:-1]==id)[0]]
        src_evt_use = np.concatenate((src_evt_use, src_evt_use_temp))
        i+=1
    return src_evt_use

def phase_fold(data_file,epoch_file,p_test,bin,net_percent,shift,label,pathout):
    epoch_info=np.loadtxt(epoch_file)
    if epoch_info.ndim == 1:
        epoch_info=np.array([epoch_info])
    # print(epoch_info[4:])
    # epoch_info = np.row_stack((epoch_info[0:1],epoch_info[4:5]))
    epoch_info = epoch_info
    # useid =np.concatenate((epoch_info[:, 2][1:2],epoch_info[:, 2][5:6]))
    # tstart=np.concatenate((epoch_info[:, 0][0:1],epoch_info[:, 0][4:5]))
    # tstop= np.concatenate((epoch_info[:, 1][0:1],epoch_info[:, 1][4:5]))
    useid =epoch_info[:, 2][5:13]
    print(useid)
    tstart=epoch_info[:, 0]
    tstop =epoch_info[:, 1]
    src_evt=np.loadtxt(data_file)
    print('counts=',len(src_evt))
    src_evt=filter_obs(src_evt,useid)
    time=src_evt[:,0]
    energy=src_evt[:,1]
    # time-=time[0]

    # for k in range(len(tstart)):
    #     plt.plot([tstart[k]-time[0],tstart[k]-time[0]],[0,20],'--',color='red')
    #     plt.plot([tstop[k]-time[0], tstop[k]-time[0]], [0, 20], '--', color='green')

    plt.hist(time,bins=100,histtype='step')
    plt.show()

    def trans(t,p_test,shift):
        ti =t
        v = 1.0 /p_test
        turns = v * ti
        turns += shift
        # 初始相位
        for i in range(len(turns)):
            turns[i] = turns[i] - int(turns[i])
        return turns

    turns=trans(time,p_test,shift)
    # turns_b=trans(time_bkg,p_test,shift)
    loc=np.zeros(bin)
    loc_b= np.zeros(bin)
    for index in turns:
        loc[int(index*bin)] += 1
    # for index in turns_b:
    #     loc_b[int(index * bin)] += 1

    x = np.array([(i / bin + 0.5 / bin) for i in range(bin)])
    AM=1-min(loc)/max(loc)
    A0=AM/(2-AM)
    print('A0={0}'.format(A0))

    src_bkg = 1 - net_percent
    bkg_y = len(time) * src_bkg
    # bkg_y_low=bkg_y-bkg_y**0.5
    # bkg_y_high = bkg_y+bkg_y**0.5

    bkg_y /= bin
    b_1sigma = poisson_conf_interval(bkg_y, interval='frequentist-confidence').T
    bkg_y_low=b_1sigma[0]
    bkg_y_high=b_1sigma[1]

    fig=plt.figure(1,(10,7.5))
    ax1 = fig.add_subplot(111)

    #plt.plot([0, 2], [bkg_y, bkg_y], '--')
    bkg_x = [0, 2]
    plt.fill_between(bkg_x, bkg_y_low, bkg_y_high,facecolor = 'blue', alpha = 0.5)
    # y3=np.concatenate((loc_b,loc_b))
    # y3*=2.1590383026593E-05/5.0169307007741E-05


    x2=np.concatenate((x,x+1))
    y2=np.concatenate((loc,loc))
    T_in_perbin = get_T_in_mbins(epoch_info, 2 * np.pi / p_test, bin, shift * 2 * np.pi)

    correct_gap = T_in_perbin / (sum(T_in_perbin) / len(T_in_perbin))
    print(correct_gap)
    y2 /= np.concatenate((correct_gap, correct_gap))
    y2_err=np.array(poisson_conf_interval(y2,interval='frequentist-confidence'))
    y2_err[0]=y2-y2_err[0]
    y2_err[1]=y2_err[1]-y2

    #label='F1'
    plt.title("#{0} P={1:.2f},C={2}".format(label,p_test,str(len(time))), fontsize = 18)
    plt.xlabel('phase',font1)
    plt.ylabel('counts/bin',font1)
    plt.tick_params(labelsize = 18)
    plt.ylim(0,(np.max(y2)+np.max(y2)**0.5)*1.05)
    plt.step(np.concatenate(([0],x2)),np.concatenate(([y2[0]],y2)),color='red')
    #plt.step(x2,y3,'--',color='green')
    print(np.size(y2))
    print(np.size(y2_err))
    #plt.errorbar(x2 - 0.5 / bin, y3, yerr = y3 ** 0.5, fmt = '.', capsize = 1, elinewidth = 0.5, ecolor = 'green')
    plt.errorbar(x2 - 0.5 / bin, y2, yerr = y2_err, fmt = '.', capsize = 1, elinewidth = 1, ecolor = 'red')

    ax2 = ax1.twinx()
    x22=x2
    y22=y2/np.mean(y2)
    yhigh=(np.max(y2)+np.max(y2)**0.5)*1.05/np.mean(y2)
    ax2.set_ylabel('Normalized flux',font1)
    ax2.plot([0,2],[1.0,1.0],'--',color='green')
    ax2.set_ylim([0,yhigh])
    ax2.tick_params(labelsize=18)

    plt.savefig(pathout+'pfold_lc_{0}002.eps'.format(label))
    plt.show()
    #plt.close()
# phase_fold('513_bkg.txt','LW_epoch.txt',5334.75593,bin = 30, net_percent = 0.9, shift = 0.2, label = 12)

path_NSC='/Users/baotong/Desktop/period/txt_all_obs_IG/'
#path_NSC='/Users/baotong/Desktop/period/txt_all_obs_G/'
#path_NSC='/Users/baotong/Desktop/period/txt_all_obs_I/'
path_LW='/Users/baotong/Desktop/period_LW/txt_all_obs/'
path_ND='/Users/baotong/Desktop/period_gc/txt_all_obs/'
path_Tuc='/Users/baotong/Desktop/period_Tuc/txt_all_obs_0.5_8/'
#path_Tuc='/Users/baotong/Desktop/period_Tuc/txt_all_obs/'
path_omg='/Users/baotong/Desktop/period_omg/txt_all_obs_0.5_8/'
path_terzan='/Users/baotong/Desktop/period_terzan5/txt_all_obs_0.5_8/'
path_M28='/Users/baotong/Desktop/period_M28/txt_all_obs_0.5_8/'
path_NGC6397='/Users/baotong/Desktop/period_NGC6397/txt_all_obs_0.5_8/'
path_NGC6752='/Users/baotong/Desktop/period_NGC6752/txt_all_obs_0.5_8/'
path_CDFS='/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep2/'
path_xmmCDFS='/Users/baotong/Desktop/CDFS/xmm_txt/'
obsid=700014
path_eSASS='/Users/baotong/eSASS/data/raw_data/47_Tuc/txt/txt_reg2_psf75_0.2_5/'
# path_eSASS='/Users/baotong/eSASS/data/raw_data/47_Tuc/txt/txt_psf75_{0}/'.format(obsid)
figurepath='/Users/baotong/Desktop/aas/pCV_GC/figure/47Tuc/'
if __name__=='__main__':
    path=path_Tuc
    period=42955.32646/2
    # dataname='481_{0}.txt'.format(obsid)
    dataname = '350.txt'
    net_p=0.8
    epoch_file = path + 'epoch_src_' + dataname
    plot_longT_V(data_file=path + dataname, bkg_file=None,epoch_file=epoch_file,)
    # epoch_file=path+'epoch_47Tuc_{0}.txt'.format(obsid)
    phase_fold(path + dataname, epoch_file, period, bin = 20, net_percent = net_p, shift = 0.0, label =dataname[0:-4],pathout=figurepath)
