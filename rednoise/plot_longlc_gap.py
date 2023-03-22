#!/bin/bash
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 19 18:13:40 2022
@author: baotong
modified version of plot lc with gap
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
import hawkeye as hawk
import rednoise

def judge_GTI(t,GTI):
    outjudge=False
    for i in range(len(GTI)):
        if GTI[i][0]<= t<= GTI[i][1]:
            outjudge= True
            break
    return [outjudge,i]
# def plot_lcwithgap_period(src_evt,epoch_info,period,len_bin=5000,ifsin=None):
#     time=src_evt[:,0];energy=src_evt[:,1]
#     obs_id = epoch_info[:, 2]
#     evt_id = src_evt[:, 2]
#     GTI=epoch_info[:,0:2]
#
#     x_all=[];y_all=[];y_err=[];z_all=[];minigap=[];gap_intP=0
#     timepoint=time[0];i=0  ## start the loop
#
#     while timepoint+(i+1)*len_bin< time[-1]:
#         t_left=timepoint+i*len_bin;t_right=timepoint+(i+1)*len_bin
#         jg_templ=judge_GTI(t_left,GTI);jg_tempr=judge_GTI(t_right,GTI)
#         if (jg_templ[0])&(jg_tempr[0]):
#             if (jg_templ[1])!=(jg_tempr[1]):
#                 print('caution, one time bin cross two obs!')
#             x_all.append(t_left / 2 + t_right / 2-period*gap_intP)
#             counts = np.where((time > t_left) & (time <= t_right))[0]
#             y_all.append(counts / len_bin)
#             z_all.append(1)
#         if (jg_templ[0])&(not jg_tempr[0]):
#             obs_point=jg_templ[1]
#             tstop=epoch_info[obs_point][1]
#             x_all.append(t_left/2+tstop/2-period*gap_intP)
#             counts = np.where((time > t_left) & (time <= tstop))[0]
#             y_all.append(counts / (tstop-t_left))
#             phase = np.mod(epoch_info[i][1], period)
#             minigap.append([epoch_info[i][1], epoch_info[i][1] + (1 - phase) * period])
#
#             z_all.append(0.5)
#         if (not jg_templ[0])&(jg_tempr[0]):
#             obs_point=jg_templ[1]
#             tstart=epoch_info[obs_point][0]
#             x_all.append(tstart/2+t_right/2-period*gap_intP)
#             counts = np.where((time > tstart) & (time <= t_right))[0]
#             y_all.append(counts / (t_right-tstart))
#             z_all.append(0.5)
#         if (not jg_templ[0])&(not jg_tempr[0]):
#             x_all.append(t_left/2+t_right/2-period*gap_intP)
#             y_all.append(0)
#             y_err.append(0)
#             z_all.append(0)
#         i+=1
#     all_data=np.column_stack((x_all,y_all,y_err,z_all))
#     for k in range(len(obs_id)):
#         phase=np.mod(epoch_info[i][1],period)
#         minigap.append([epoch_info[i][1],epoch_info[i][1]+(1-phase)*period])
#
#     return None

def poisson_conf(loop,N,threshold):
    num=0
    for i in range(loop):
        a = np.random.random(N)
        if np.max(a) - np.min(a) < threshold:
            num += 1
    return num/loop

def plot_lc(source_id,period=None):
    # path='/Users/baotong/Desktop/period_NGC6397/txt_all_obs_0.5_8/'
    path='/Users/baotong/Desktop/period_Tuc/txt_startover/txt_all_obs_p{0}/'.format(90)

    (evt_file,epoch_file)=rednoise.load_data(source_id,ecf=90,ifobsid=[2735,2736],ifexpT=None)
    # expvalue=[2.63e7,2.63e7,2.63e7,2.63e7,1.61e7,1.97e7,9.7e6]

    obs_id=epoch_file[:,2]
    evt_id=evt_file[:,2]
    evt_list=evt_file[:,0]
    period =period
    # bin_len=period/8
    bin_len=2000.
    x_all = [];
    y_all = [];
    xerr_all = [];
    yerr_all = []
    for i in range(len(obs_id)):
        x=[];y=[];xerr=[];yerr=[]
        time=evt_list[np.where(evt_id==obs_id[i])]
        if len(time)<2:
            print('attention')
            continue

        k=0;
        while (time[0]+k*bin_len)<time[-1]:
            evt_temp=time[np.where((time>(time[0]+k*bin_len))&(time<(time[0]+(k+1)*bin_len)))]
            x.append(time[0]+(k+0.5)*bin_len)
            if (time[0]+(k+1)*bin_len)<time[-1]:
                exp_temp=bin_len
            else:
                exp_temp=time[-1]-(time[0]+k*bin_len)
            y.append(len(evt_temp)/exp_temp)
            xerr.append(0.5*bin_len)
            yerr.append(np.sqrt(len(evt_temp))/exp_temp)
            k+=1
        # print(len(x))

        x_all.append(x);y_all.append(y);xerr_all.append(xerr);yerr_all.append(yerr)

    ## plot the max point phase ##

    plt.figure(2)
    for i in range(0, len(x_all)):
        maxf_time=x_all[i][np.argmax(y_all[i])]

        phase=np.mod(maxf_time,period)/period
        plt.scatter(i,phase)
    phase0=np.mod(x_all[0][np.argmax(y_all[0])],period)/period
    plt.plot([0,len(x_all)],[phase0-0.5*bin_len/period,phase0-0.5*bin_len/period],'--',color='grey')
    plt.plot([0,len(x_all)],[phase0+0.5*bin_len/period,phase0+0.5*bin_len/period],'--',color='grey')
    plt.ylim(0,1.01)
    plt.show()
    return_time=x_all[0];return_flux=y_all[0]
    for i in range(1,len(x_all)):
        return_time=np.concatenate((return_time,x_all[i]))
        return_flux=np.concatenate((return_flux,y_all[i]))
    ## 把数组拉平，称为return_time
    ##=== MJD=== ##
    plt.figure(1,figsize=(12,4))
    plt.errorbar(x_all[0], np.array(y_all[0])*25.423, xerr=xerr_all[0],
                 yerr=np.array(yerr_all[0])*25.423,fmt='.', capsize=1)

    ### sine function ####
    def fmax(x, a, b):
        return a * np.sin(x * 2 * np.pi / period-0.0*2*np.pi) + b

    xsin=np.linspace(x_all[0][0],x_all[-1][-1],10000000)
    ysin=0.001*np.sin(2*np.pi/period*xsin)+0.01
    # plt.plot(xsin,ysin)
    for i in range(1,len(x_all)):
        if (np.mod(x_all[i][0],period)/period)>(np.mod(x_all[i-1][-1],period)/period):
            gap=(np.mod(x_all[i][0],period)/period-np.mod(x_all[i-1][-1],period)/period)*period
        else:
            gap = (np.mod(x_all[i][0], period) / period - np.mod(x_all[i - 1][-1], period) / period) * period+period
        # print(gap)
        x_all[i] =x_all[i]- x_all[i][0] + x_all[i - 1][-1] + gap
        plt.errorbar(x_all[i],np.array(y_all[i])*25.423,xerr=xerr_all[i],
                     yerr=np.array(yerr_all[i])*25.423,fmt='.', capsize=1)
        plt.fill_between([x_all[i-1][-1],x_all[i][0]],np.max(return_flux*25.423),facecolor='yellow',alpha=0.2)

    x_all_list=x_all[0];y_all_list=y_all[0]
    for i in range(1,len(x_all)):
        x_all_list=np.concatenate((x_all_list,x_all[i]))
        y_all_list= np.concatenate((y_all_list, y_all[i]))

    fita, fitb = optimize.curve_fit(fmax, x_all_list, y_all_list, [0.005, 0.01])
    xperiod=x_all[0][0]+period*np.arange(0,int((x_all_list[-1]-x_all_list[0])/period),1)
    # for i in range(len(xperiod)):
    #     plt.plot([xperiod[i],xperiod[i]],[0,np.max(return_flux)],'--',c='grey')
    # plt.plot(xsin[np.where(xsin<x_all[-1][-1])], fmax(xsin[np.where(xsin<x_all[-1][-1])], fita[0], fita[1]),'--',c='r')
    print(xsin)
    plt.plot(xsin[np.where(xsin<x_all[-1][-1])], fmax(xsin[np.where(xsin<x_all[-1][-1])] ,0.003,0.003)*25.423,'-',c='grey')
    plt.ylim(0,2)
    plt.xlabel('Time',hawk.font1)
    plt.ylabel(r'Photon flux ($\rm 10^{-4}~ph~s^{-1}~cm^{-2}$)',hawk.font1)
    plt.tick_params(labelsize=16)
    figurepath='/Users/baotong/Desktop/aas/pXS_Tuc_mod1/figure/'
    # plt.savefig(figurepath+f'{source_id}_long_gap_org.pdf',bbox_inches='tight', pad_inches=0.1)
    plt.show()
    return [return_time,return_flux]



def plot_phase_obs(source_id,period):
    path='/Users/baotong/Desktop/period_Tuc/txt_startover/txt_all_obs_p{0}/'.format(90)
    (evt_file,epoch_file)=rednoise.load_data(source_id,ecf=90,ifobsid=None,ifexpT=10000)

    obs_id=epoch_file[:,2]
    evt_id=evt_file[:,2]
    evt_list=evt_file[:,0]
    period =period
    binnum=10
    bin_len=period/(binnum+0.0)
    # bin_len=1000.
    x_all = [];
    y_all = [];
    xerr_all = [];
    yerr_all = []
    for i in range(len(obs_id)):
        x = [];
        y = [];
        xerr = [];
        yerr = []
        time = evt_list[np.where(evt_id == obs_id[i])]
        if len(time) < 2:
            print('attention')
            continue

        k = 0;
        while (time[0] + k * bin_len) < time[-1]:
            evt_temp = time[np.where((time > (time[0] + k * bin_len)) & (time < (time[0] + (k + 1) * bin_len)))]
            x.append(time[0] + (k + 0.5) * bin_len)
            if (time[0] + (k + 1) * bin_len) < time[-1]:
                exp_temp = bin_len
            else:
                exp_temp = time[-1] - (time[0] + k * bin_len)
            y.append(len(evt_temp) / exp_temp)
            xerr.append(0.5 * bin_len)
            yerr.append(np.sqrt(len(evt_temp)) / exp_temp)
            k += 1
        # print(len(x))

        x_all.append(x);
        y_all.append(y);
        xerr_all.append(xerr);
        yerr_all.append(yerr)
    (fig,axes)=plt.subplots(nrows=6,ncols=1,sharex=True)

    for i in range(len(x_all)):
        max_index=np.argmax(y_all[i])
        ##选一段时间来plot
        shift_index=0
        shift1=binnum/2;shift2=binnum/2
        if np.mod(binnum,2)>0:
            shift1=int(binnum/2);shift2=int(binnum/2)+1
        while shift_index <=shift2:
            if (max_index-shift1+shift_index)>=0 and (max_index+shift2+shift_index)<=len(x_all[i]):
                break
            else:shift_index+=1.
        if shift_index>=shift2:
            shift_index =0
            while shift_index >=-shift1:
                if (max_index - shift1 + shift_index) >= 0 and (max_index + shift2 + shift_index) <= len(
                        x_all[i]):
                    break
                else:
                    shift_index -= 1.

        phase_obs=np.mod(np.array(x_all[i]),period)/period
        i1=int(max_index-binnum/2+shift_index);i2=int(max_index+binnum/2+shift_index)
        # if i1<=0 or i2>=len(x_all[i]):
        #     i1=0;i2=i1+int(period/bin_len)
        axes[i].errorbar(phase_obs[i1:i2],
                     y_all[i][i1:i2],
                     yerr=yerr_all[i][i1:i2],
                     xerr=np.zeros(i2-i1)+bin_len/period/2,fmt='.')
        # axes[i].set_yscale('log')
    plt.show()
if __name__=="__main__":
    plot_lc('185',period=8517.19)
    # plot_phase_obs('290',period=46082.9493)
    # print(poisson_conf(10000,6,0.15))