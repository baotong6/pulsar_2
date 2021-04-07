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
import datetime
from scipy.interpolate import lagrange
from scipy import interpolate
from scipy.fftpack import fft,ifft
import scipy.signal as ss
import scipy.stats as stats
import random
import scipy.stats as stats
import matplotlib.pyplot as plt
import matplotlib.style as style
from IPython.core.display import HTML
import read_data as data
#import read_data_lw as data
#
#import fits_to_txt
starttime=datetime.datetime.now()
def sim_time(dataname):
    time = data.get_data(dataname)[0]
    energy = data.get_data(dataname)[1]
    obs_ID = data.get_data(dataname)[2]
    t = data.get_data(dataname)[3]
    E = data.get_data(dataname)[4]
    dict = data.get_data(dataname)[5]

    t = np.array(t)

    sim_t = t - t
    for i in range(len(t)):
        sim_t[i] = np.random.random(len(t[i])) * (t[i][-1] - t[i][0]) + t[i][0]
    return sim_t

def get_Z2(dataname,freq):
    time=data.get_data(dataname)[0]
    energy=data.get_data(dataname)[1]
    obs_ID=data.get_data(dataname)[2]
    t=data.get_data(dataname)[3]
    E=data.get_data(dataname)[4]
    dict=data.get_data(dataname)[5]

    sss=0
    for index in t:
        for i in range(len(index)):
            sss+=1
    if sss!=len(time):
        print('error')
        print('sss=%d %sss')
        print('len(time)=%d len(time)')
    #time=t[0]
    # #
    delete_i=[]
    use_ID_o=[]
    #dict[242,2951][0]

    for i in dict:
        if dict[i][-2]<5000:
            delete_i.append(i)
        j=0
        while j <len(dict[i][1]):
        ##filter energy
            if dict[i][1][j]<2000:
                dict[i][1]=np.delete(dict[i][1],j)
                dict[i][0]=np.delete(dict[i][0],j)

            else:j+=1

    for item in delete_i:
        dict.pop(item)
    for i in dict:
        use_ID_o.append(i)

    epoch0 = np.array([2943,3663])
    epoch0 = np.array([3392, 3393])
    epoch1=np.array([2943,3663,3392,3393,3665])
    epoch2=np.array([5950,5951,5952,5953,5954])
    epoch3=np.array([9169,9170,9171,9172])
    epoch4=np.array([10556,11843])
    epoch_g1=np.array([13847,14427,13848,13849,13846,14438,13845])
    #above is epoch in ACIS-I

    # epoch1=np.array([13847,14427,13848,13849,13846,14438,13845])
    # epoch2=np.array([14461,13853,13841,14465,14466,13842,13839,13840,14432,13838,13852,14439])
    # epoch3=np.array([14462,14463,13851,15568,13843,15570,14468])
    #epoch_gc = np.array([14897, 17236, 17239, 17237, 18852, 17240, 17238, 20118, 17241, 20807, 20808])
    #use_ID_o=np.concatenate((epoch2,epoch3,epoch4))
    #use_ID_o=epoch_gc
    #use_ID_o=epoch_g1
    #use_ID_o=np.array([3549])
    #use_ID_o = epoch1
    #use_ID_o=epoch3
    #print(dict[use_ID[0]][0])
    #print(dict)
    #use_ID=use_ID_o
    not_in_id=[]
    for i in range(len(use_ID_o)):
        if use_ID_o[i] in dict:
            continue
        else:
            not_in_id.append(i)
    use_ID=np.delete(use_ID_o,not_in_id)

    if len(use_ID)<1:
        return [0,0]
    elif len(use_ID)==1:
        time=dict[use_ID[0]][0]
    else:
        time=np.concatenate((dict[use_ID[0]][0],dict[use_ID[1]][0]))

        for i in range(2,len(use_ID)):
            time=np.concatenate((time,dict[use_ID[i]][0]))
    plt.hist(time-time[0],bins=50,histtype='step')
    plt.show()
    N = len(time)

    def turns(t,f):
        ti=t-t[0]
        #print ti
        v=f
        #p_test = 1.0/5.500550055005501e-06
        p=1.0/v
        # pdot=-vdot/(v*v)
        # vddot = 2.0 * pdot * pdot / (p * p * p)
        freq = v # + ti * vdot + ti * ti / 2.0 * vddot
        preq = 1.0 / freq
        turns=v*ti
        INT_turns=np.trunc(turns)
        turns=turns-INT_turns
        turns = 2.0*np.pi*turns #+ vdot * ti * ti / 2.0 + vddot * ti * ti * ti / 6.0
        #print turns
        #turns=list(turns)
        #turns=np.array(turns[0])
        #初始相位

        return turns
    #print time
    Z2=[]

    for fi in freq:
        Z2.append((2.0 / N)*(sum(np.cos(turns(time,fi)))**2+sum(np.sin(turns(time,fi))**2)))
    cts=len(time)
    return [Z2,cts,time[-1]-time[0]]

def get_Z2_sim(dataname):
    sim_t=sim_time(dataname)
    simtime=data.tran_t_time(sim_t)
    Zr_temp = []
# simtime=np.reshape(sim_time(dataname),(1,len(time)))
    for fi in freq:
        Zr_temp.append((2.0 / N)*(sum(np.cos(turns(simtime,fi)))**2+sum(np.sin(turns(simtime,fi))**2)))
    Zr_temp=np.array(Zr_temp)
    return Zr_temp

# Zr=get_Z2_sim(dataname)
# for i in range(10):
#     Zr_temp=get_Z2_sim(dataname)
#     Zr=(Zr+Zr_temp)/2.0

def make_period_range(pmin,pmax,expT):
    P=[pmin]
    while P[-1]<pmax:
        dP=0.1*P[-1]**2/(expT-P[-1])
        P.append(P[-1]+dP)
    return np.array(P)
#print(len(make_period_range(300,3600,1618692)))
def get_fig_Z2(dataname):
    border=10000
    vary=np.array([i for i in range(0,border)])
    freq=1/4.+vary*1.e-5
    #freq=1/100000.+vary*1.e-6
    #P = make_period_range(3600, 8000, 2e8)
    #freq=1.0/P
    #print(len(freq))

    # border=10000
    # vary=np.array([i for i in range(0,border)])
    # freq=1.e-4+vary*1.e-6
    # p99_Z2 = 23.01325331
    # p90_Z2 = 18.31207802

    p99_Z2=27.
    p90_Z2=22.

    (Z2,cts,expT)=get_Z2(dataname,freq)
    if Z2==0 or cts==0:
        return 'none'
    Z2 = np.array(Z2)
    #(Z2_standard,cts_standard)=get_Z2('Cannonball.txt',freq)
    #Z2_standard=np.array(Z2_standard)

    #plt.subplot(211)
    plt.title(dataname[0:-4]+',cts={0}'.format(cts))
    plt.semilogx(freq,[p99_Z2 for i in range(len(Z2))],'--',color='black')
    plt.semilogx(freq,[p90_Z2 for i in range(len(Z2))],'--',color='red')
    plt.step(freq,Z2,color='black')
    #plt.step(freq,Z2_standard*max(Z2)/max(Z2_standard),color='pink')
    plt.text(0.5,p99_Z2+1,"99%")
    plt.text(0.5,p90_Z2+1,"90%")

    # plt.subplot(212)
    # plt.semilogx(freq,[p99_Z2 for i in range(len(Z2))],'--',color='black')
    # plt.semilogx(freq,[p90_Z2 for i in range(len(Z2))],'--',color='red')
    # plt.step(freq,Z2-Z2_standard*max(Z2)/max(Z2_standard))

    print('running')
    #plt.savefig('/Users/baotong/Desktop/li_pr/result_final/fig_NSC_I/{0}_Z2.eps'.format(dataname[0:-4]))
    #plt.savefig(path + 'fig_Z_high_res/' + dataname[0:-4] + "_Z2.jpg")
    plt.show()
    #plt.savefig(path+'fig_Z_high_res/'+dataname[0:-4]+"_Z2.jpg")
    #plt.close()

    # 卡方分布——画图

    # PLOTTING CONFIG 绘图配置

    # PDF  概率密度函数

    # # CDF 累积概率密度函数
    # plt.plot(np.linspace(0, 20, 100), stats.chi2.cdf(np.linspace(0, 20, 100), df=2))  # 绘制累积概率密度函数
    #
    # weights1 = np.ones_like(Z2) / float(len(Z2))
    # plt.subplot(211)
    # plt.hist(Z2,bins=200,weights=weights1, histtype='step')
    # plt.plot(np.linspace(0, 40, 1000), stats.chi2.pdf(np.linspace(0, 40, 1000), df=2))  # 绘制0到20的卡方分布曲线,给定自由度为2
    # plt.fill_between(np.linspace(0, 40, 1000), stats.chi2.pdf(np.linspace(0, 40, 1000), df=2), alpha=0.15)  # 填充曲线
    # plt.subplot(212)
    # plt.hist(Z2_standard,bins=200,weights=weights1,histtype='step')
    # plt.plot(np.linspace(0, 40, 1000), stats.chi2.pdf(np.linspace(0, 40, 1000), df=2))  # 绘制0到20的卡方分布曲线,给定自由度为2
    # plt.fill_between(np.linspace(0, 40, 1000), stats.chi2.pdf(np.linspace(0, 40, 1000), df=2), alpha=0.15)  # 填充曲线
    # #plt.step(a[1],a[0])
    # plt.show()
    # print(dataname)

    #plt.savefig(path+'fig_Z/'+dataname[0:-4]+"_Z2.eps")

item='1502'
get_fig_Z2(str(item) + '.txt')

# endtime=datetime.datetime.now()
# print((endtime-starttime).seconds)

# cand_ID=fits_to_txt.cand_ID
# for item in cand_ID[1:12]:
#     get_fig_Z2(str(int(item)) + '.txt')


# src_in_G=np.loadtxt(path+'txt_G/'+'G_src.txt')
# src_in_G=src_in_G.astype(int)
# for item in src_in_G:
#     get_fig_Z2(str(int(item)) + '.txt')

# cand_ID=np.arange(1,518)
# for item in cand_ID:
#     get_fig_Z2(str(int(item)) + '.txt')