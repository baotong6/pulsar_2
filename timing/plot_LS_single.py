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
from scipy.fftpack import fft,ifft
import scipy.signal as ss
import scipy.stats as stats
import random
from astropy.timeseries import LombScargle
#path='/Volumes/halo/chandra_obs_hr_m04_for_baotong/'
# path='/Volumes/pulsar/xmm_obs_16/'
# path='/Volumes/pulsar/xmm_obs_16/0112870101/'
# path='/Users/baotong/Desktop/period_LW/xmm_CV/'
#dataname='2560.txt'
#dataname="Muno2.txt"
def get_fig_LS(dataname,exptime):
    def get_LS(freq,len_bin,dataname):
        time=np.loadtxt(path+'txt/'+dataname+'.txt')[:,0]

        def get_hist(t,len_bin):
            t_test=t
            t_test-=t_test[0]
            a=[0 for i in range(int((t_test[-1]-t_test[0])/len_bin)+1)]
            for i in range(len(t_test)):
                # print int(t_test[i]/len_bin)
                # print int(obs_d[-1][1]/len_bin)
                a[int(t_test[i]/len_bin)]+=1
            return a

        cts = len(time)
        y=get_hist(time,len_bin)
        y=np.array(y)
        x=np.arange(time[0],time[-1],len_bin)
        ### plot candiate period in lc ###
        # a=[11758.1*i for i in range(4)]
        # plt.step(x, y)
        # for i in range(len(a)):
        #     plt.plot([a[i],a[i]],[0,60],'--',color='green')
        # plt.show()

        ### plot candiate period in lc ###

        # plt.scatter(x,window)
        # plt.show()

        LS = LombScargle(x, y, dy=1, normalization='psd', fit_mean=True, center_data=True).power(freq,
                                                                                                 method='cython')

        return [LS,cts]

    freq_1=[1/7000.,1e-9,50000]
    #freq_1=[1e-4,1e-6,10000]
    freq=freq_1[0]+freq_1[1]*np.arange(freq_1[2])

    ### uneven sampled freq ###
    p_unsamp=[]
    freq_unsamp=[]
    p_unsamp.append(5e-4*exptime)
    while p_unsamp[-1]<exptime:
        if p_unsamp[-1]<0.3*exptime:
            p_unsamp.append(p_unsamp[-1]+p_unsamp[-1]**2/(2*exptime*2))
        elif 0.3*exptime<p_unsamp[-1]<0.5*exptime:
            p_unsamp.append(p_unsamp[-1] + p_unsamp[-1] ** 2 / (3 * exptime * 2))
        else:
            p_unsamp.append(p_unsamp[-1] + p_unsamp[-1] ** 2 / (4* exptime * 2))
    p_unsamp=np.array(p_unsamp)
    freq_unsamp=1.0/p_unsamp

    ### uneven sampled freq ###

    def get_conf90(Np,Ns):
        return(-np.log(1-0.99**(1.0/(Ns*Np))))

    print('running')
    #res_normal=get_LS(freq_unsamp,len_bin=20.,dataname=dataname)
    res_normal = get_LS(freq, len_bin = 20., dataname = dataname)
    LS=res_normal[0]
    cts=res_normal[1]

    conf90=get_conf90(len(freq),16)

    plt.figure(1,(6,6))
    #plt.subplot(211)
    plt.title(dataname+',cts={0}'.format(cts))
    plt.semilogx()
    plt.step(freq,LS)
    plt.plot([freq[0],freq[-1]],[conf90,conf90],'--')
    # plt.subplot(212)
    # plt.title('window_function')
    # plt.semilogx()
    #plt.step(freq, LS_W)

    #plt.savefig(path+'fig/un_uniform_pn_bin5/'+dataname+'.eps')

    plt.show()
    #plt.close()


# evt_list_ID=np.arange(1,49)
# evt_list_ID=np.delete(evt_list_ID,[20,31,36,37,38,39,40,41,45])
# evt_list_name=[]
# for i in range(len(evt_list_ID)):
#     if evt_list_ID[i]<10:
#         evt_list_name.append('src'+'0'+str(evt_list_ID[i]))
#     else:
#         evt_list_name.append('src'+str(evt_list_ID[i]))
#,1906,3236,1671,2148,2176,1612,3107,3401,2276]
         # 2338,3051,2942,2652,2187,2962,2798,2560,1537,1502,3594,1937,
         # 3021,3564,2888,3283,2841,2178,3006,3155,2226,3231,3412,
         # 566,2287,2380,2198,914,3178]

# evt_list_name=['src1pn','src2pn','src3pn','src4pn','src5pn',
#                'src6pn','src7pn','src8pn','src9pn','src10pn',
#                'src11pn','src12pn','src14pn','src15pn','src15pn','src16pn']
# exptime=[22713,58317,63086,120500,114000,23000,50755,29617,58319,49200,
#          31569,82900,39620,58211,23861,98920]

path='/Users/baotong/Desktop/period_LW/xmm_CV/0099020301/'
evt_list_name=['OY_CAR_pn']
exptime=[49582]
for i in range(len(evt_list_name)):
    get_fig_LS(evt_list_name[i],exptime[i])
# dataname='src22mos1'
