# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import string
import plot_result_eclipse as sim_e
import plot_result_sin as sim_sin
import read_csv as data
path_txt='/Users/baotong/Desktop/period_LW/txt_all_obs/'
ID_num=np.linspace(1,847,847)#源的数目
# ID_num=np.linspace(1,3619,3619)
# [87.48, 180.45, 274.98, 369.05, 461.18, 555.41, 649.47, 744.78, 836.88]
cts = []
for id in ID_num:
    file_name = str(int(id)) + '.txt'
    file_info = np.loadtxt(path_txt + file_name)
    cts.append(len(file_info[:, 0]))
Lum_LW=np.array([31.26,32.36,31.82,32.35,31.75,31.97
                        ,31.83,31.25,31.54,31.33,31.98
                        ,31.88,31.93,32.11,31.93
                        ,32.05,31.62,32.92,31.73,32.78,31.97,30.88])
cts_LW=np.array([202,902,394,784,293,437,335,121,347,211,438,760,512,823,
                     487,535,214,3402,263,1963,307,138])
cts=np.array(cts)
# plt.hist(cts,bins=100)
# cts=cts[np.where(cts>100)]
# # print(len(np.where(cts>100)[0]))
# #plt.scatter(cts_LW,Lum_LW)
#
# plt.show()


def get_DR_DNe():
    ##for DNe
    # cts_bin = np.linspace(60, 900, 10)#按光子数取bin
    # cts_bin=np.concatenate((cts_bin,[1000,1500,1800]))

    cts_bin=np.linspace(75,3575,36)
    cts_bin=np.concatenate(([35],cts_bin))

    cts_hist=plt.hist(cts,histtype = 'step',bins=cts_bin)[0]
    cts_hist*=0.9
    #plt.savefig('/Users/baotong/Desktop/aas/V63/figure/cts_hist.eps')
    plt.show()
    detect_rate_p1=sim_e.get_info(sim_e.path_1,5258.)[0][0]
    detect_rate_p2=sim_e.get_info(sim_e.path_2,15258.)[0][0]
    detect_rate_p1.extend(np.zeros(24) + 1.)
    detect_rate_p2.extend(np.zeros(24) + 1.)
    print(cts_hist)
    a=(cts_hist[1:6]*detect_rate_p1[1:6]*0.18*0.83+cts_hist[1:6]*detect_rate_p2[1:6]*0.2*0.17)*0.53
    de_eclipse=(np.sum(cts_hist[1:6]*detect_rate_p1[1:6]*0.18*0.92)+np.sum(cts_hist[1:6]*detect_rate_p2[1:6]*0.2*0.08))*0.77
    # 0.75 and 0.25分别代表低于period gap和高于period gap的DN比例
    print(de_eclipse)
    return de_eclipse
#get_DR_DNe()

def get_DR_polar():
    # for polar
    DR_p1=[]
    DR_p2=[]
    cts_bin=np.linspace(75,3575,36)
    cts_bin=np.concatenate(([35],cts_bin))
    # cts_bin = np.concatenate((cts_bin,[1000,1500,1800]))
    path_p1 = '/Users/baotong/Desktop/period_LW/simulation/simulation_LW_5540/'
    path_p2 = '/Users/baotong/Desktop/period_LW/simulation/simulation_LW_45540/'
    detect_rate_p1 = sim_sin.get_info(path_p1, 5540.)[0]
    detect_rate_p2 = sim_sin.get_info(path_p2, 45540.)[0]
    for i in range(len(detect_rate_p1)):
        DR_p1.append(np.mean(detect_rate_p1[i]))
        DR_p2.append(np.mean(detect_rate_p2[i]))
    DR_p1.extend(np.zeros(30)+1.)
    DR_p2.extend(np.zeros(30)+1.)
    print(DR_p1)
    print(DR_p2)
    cts_hist = plt.hist(cts, histtype = 'step', bins = cts_bin)[0]*0.9
    plt.close()
    alpha=0.18#for fraction
    beta = 0.25  # ??? for polar intrinsic

    det_n1=cts_hist*DR_p1*0.82*alpha*beta;det_n2=cts_hist*DR_p1*0.18*alpha*beta
    de_sin_polar=np.sum(det_n1[1:])+np.sum(det_n2[1:])
    Lum_LW=np.array([31.26,32.36,31.82,32.35,31.75,31.97
                        ,31.83,31.25,31.54,31.33,31.98
                        ,31.88,31.93,32.11,31.93
                        ,32.05,31.62,32.92,31.73,32.78,31.97,31.94,30.88])
    cts_LW=np.array([202,902,394,784,293,437,335,121,347,211,438,760,512,823,
                     487,535,214,3402,263,1963,307,1039,138])
    cts_LW_ex_IP=np.array([902,394,784,293,437,335,121,347,211,438,760,512,823,
                     487,535,214,3402,263,307])
    plt.scatter(cts_LW,Lum_LW)
    # plt.hist(det_n1+det_n2,bins=20,histtype = 'step')
    DET=det_n1+det_n2
    DET=np.concatenate((DET,[0.]))

    a=plt.hist(cts_LW_ex_IP,bins=cts_bin,histtype='step')
    plt.close()
    plt.figure(1,(8,6))
    font1 = {'family': 'Normal',
             'weight': 'normal',
             'size': 15, }

    plt.ylim(0.,4.5)
    plt.tick_params(labelsize=15)
    plt.step(cts_bin, DET, where = 'post')
    plt.step(cts_bin,np.concatenate((a[0],[0.])),where = 'post')


    plt.legend(['predicted','detected'])
    plt.xlabel('Counts',font1)
    plt.ylabel('Number of sources',font1)
    plt.savefig('/Users/baotong/Desktop/aas/mod_MN_pCV/figure/sim_LW/est_obs.eps')
    plt.show()
    print(de_sin_polar)
    return de_sin_polar
get_DR_polar()
def get_DR_IP():
    # for polar
    DR_p1 = []
    cts_bin = np.linspace(75, 3575, 36)
    cts_bin = np.concatenate(([35], cts_bin))
    path_p1 = '/Users/baotong/Desktop/period_LW/simulation/simulation_LW_554/'
    detect_rate_p1 = sim_sin.get_info(path_p1, 554.)[0]
    print(detect_rate_p1)
    for i in range(len(detect_rate_p1)):
        #DR_p1.append(np.mean(detect_rate_p1[i]))

        DR_p1.append(detect_rate_p1[i][0])
    DR_p1.extend(np.zeros(30) + 1.)
    print(DR_p1)
    cts_hist = plt.hist(cts, histtype = 'step', bins = cts_bin)[0]*0.9
    plt.close()
    #print(cts_hist)
    alpha = 0.03  # for fraction
    beta = 0.5  # ??? for IP intrinsic
    det_n1 = cts_hist * DR_p1 * alpha * beta;
    print(np.sum(det_n1[1:]))

#get_DR_IP()



