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
import read_data as data
#import read_data_gc as data
#import read_data_lw as data

path='/Users/baotong/Desktop/period/'
#dataname="src2338_I_lc.txt"
dataname='1671.txt'
#epoch_file=path+'txt_90/'+'SgrA_I_epoch.txt'
#epoch_file=path+'txt_2_8k/'+'SgrA_I_epoch.txt'
epoch_file=path+'txt_all_obs_IG/'+'epoch_src_'+dataname

#epoch_file = path + 'txt_G/' + 'SgrA_G_epoch.txt'
#epoch_file =path + 'txt_merge_2_8k/' + 'SgrA_IG_epoch.txt'

p_test=48016.90195
bin=20
shift=0.0
net_percent=0.716
src_bkg=1-net_percent

##tot cts/bkg cts
#dataname='Magnetar_S_lc.txt'
time=data.get_data(dataname)[0]
energy=data.get_data(dataname)[1]
obs_ID=data.get_data(dataname)[2]
t=data.get_data(dataname)[3]
E=data.get_data(dataname)[4]
dict=data.get_data(dataname)[5]
#time=t[0]
N=len(time)
# time=np.concatenate((dict[242][0],dict[2951][0],dict[2952][0],dict[2953][0],dict[2954][0],
#                      dict[2943][0],dict[3663][0],dict[3392][0],dict[3393][0],dict[3665][0]))
delete_i=[]
use_ID=[]
#dict[242,2951][0]
for i in dict:
    if dict[i][-2]<5000:
        delete_i.append(i)
for item in delete_i:
    dict.pop(item)
for i in dict:
    use_ID.append(i)
#print use_ID
#use_ID=[10556,11843]
#print dict[10556]

epoch0=[3392,3393]
epoch_01=[2943,3663,3665]
epoch1=[2943,3663,3392,3393,3665]
epoch2=[5950,5951,5952,5953,5954]
epoch3=[9169,9170,9171,9172,9174, 9173,10556,11843,13016,13017,14941,14942]
epoch4=[3392,3393]
#
# epoch_gc = np.array([14897, 17236, 17239, 17237, 18852, 17240, 17238, 20118, 17241, 20807, 20808])
#use_ID = epoch_gc
# epoch1 = np.array([13847, 14427, 13848, 13849, 13846, 14438, 13845])
# epoch2 = np.array([14461, 13853, 13841, 14465, 14466, 13842, 13839, 13840, 14432, 13838, 13852, 14439])
# epoch3 = np.array([14462, 14463, 13851, 15568, 13843, 15570, 14468])
#use_ID = np.array([3549])
#use_ID=np.concatenate((epoch2,epoch3))
if len(use_ID) < 1:
    print('error')
elif len(use_ID) == 1:
    time = dict[use_ID[0]][0]
else:
    time = np.concatenate((dict[use_ID[0]][0], dict[use_ID[1]][0]))

    for i in range(2, len(use_ID)):
        time = np.concatenate((time, dict[use_ID[i]][0]))

vdot=0
#根据表来查
num_p=1.0
# p_test=49405.45
def trans(t,num_p,p_test,resolution,border):
    ti=t
    #print ti
    vary=[i for i in range(-border,border)]
    vary=np.array(vary)
    #p_test = 1.0/5.500550055005501e-06
    #p_test=4945.45
    p_test=num_p*p_test
    v = 1.0 / (p_test + vary *resolution)
    p=1.0/v
    # pdot=-vdot/(v*v)
    # vddot = 2.0 * pdot * pdot / (p * p * p)
    freq = v # + ti * vdot + ti * ti / 2.0 * vddot
    preq = 1.0 / freq
    turns = v * ti #+ vdot * ti * ti / 2.0 + vddot * ti * ti * ti / 6.0
    turns += shift
    #初始相位
    for i in range(len(turns)):
        turns[i]=turns[i] - int(turns[i])
    turns = num_p * turns
    return turns

def getIndexes(y_predict, y_data):
    n = y_data.size
    # SSE为和方差
    SSE = ((y_data - y_predict) ** 2).sum()
    # MSE为均方差
    MSE = SSE / n
    # RMSE为均方根,越接近0，拟合效果越好
    RMSE = np.sqrt(MSE)

    # 求R方，0<=R<=1，越靠近1,拟合效果越好
    u = y_data.mean()
    SST = ((y_data - u) ** 2).sum()
    SSR = SST - SSE
    R_square = SSR / SST
    return SSE, MSE, RMSE, R_square

def func(x,a,b,c):
    return a*np.sin(2*np.pi*x+b)+c

border=100
resolution=0.1
def phasefold(time,num_p,bin,p_test):
    turns_n=[]
    for index in time:
        turns_n.append(trans(index,num_p,p_test,border=100,resolution=0.1))
    turns_n=np.transpose(turns_n)
    #bin=20
    #print len(turns_n[0])
    sort_turn=[]
    for index in turns_n:
        sort_turn.append(np.sort(index))
    #print len(sort_turn)
    loc=np.zeros([2*border,bin])
    for i in range(2*border):
        for index in sort_turn[i]:
            loc[i][int(index / (num_p / bin))] += 1

    x = np.array([(num_p * i / bin + 0.5 / bin) for i in range(bin)])
    return [x,loc]

def error_sin(x,loc):
    error_MSE=[]
    cof=[]
    for i in range(len(loc)):
        ydata=loc[i]
        popt, pcov = curve_fit(func, x, ydata,bounds=(0, [1000., 2*np.pi, 3000.]))
        y_predict = func(x, *popt)
        indexes_1 = getIndexes(y_predict, ydata)
        error_MSE.append(indexes_1[1])
        cof.append(popt)
    return [error_MSE,cof]

def get_T_in_mbins(epoch_file,w,m,fi):
    T=2*np.pi/w
    T_in_perbin = np.zeros(m)
    # 每个bin的总积分时间
    tbin = T/m
    # 每个bin的时间长度
    epoch_info = np.loadtxt(epoch_file)
    t_start = epoch_info[:, 0]
    t_end = epoch_info[:, 1]
    ID = epoch_info[:, 2]
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
T_in_perbin=get_T_in_mbins(epoch_file,2*np.pi/p_test,bin,shift*2*np.pi)

x=phasefold(time,num_p=1.0,bin=bin
            ,p_test=p_test)[0]
loc=phasefold(time,num_p=1.0,bin=bin,p_test=p_test)[1]
error_MSE = error_sin(x, loc)[0]
cof = error_sin(x, loc)[1]
print('done')
# plt.figure(1)
# plt.step(x, loc[error_MSE.index(min(error_MSE))])
# plt.plot(x, func(x, *cof[error_MSE.index(min(error_MSE))]))
# plt.step(x, loc[border])
# plt.plot(x, func(x, *cof[border]))
# plt.legend(["min","min_sin","original","original_sin"])
# # print p_test+(error_MSE.index(min(error_MSE))-border)*resolution/num_p
# print(error_MSE.index(min(error_MSE)))
# plt.figure(1,(15,10))
plt.figure(1)
x1=x+1
x2=np.concatenate((x,x1))
y2=np.concatenate((loc[border],loc[border]))
plt.xlabel('phase')
plt.ylabel('counts/bin')
bkg_y=len(time)*src_bkg
bkg_y/=bin
print(bkg_y)
plt.plot([0,2],[bkg_y,bkg_y],'--')

#plt.step(x2,y2/(p_test/50.))
correct_gap=T_in_perbin/(sum(T_in_perbin)/len(T_in_perbin))
print(correct_gap)
y2/=np.concatenate((correct_gap,correct_gap))
plt.step(x2,y2,color='red')
#plt.scatter(x2,y2)
plt.errorbar(x2-0.5/bin, y2, yerr=y2**0.5,fmt='.',capsize=1,elinewidth=1,ecolor='red')
#plt.plot(x2,func(x2,*cof[border])/(p_test/50))
per=p_test
plt.ylim(ymin=bkg_y-10)
plt.title("{0},P={1}".format(dataname[0:-4],str(per)),fontsize=18)
#plt.text(0.5,max(loc[border])-5,"p="+str(per)+'\n'+"counts="+str(len(time)),fontsize=15)
plt.text(0.8,bkg_y-5,"counts="+str(len(time)))
#plt.savefig('/Users/baotong/Desktop/li_pr/result_final/fig_NSC_I/pfold_lc_{0}.eps'.format(dataname[0:-4]))
plt.show()

'''
p_possible=np.loadtxt(path+dataname[3:7]+"/p_possible.txt")
p_possible=np.array(p_possible)
#p_possible=np.array([p_possible])
p_possible=np.concatenate((p_possible,[29879.7]))
print p_possible
for p_test in p_possible:


    x=phasefold(time,num_p=2.0,bin=100,p_test=p_test)[0]
    loc=phasefold(time,num_p=2.0,bin=100,p_test=p_test)[1]
    error_MSE=error_sin(x,loc)[0]
    cof=error_sin(x,loc)[1]
    plt.figure(1)
    plt.step(x,loc[error_MSE.index(min(error_MSE))])
    plt.plot(x,func(x, *cof[error_MSE.index(min(error_MSE))]))
    plt.step(x,loc[border])
    plt.plot(x,func(x,*cof[border]))
    plt.legend(["min","min_sin","original","original_sin"])
    plt.savefig(path+data.dataname[3:7]+"/"+str(int(p_test))+"compare_phase_2p.eps")
    plt.close(1)

    plt.figure(2)
    plt.step(x,loc[error_MSE.index(min(error_MSE))])
    plt.plot(x,func(x, *cof[error_MSE.index(min(error_MSE))]))
    per=p_test+((error_MSE.index(min(error_MSE))-border)*resolution)/num_p
    plt.text(0.5,max(loc[error_MSE.index(min(error_MSE))]),"p="+str(per))
    plt.savefig(path+data.dataname[3:7]+"/"+str(int(p_test))+"min_phase_2p.eps")
    plt.close(2)

    plt.figure(3)
    plt.step(x,loc[error_MSE.index(min(error_MSE))])
    plt.plot(x,func(x, *cof[error_MSE.index(min(error_MSE))]))
    per=p_test
    plt.text(0.5,max(loc[error_MSE.index(min(error_MSE))]),"p="+str(per))
    plt.savefig(path+data.dataname[3:7]+"/"+str(int(p_test))+"original_phase_2p.eps")
    plt.close(3)


    x=phasefold(time,num_p=1.0,bin=50,p_test=p_test)[0]
    loc=phasefold(time,num_p=1.0,bin=50,p_test=p_test)[1]
    error_MSE=error_sin(x,loc)[0]
    cof=error_sin(x,loc)[1]
    plt.figure(4)
    plt.step(x,loc[error_MSE.index(min(error_MSE))])
    plt.plot(x,func(x, *cof[error_MSE.index(min(error_MSE))]))
    plt.legend(["min","min_sin"])
    per=p_test+((error_MSE.index(min(error_MSE))-border)*resolution)/num_p
    plt.text(0.5,max(loc[error_MSE.index(min(error_MSE))]),"p="+str(per))
    plt.savefig(path+data.dataname[3:7]+"/"+str(int(p_test))+"min_phase.eps")
    plt.close(4)

    x1=x+1
    x2=np.concatenate((x,x1))
    y2=np.concatenate((loc[error_MSE.index(min(error_MSE))],loc[error_MSE.index(min(error_MSE))]))
    plt.figure(5)
    plt.step(x2,y2)
    plt.plot(x2,func(x2, *cof[error_MSE.index(min(error_MSE))]))
    per=p_test+((error_MSE.index(min(error_MSE))-border)*resolution)/num_p
    plt.text(0.5,max(loc[error_MSE.index(min(error_MSE))]),"p="+str(per))
    plt.savefig(path+data.dataname[3:7]+"/"+str(int(p_test))+"min_fake_p2.eps")
    plt.close(5)


    x=phasefold(time,num_p=1.0,bin=50,p_test=p_test)[0]
    loc=phasefold(time,num_p=1.0,bin=50,p_test=p_test)[1]
    error_MSE=error_sin(x,loc)[0]
    cof=error_sin(x,loc)[1]
    plt.figure(6)
    plt.step(x,loc[border])
    plt.plot(x,func(x,*cof[border]))
    plt.legend(["original","original_sin"])
    per=p_test
    plt.text(0.5,max(loc[border]),"p="+str(per))
    plt.savefig(path+data.dataname[3:7]+"/"+str(int(p_test))+"original_phase.eps")
    plt.close(6)


    # x=phasefold(time,num_p=1.0,bin=30)[0]
    # loc=phasefold(time,num_p=1.0,bin=30)[1]
    x1=x+1
    x2=np.concatenate((x,x1))
    y2=np.concatenate((loc[border],loc[border]))
    plt.figure(7)
    plt.step(x2,y2)
    plt.plot(x2,func(x2,*cof[border]))
    per=p_test
    plt.text(0.5,max(loc[border]),"p="+str(per))
    plt.savefig(path+data.dataname[3:7]+"/"+str(int(p_test))+"original_fake_p2.eps")
    plt.close(7)

    print error_MSE.index(min(error_MSE))
    #原始数据

    # plt.legend(["original","original_sin"])
    # plt.show()
'''