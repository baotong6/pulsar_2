#!/usr/bin/env python
# coding: utf-8

# In[34]:


import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import pandas as pd
import sys
import os
from scipy import interpolate

# In[2]:

font1 = {'family': 'Normal',
         'weight': 'normal',
         'size': 15, }
path="/Users/baotong/Desktop/period_Tuc/"
type=['47Tuc','terzan5','M28','omg_cen']


# In[45]:


def add_pos_for_excel(label):
    if label=='omg_cen':dir=label[0:-4]
    else: dir=label
    res=pd.read_excel(path+'result_0.5_8_all.xlsx',label)
    #srclis=fits.open('/Users/baotong/Desktop/period_{0}/{1}_p50_i5_src_1_2_4_8.fits'.format(dir,label))
    srclis=fits.open('/Users/baotong/Desktop/period_Tuc//xray_properties-592.fits'.format(dir))
    #for 47Tuc
    ra=srclis[1].data['RAdeg']
    dec=srclis[1].data['DEdeg']
    seq=np.linspace(1,len(ra),len(ra))
    seq=seq.astype(np.int)
    result_out=np.column_stack((ra[res['seq']-1],dec[res['seq']-1]))
    np.savetxt(path+'radec_{0}.txt'.format(label),result_out,fmt="%10.5f %10.5f")


# In[46]:


# add_pos_for_excel(type[0])

# # 对已有的四个星团的变源src list做一个总region

# In[52]:


def make_region(label):
    res=pd.read_excel(path+'result_0.5_8_all.xlsx',label)
    os.system('rm {0}'.format(path)+'all_pCV_{0}.reg'.format(label))
    with open(path+'all_pCV_{0}.reg'.format(label),'a+') as f2:
        f2.writelines('fk5'+'\n')
        
    ra=res['RA']
    dec=res['DEC']
    seq=res['seq']
    for i in range(len(ra)):
        with open(path+'all_pCV_{0}.reg'.format(label),'a+') as f2:
            reg = 'circle(' + str(ra[i]) + ',' + str(dec[i]) + ',' + str('2"') + ')'+" # color=green width=2 text={"+str(seq[i])+"}"
            f2.writelines(reg+'\n')

# make_region(type[1])
# # 看一下epoch里的观测时间如何，方便plot long-term 时选择x轴范围

def print_MJD():
    path='/Volumes/pulsar/M28/merge_data/spectra/aprates/'
    EPOCH=np.loadtxt(path + 'M28_epoch.txt')
    obs_time=(EPOCH[:, 0] + EPOCH[:, 1]) / 2
    time = obs_time / 86400 + 2449352.5 - 2400000.5
    print(time)
#print_MJD()

def plot_profile_distb():
    # label_all=['47Tuc','terzan5','M28','omg_cen']
    label_all=['47Tuc','terzan5','M28','omg_cen','NGC6397','NGC6752']
    pos_all=[[6.0236250,-72.0812833,3.17*60,3.17/8.8*60],  #47Tuc
             [267.0202083,-24.7790556,0.72*60,0.72/3.4*60],  #terzan5
             [276.1363750,-24.8702972,1.97*60,1.97/8.2*60], # M28
             [201.69700,-47.47947 , 5*60,5/2.1*60],         #omega_cen
            [265.17539,-53.67433,2.9*60,2.9/58*60],        #NGC 6397
             [287.71713,-59.98455,1.91,1.91/11.24*60]]      #NGC 6752
    distance_all=[0]
    period_all=[0];L_all=[0];Lmin_all=[0];Lmax_all=[0]
    for i in range(5):
        label=label_all[i]
        pos=pos_all[i]
        res=pd.read_excel(path+'result_0.5_8_all.xlsx',label)
        ra=np.array(res['RA'])
        dec=np.array(res['DEC'])
        seq=np.array(res['seq'])
        period=np.array(res['P_out'])
        L=np.array(res['L'])
        Lmin=np.array(res['Lmin'])
        Lmax=np.array(res['Lmax'])
        type=np.array(res['type'])
        i=0
        while i<len(seq):
            if period[i]<4500. or type[i]!='CV':
            #if period[i]<4500:
                ra=np.delete(ra,i)
                dec=np.delete(dec, i)
                seq=np.delete(seq, i)
                period=np.delete(period, i)
                type=np.delete(type,i)
                L = np.delete(L, i);Lmin = np.delete(Lmin, i);Lmax = np.delete(Lmax, i)
            else:i+=1;
        #print(seq)
        #print((((ra-pos[0])*3600)**2+((dec-pos[1])*3600)**2)**0.5)
        distance=(((ra-pos[0])*3600)**2+((dec-pos[1])*3600)**2)**0.5/pos[-1]
        distance_all=np.concatenate((distance_all,distance))
        period_all=np.concatenate((period_all,period))
        L_all=np.concatenate((L_all,L));Lmin_all=np.concatenate((Lmin_all,L));Lmax_all=np.concatenate((Lmax_all,L))
    #print(distance_all)
    distance_all=distance_all[1:]
    period_all=period_all[1:]
    L_all=L_all[1:];Lmin_all=Lmin_all[1:];Lmax_all=Lmax_all[1:]
    period_plot=[]
    num_density=[]
    num_periodbin = []
    period_all=np.array(period_all)
    print(len(period_all))
    #plt.ylim(1e-4,1e2)
    #bins=np.linspace(0,2,16)
    def plot_P_Lum():
        label_all = ['47Tuc', 'terzan5', 'M28', 'omg_cen','NGC6397','NGC6752']
        def read_P_L(label):
            res = pd.read_excel(path + 'result_0.5_8_all.xlsx', label)
            type=np.array(res['type'])
            period=np.array(res['P_out'])
            L=np.array(res['L'])
            Lmin=np.array(res['Lmin'])
            Lmax=np.array(res['Lmax'])
            distance=np.array(res['distance'])
            i = 0
            while i < len(type):
                if (type[i] != 'CV'):
                    #print(type[i])
                    # if str(seq[i])[-3:]=='001':
                    period = np.delete(period, i)
                    type = np.delete(type, i)
                    L = np.delete(L, i);
                    Lmin = np.delete(Lmin, i);
                    Lmax = np.delete(Lmax, i)
                    distance=np.delete(distance,i)
                else:
                    i += 1;
            return (period,L,distance)

        PL_1 = read_P_L('47Tuc')
        PL_2 = read_P_L('terzan5')
        PL_3 = read_P_L('M28')
        PL_4 = read_P_L('omg_cen')
        PL_5= read_P_L('NGC6397')
        # PL_6 = read_P_L('NGC6752')
        plt.figure(1,(10,8))

        plt.subplot(211)
        plt.semilogy()
        plt.xlim(-0.2,15)
        plt.ylim(5e30,1e33)
        plt.scatter(PL_1[0] / 3600., PL_1[1],marker='^',color='w',linewidths=2,edgecolors='r')
        plt.scatter(PL_2[0] / 3600., PL_2[1],marker='X',color='green',s=50)
        plt.scatter(PL_3[0] / 3600., PL_3[1],marker='s')
        plt.scatter(PL_4[0] / 3600., PL_4[1],marker='*')
        plt.scatter(PL_5[0] / 3600., PL_5[1], marker='v',color='black')

        plt.legend(label_all)
        P_min = 4944.0 / 3600.
        P_gap = [7740.0 / 3600., 11448.0 / 3600.]
        x = P_gap
        y = [0, 1e36]
        plt.text(7900 / 3600., 5e36, 'period gap')
        plt.fill_between(x, y[1], facecolor='yellow', alpha=0.2)
        plt.semilogy()

        plt.plot([P_min, P_min], [0, 1e36], '--')
        plt.text(P_min - 0.5, 1e36, 'period minum')
        plt.xlabel('Period(h)',font1)
        plt.ylabel('Luminosity (erg/s)',font1)
        plt.tick_params(labelsize=15)

        plt.subplot(212)
        plt.xlim(-0.2,15)
        plt.semilogy()
        plt.ylim(1e-2,12)
        plt.plot([0,15],[1,1],'--',color='grey')
        plt.scatter(PL_1[0] / 3600., PL_1[2], marker='^', color='w', linewidths=2, edgecolors='r')
        plt.scatter(PL_2[0] / 3600., PL_2[2], marker='X', color='green', s=50)
        plt.scatter(PL_3[0] / 3600., PL_3[2], marker='s')
        plt.scatter(PL_4[0] / 3600., PL_4[2], marker='*')
        plt.scatter(PL_5[0] / 3600., PL_5[2], marker='v',color='black')

        P_min = 4944.0 / 3600.
        P_gap = [7740.0 / 3600., 11448.0 / 3600.]
        x = P_gap
        y = [0, 12]
        plt.fill_between(x, y[1], facecolor='yellow', alpha=0.2)
        plt.plot([P_min, P_min], [0, 12], '--')
        plt.xlabel('Period(h)',font1)
        plt.ylabel('Distance (/half-light radius)',font1)
        plt.tick_params(labelsize=15)
        figurepath='/Users/baotong/Desktop/aas/pCV_GC/figure/'
        plt.savefig(figurepath+'L_PandR_P.eps')
        plt.show()

    plot_P_Lum()

    bins=np.linspace(0,5,15)
    for i in range(len(bins)-1):
        period_bin=period_all[np.where((distance_all<bins[i+1])&(distance_all>=bins[i]))]
        num_density.append(len(period_bin)/(np.pi*(bins[i+1]**2-bins[i]**2)))
        #print(len(period_bin),np.pi*(bins[i+1]**2-bins[i]**2))
        period_acm=np.mean(period_bin)
        period_plot.append(period_acm)
        num_periodbin.append(len(period_bin))

    # num_periodbin=np.array(num_periodbin)
    # x=bins[0:-1]+(bins[1]-bins[0])/2.
    # y=np.array(num_density)
    # func = interpolate.interp1d(x, y, kind='cubic')
    # x_new = np.linspace(min(x), max(x),100)
    # y_new = func(x_new)
    # #plt.semilogy()
    # #plt.ylim(1e-2,1e3)
    # #---num_density---#####
    # plt.scatter(x, y, color='red')
    # plt.plot(x_new, y_new,color='green',linestyle='--')
    #----distance distribution-----#####
    # plt.hist(distance_all,bins=20,histtype='step')

    # ##---period_mean---#####
    # plt.scatter(x, period_plot, s=num_periodbin * 10, marker='o', edgecolors='r')  # period_mean
    # num_periodbin=num_periodbin.astype('str')
    # print(num_periodbin)
    # for i in range(len(x)):
    #     plt.text(x[i],period_plot[i]+1000,num_periodbin[i])
    # ##---------#####

    #plt.scatter(period_all,distance_all)

    #----period_hist_r------#####
    ###按照距离取bin
    # plt.scatter(x,period_plot)  # period_mean
    #
    # P_gap = [7740.0, 11448.0]
    # plt.plot([P_gap[0], P_gap[0]],[0,10], '--', lw=2., color='yellow')
    # plt.plot([P_gap[1], P_gap[1]],[0,10], '--', lw=2., color='yellow')
    # print(len(period_all))
    # plt.hist(period_all,bins=10,histtype='step',lw=3.,color='purple')  # all period
    # plt.plot([pos[-1], pos[-1]], [1000,50000], '--')

    #
    ## #-------period_hist_p--------###
    ## #按照周期取bin
    i=0
    period_bar=[];distance_bar=[]
    Z = zip(distance_all,period_all)
    Z2= zip(period_all,distance_all)
    Z2 = sorted(Z2, reverse=False)
    Z=sorted(Z, reverse=False)
    #distance_all, period_all=zip(*Z)
    period_all, distance_all = zip(*Z2)

    # leng=5  ##每个bin里的源个数
    # #
    # while i<len(period_all):
    #     if i+2*leng<len(period_all):
    #         period_bar.append(np.mean(period_all[i:i+leng]))
    #         distance_bar.append(np.mean(distance_all[i:i+leng]))
    #     else:
    #         period_bar.append(np.mean(period_all[i:]))
    #         distance_bar.append(np.mean(distance_all[i:]))
    #     i+=leng
    # period_bar=np.array(period_bar);distance_bar=np.array(distance_bar)
    # plt.xlabel('Period mean (hour)')
    # plt.ylabel('Distance mean (/half-light radius)')
    # plt.scatter(period_bar/3600.,distance_bar,color='red')

    #
    plt.show()

plot_profile_distb()





