import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import pandas as pd
import sys
import os
from scipy import interpolate
from scipy.optimize import curve_fit
import astropy.units as u
import astropy.constants as c
from scipy import interpolate
import stingray as sr
import useful_functions as func

font1=func.font1
figurepath='/Users/baotong/Desktop/aas/AGN_CDFS/figure/'

def plot_scatter_single(filename,qpo_P,threshold=0.0027):
    temp=np.loadtxt(filename)
    FP=temp[:,0];period=temp[:,1];
    if qpo_P:
        id=np.where((FP<threshold)&(np.abs(period-qpo_P)<1/16*qpo_P))[0]
        DR = len(id) / 1000
    else:
        id = np.where((FP < threshold))[0]
        DR = len(id) / 1000  ##这里的DR为fDR
    plt.scatter(period,FP,marker='x',color='green')
    plt.loglog()
    plt.show()
    print(DR)
    return DR
def run_test():
    epoch1 = '/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep1/CDFS_epoch_ep1.txt'
    epoch2 = '/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep2/CDFS_epoch_ep2.txt'
    epoch3 = '/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep3/CDFS_epoch_ep3.txt'
    epoch4 = '/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep4/CDFS_epoch_ep4.txt'
    epoch_all=[epoch1,epoch2,epoch3,epoch4]
    CR_all=[5e-4,6e-4,7e-4,8e-4,9e-4,1e-3]
    period_all=[1800,3600,7200]
    ep=3;CR=CR_all[5];period=period_all[2]
    path = '/Users/baotong/Desktop/CDFS/simulation/EP{0}/'.format(int(ep + 1))
    # filename=path+'CR_{0}_P_{1}_REJ1034.txt'.format("%.0e"%CR,str(int(period)))

    filename = path + 'CR_{0}_P_{1}_REJ1034.txt'.format("%.0e" % CR, str(int(period)))
    # filename = path + 'CR_{0}_noQPO_REJ1034.txt'.format("%.0e" % CR)
    # filename = path + 'CR_1e-03_P_1800_REJ1034_three_way.txt'.format("%.0e" % CR, str(int(period)))
    plot_scatter_single(filename,qpo_P=None,threshold=0.0027)


def read_one_ep(CR_all,period_all,ep,threshold=0.01):
    epoch='/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep{0}/CDFS_epoch_ep{0}.txt'.format(ep)
    path = '/Users/baotong/Desktop/CDFS/simulation/EP{0}/'.format(ep)
    DR=np.zeros((len(CR_all),(len(period_all))))
    for i in range(len(CR_all)):
        CR = CR_all[i];
        for j in range(len(period_all)):
            # print(i,j)
            CR = CR_all[i];
            period = period_all[j]
            filename = path + 'CR_{0}_P_{1}_REJ1034.txt'.format("%.0e" % CR, str(int(period)))
            temp = np.loadtxt(filename)
            FP = temp[:, 0];
            P_det = temp[:, 1];
            id = np.where((FP < threshold) & (np.abs(P_det - period) < 1 / 16 * period))[0]
            DR[i][j]=len(id)
    return DR
def plot_sim_one_ep(ep,CR_all,period_all,threshold=0.01,):
    DR=read_one_ep(CR_all,period_all,ep,threshold)
    for j in range(len(period_all)):
        #用百分制
        plt.plot(CR_all*1e5,DR[:,j]/10,marker='v', linestyle='-')
    plt.xlabel('Photon flux ($10^{-5}$ counts $s^{-1}$ )', font1)
    plt.ylabel('Detection rate(%)', font1)
    plt.tick_params(labelsize=16)
    plt.legend(['P=0.5h','P=1h','P=2h'])
    plt.show()

def plot_fDR_confirm(CR_all,ep):
    path = '/Users/baotong/Desktop/CDFS/simulation/EP{0}/'.format(ep)
    threshold_range = np.linspace(0, 1, 50)
    FP_ideal = np.zeros((len(CR_all),len(threshold_range)))
    for i in range(len(CR_all)):
        CR = CR_all[i];
        filename=path + 'CR_{0}_noqpo_REJ1034.txt'.format("%.0e" % CR)
        temp=np.loadtxt(filename)
        FP=temp[:,0]
        for j in range(len(threshold_range)):
            FP_ideal[i][j]=len(np.where(FP<threshold_range[j])[0])

    FP_ideal/=1000
    for k in range(len(CR_all)):
        plt.plot(threshold_range,FP_ideal[k],'-.',color='green')
    plt.legend(['$\\bar{\lambda}$=5,6,7,8,9,10 x $10^{-4}$ counts$~$$s^{-1}$'],fontsize=14)
    plt.plot([0,1],[0,1],'-',color='red')
    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.xlabel('The Value of FAP',font1)
    plt.ylabel('Cumulative distribution function',font1)
    plt.tick_params(labelsize=16)
    plt.savefig(figurepath+'cdf_FAP.eps',bbox_inches='tight', pad_inches=0.0)
    plt.show()

def read_source_info():
    path='/Users/baotong/Desktop/CDFS/'
    info = pd.read_excel(path + 'source_infomation.xlsx')

def compute_expect_Det(CR_all,period_all,k_num,threshold):
    path='/Users/baotong/Desktop/CDFS/'
    info = pd.read_excel(path + 'source_infomation.xlsx')
    for k in k_num:
        DR = read_one_ep(CR_all, period_all, k, threshold)
        DR=DR/1000.
        # label='CR_ep'+str(k)
        label='CR'
        CR=info[label]
        CR=np.array(CR)
        CR_use=CR[np.where((CR>=CR_all[0])&(CR<=CR_all[-1]))]
        ##在simulation范围内的源的expect探测数目##
        print(len(CR_use))
        DR_meanP=DR.mean(axis=1)
        ##插个值##
        x=CR_all;y=DR_meanP
        kind=["nearest","zero","slinear","quadratic","cubic"]
        f = interpolate.interp1d(x, y, kind='cubic')
        DR_use=f(CR_use)

        CR_high_use=CR[np.where(CR>CR_all[-1])]
        print(len(CR_high_use))
        ##在simulation范围之上的源的expect探测数目##
        z1 = np.polyfit(x, y, 1)
        p1 = np.poly1d(z1)
        DR_high_use=p1(CR_high_use)
        DR_high_use[np.where(DR_high_use>1)]=1
        # print(DR_high_use)

        DR_use=np.concatenate((DR_use,DR_high_use))
        print(np.sum(DR_use))


def plot_all_ep(CR_all,period_all,k_num,threshold):
    figlabel = [[0, 0], [0, 1], [1, 0], [1, 1]]
    fig, axes = plt.subplots(2, 2,figsize=(15,10))
    for i in range(len(k_num)):
        k = k_num[i];
        DR= read_one_ep(CR_all,period_all,k, threshold)
        ax_temp = axes[figlabel[i][0], figlabel[i][1]]
        x = CR_all * 1e5
        for j in range(len(period_all)):
            y1=DR[:,j]/10
            ax_temp.plot(x, y1, marker='v', linestyle='-')

        ax_temp.text(50, 40, 'Epoch{0}'.format(k), font1)
        ax_temp.legend(['DR P=0.5h','DR P=1h','DR P=2h'],loc='center left')
        ax_temp.set_ylim(0,45)
        if i<3:ax_temp.legend_.remove()
        # ax_temp.set_title('Epoch {0}: LS detection results'.format(k),font1)
        if (i==2 or i==3):ax_temp.set_xlabel('Photon flux ($10^{-5}$ counts$~$$s^{-1}$ )',font1)
        if (i==0 or i==2):ax_temp.set_ylabel('Detection rate (%)',font1)
        ax_temp.tick_params(labelsize=16)
    plt.savefig(figurepath + 'DR_thres_{0}.eps'.format(1 - threshold), bbox_inches='tight', pad_inches=0.0)
    plt.show()

if __name__=='__main__':
    run_test()
    # read_source_info()
    CR_all=np.array([5e-4,6e-4,7e-4,8e-4,9e-4,1e-3])
    period_all=np.array([1800,3600,7200])
    # compute_expect_Det(CR_all, period_all, k_num=[1,2,3,4], threshold=0.0027)
    # plot_fDR_confirm(CR_all,ep=4)
    # plot_all_ep(CR_all, period_all, k_num=[1,2,3,4], threshold=0.0027)
    # plot_sim_one_ep(1, threshold=0.0027)
    # plot_sim_one_ep(2,threshold=0.0027)
    # plot_sim_one_ep(3, threshold=0.0027)
    # plot_sim_one_ep(4, threshold=0.0027)

