import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.patches import ConnectionPatch
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
from CDFS.CDFS_startover import useful_functions as func
from CDFS.CDFS_startover import sim_psd as sim


import matplotlib.patches as mpatches
font1=func.font1
figurepath='/Users/baotong/Desktop/aas/AGN_CDFS_mod1/figure/'

def plot_scatter_single(filename,qpo_P,threshold=0.0027):
    temp=np.loadtxt(filename)
    FP=temp[:,0];period=temp[:,1];
    if qpo_P:
        id=np.where((FP<threshold)&(np.abs(period-qpo_P)<1/10*qpo_P))[0]
        DR = len(id) / 1000
    else:
        id = np.where((FP < threshold))[0]
        DR = len(id) / 1000  ##这里的DR为fDR
    plt.scatter(period,FP,marker='x',color='green')
    plt.loglog()
    plt.show()
    # print(DR)
    return DR

def run_test():
    epoch1 = '/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep1/CDFS_epoch_ep1.txt'
    epoch2 = '/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep2/CDFS_epoch_ep2.txt'
    epoch3 = '/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep3/CDFS_epoch_ep3.txt'
    epoch4 = '/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep4/CDFS_epoch_ep4.txt'
    epoch_all=[epoch1,epoch2,epoch3,epoch4]
    CR_all=[5e-4,6e-4,7e-4,8e-4,9e-4,1e-3]
    period_all=[3600,5400,7200]
    ep=1;CR=CR_all[5];period=period_all[2]
    path = '/Users/baotong/Desktop/CDFS/simulation/EP{0}/'.format(int(ep + 1))
    # filename=path+'CR_{0}_P_{1}_REJ1034.txt'.format("%.0e"%CR,str(int(period)))

    filename = path + 'CR_{0}_P_{1}_A002_R15.txt'.format("%.0e" % CR, str(int(period)))
    # filename = path + 'CR_{0}_noQPO_REJ1034.txt'.format("%.0e" % CR)
    # filename = path + 'CR_1e-03_P_1800_REJ1034_three_way.txt'.format("%.0e" % CR, str(int(period)))
    plot_scatter_single(filename,qpo_P=None,threshold=0.0027)

def read_one_ep(CR_all,period_all,ep,threshold=0.01,label='REJ1034'):
    epoch='/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep{0}/CDFS_epoch_ep{0}.txt'.format(ep)
    path = '/Users/baotong/Desktop/CDFS/simulation/EP{0}/'.format(ep)
    DR=np.zeros((len(CR_all),(len(period_all))))
    for i in range(len(CR_all)):
        CR = CR_all[i];
        for j in range(len(period_all)):
            # print(i,j)
            CR = CR_all[i];
            period = period_all[j]
            # filename = path + 'CR_{0}_P_{1}_REJ1034_R_5_Q_15.txt'.format("%.0e" % CR, str(int(period)))
            # filename = path + 'CR_{0}_P_{1}_REJ1034_const1Ms.txt'.format("%.0e" % CR, str(int(period)))
            # filename = path + 'CR_{0}_P_{1}_REJ1034.txt'.format("%.0e" % CR, str(int(period)))
            filename = path + 'CR_{0}_P_{1}_{2}.txt'.format("%.0e" % CR, str(int(period)),label)
            # filename = path + 'CR_{0}_P_{1}_A002_R15.txt'.format("%.0e" % CR, str(int(period)))
            temp = np.loadtxt(filename)
            # print(temp[0])
            FP = temp[:, 0];
            P_det = temp[:, 1];
            id = np.where((FP < threshold) & (np.abs(P_det - period) < 1 / 16 * period))[0]
            DR[i][j]=len(id)
    return DR

def read_one_ep_noqpo(CR_all,ep,threshold=0.01,label='REJ1034'):
    path = '/Users/baotong/Desktop/CDFS/simulation/EP{0}/'.format(ep)
    threshold_range = np.linspace(0, 1, 50)
    FP_ideal = np.zeros((len(CR_all),len(threshold_range)))
    fDR=np.zeros(len(CR_all))
    for i in range(len(CR_all)):
        CR = CR_all[i];
        filename=path + 'CR_{0}_noqpo_{1}.txt'.format("%.0e" % CR,label)
        temp=np.loadtxt(filename)
        FP=temp[:,0]
        id = np.where((FP < threshold))[0]
        fDR[i]=len(id)
    return fDR

def plot_sim_one_ep(ep,CR_all,period_all,threshold=0.01,label='REJ1034'):
    DR=read_one_ep(CR_all,period_all,ep,threshold,label=label)
    for j in range(len(period_all)):
        #用百分制
        plt.plot(CR_all*1e5,DR[:,j]/10,marker='v', linestyle='-')
    plt.xlabel('Photon flux ($10^{-5}$ counts $s^{-1}$ )', font1)
    plt.ylabel('Detection rate(%)', font1)
    plt.tick_params(labelsize=16)
    # plt.legend(['P=0.5h','P=1h','P=2h'])
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
        plt.plot(threshold_range,FP_ideal[k],'-.')

    plt.legend(['$\\bar{x}$=5 x $10^{-4}$ cts$~$$s^{-1}$',
                '$\\bar{x}$=6 x $10^{-4}$ cts$~$$s^{-1}$',
                '$\\bar{x}$=7 x $10^{-4}$ cts$~$$s^{-1}$',
                '$\\bar{x}$=8 x $10^{-4}$ cts$~$$s^{-1}$',
                '$\\bar{x}$=9 x $10^{-4}$ cts$~$$s^{-1}$',
                '$\\bar{x}$=10 x $10^{-4}$ cts$~$$s^{-1}$'], fontsize=10)

    plt.plot([0,1],[0,1],'-',color='grey',linewidth=2)
    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.xlabel('FAP',font1)
    plt.ylabel('Cumulative distribution function',font1)
    plt.tick_params(labelsize=16)
    # plt.savefig(figurepath+'cdf_FAP.eps',bbox_inches='tight', pad_inches=0.0)
    # plt.savefig(figurepath + 'cdf_FAP.pdf', bbox_inches='tight', pad_inches=0.0)
    plt.show()

def read_source_info():
    path='/Users/baotong/Desktop/CDFS/'
    info = pd.read_excel(path + 'source_infomation.xlsx')

def compute_expect_Det(CR_all,period_all,k_num,threshold):
    path='/Users/baotong/Desktop/CDFS/'
    info = pd.read_excel(path + 'source_infomation.xlsx')
    for k in k_num:
        DR = read_one_ep(CR_all, period_all, k, threshold,label='REJ1034')
        # DR = read_one_ep(CR_all, period_all, k, threshold, label='A002_R15')
        DR=DR/1000.
        # DR_meanP = DR.mean(axis=1) ##取了三个周期的平均
        DR_meanP=DR[:,0]

        # print(DR_meanP)
        # label='CR_ep'+str(k)
        label='CR'
        CR=info[label]
        CR=np.array(CR)
        CR_use=CR[np.where((CR>=CR_all[0])&(CR<=CR_all[-1]))]
        CR_use=np.sort(CR_use)
        ##在simulation范围内的源的expect探测数目##
        print('Number of used source(low CR)=',len(CR_use))
        ##插个值##
        x=CR_all;y=DR_meanP
        kind=["nearest","zero","slinear","quadratic","cubic"]
        f = interpolate.interp1d(x, y, kind='cubic')
        DR_use=f(CR_use)

        CR_high_use=CR[np.where(CR>CR_all[-1])]
        CR_high_use=np.sort(CR_high_use)
        print('Number of used source(beyond CR)=',len(CR_high_use))
        ##在simulation范围之上的源的expect探测数目##
        z1 = np.polyfit(x, y, 1)
        p1 = np.poly1d(z1)
        DR_high_use=p1(CR_high_use)
        DR_high_use[np.where(DR_high_use>1)]=1

        CR_use=np.concatenate((CR_use,CR_high_use))
        DR_use=np.concatenate((DR_use,DR_high_use))
        print(DR_use)
        print(np.sum(DR_use))

def plot_all_ep(CR_all,period_all,k_num,threshold):
    figlabel = [[0, 0], [0, 1], [1, 0], [1, 1]]
    fig, axes = plt.subplots(2, 2,figsize=(15,10))
    for i in range(len(k_num)):
        k = k_num[i];
        DR= read_one_ep(CR_all,period_all,k, threshold)
        ax_temp = axes[figlabel[i][0], figlabel[i][1]]
        x = CR_all * 1e4
        for j in range(len(period_all)):
            y1=DR[:,j]/10
            ax_temp.plot(x, y1, marker='v', linestyle='-')

        ax_temp.text(5, 25, 'Epoch{0}'.format(k), font1)
        ax_temp.legend(['P=1h','P=1.5h','P=2h'],loc='center left')
        ax_temp.set_ylim(0,30)
        if i<3:ax_temp.legend_.remove()
        # ax_temp.set_title('Epoch {0}: LS detection results'.format(k),font1)
        if (i==2 or i==3):ax_temp.set_xlabel('Count rate ($10^{-4}$ counts$~$$s^{-1}$ )',font1)
        if (i==0 or i==2):ax_temp.set_ylabel('Detection efficiency (%)',font1)
        ax_temp.tick_params(labelsize=16)
    plt.savefig(figurepath + 'DR_thres_{0}.pdf'.format(1 - threshold), bbox_inches='tight', pad_inches=0.0)
    plt.show()

def plot_all_ep_two_model(CR_all,period_all,k_num,threshold):
    figlabel = [[0, 0], [0, 1], [1, 0], [1, 1]]
    fig, axes = plt.subplots(2, 2,figsize=(15,10))
    for i in range(len(k_num)):
        k = k_num[i];
        DR_A= read_one_ep(CR_all,period_all,k, threshold,label='REJ1034')
        DR_B= read_one_ep(CR_all, period_all, k, threshold, label='A002_R15')
        fDR_A=read_one_ep_noqpo(CR_all, k, threshold, label='REJ1034')
        ax_temp = axes[figlabel[i][0], figlabel[i][1]]
        x = CR_all * 1e4
        ax_temp.text(2.8, 60, 'Epoch {0}'.format(k), font1)
        ##==for legend==**
        colorlist=['red','dodgerblue','green','orange']
        if i==3:
            ax_temp.plot([5, 5], [0, 0], linestyle='-', color=colorlist[0])
            ax_temp.plot([5, 5], [0, 0], linestyle='-', color=colorlist[1])
            ax_temp.plot([5, 5], [0, 0], linestyle='-', color=colorlist[2])
            ax_temp.plot([5, 5], [0, 0], linestyle='-', color=colorlist[3])
            ax_temp.plot([5, 5], [5, 5], marker='v', linestyle='-', color='black')
            ax_temp.plot([5, 5], [5, 5], marker='s', linestyle='-.', color='black')
            ax_temp.plot([5, 5], [5, 5], linestyle='dotted', color='grey')
            ax_temp.legend(['P=1h', 'P=1.5h', 'P=2h', 'P=5h','Model A','Model B','fDR'], loc='center left')

        ax_temp.plot(x, fDR_A/10, linestyle='dotted', color='grey')

        # axins = inset_axes(ax_temp, width="40%", height="30%", loc='lower left',
        #                    bbox_to_anchor=(0.3, 0.1, 1, 1),
        #                    bbox_transform=ax_temp.transAxes)
        axins = ax_temp.inset_axes((0.3, 0.67, 0.4, 0.3))
        (xlim0,xlim1,ylim0,ylim1)=(2.8,4.2,0.05,4)
        axins.plot(x, fDR_A/10, linestyle='dotted', color='grey')
        axins.set_yscale('log')
        # axins.set_yticks(np.log10([0.1,0.5,1,2,3]))
        axins.set_yticks([0.1,0.5,1,2])
        axins.set_yticklabels([r'$0.1$',r'$0.5$', r'$1$', r'$2$'])
        # axins.set_yticklabels(['0.1','0.1','0.1','0.1','0.1'])
        axins.set_xlim(xlim0,xlim1)
        axins.set_ylim(ylim0, ylim1)
        if i==3:axins.set_ylim(ylim0, ylim1+1)


        for j in range(len(period_all)):
            y1=DR_A[:,j]/10
            y2=DR_B[:,j]/10
            ax_temp.plot(x, y1+0.1, marker='v', linestyle='-',color=colorlist[j])
            ax_temp.plot(x, y2+0.1, marker='s', linestyle='-.',color=colorlist[j])

            axins.plot(x, y1+0.1, marker='v', linestyle='-',color=colorlist[j])
            axins.plot(x, y2+0.1, marker='s', linestyle='-.',color=colorlist[j])

        ax_temp.set_ylim(-5,65)
        ax_temp.set_xlim(2.5, 10.5)
        # if i < 2: ax_temp.legend_.remove()
        # ax_temp.set_title('Epoch {0}: LS detection results'.format(k),font1)
        if (i == 2 or i == 3): ax_temp.set_xlabel('Count rate ($10^{-4}$ counts$~$$s^{-1}$ )', font1)
        if (i == 0 or i == 2): ax_temp.set_ylabel('Detection efficiency (%)', font1)
        ax_temp.tick_params(labelsize=16)
        axins.tick_params(labelsize=14)

        ###=====画方框========##
        tx0 = xlim0
        tx1 = xlim1
        ty0 = ylim0
        ty1 = ylim1
        sx = [tx0, tx1, tx1, tx0, tx0]
        sy = [ty0, ty0, ty1, ty1, ty0]
        ax_temp.plot(sx, sy, "black")

        # 画两条线
        # xy = (xlim0, ylim0)
        # xy2=(xlim0, ylim0)
        # con = ConnectionPatch(xyA=xy2, xyB=xy, coordsA="data", coordsB="data",
        #                       axesA=axins, axesB=ax_temp)
        # axins.add_artist(con)
        #
        # xy = (xlim1, ylim0)
        # xy2 = (xlim1, ylim0)
        # con = ConnectionPatch(xyA=xy2, xyB=xy, coordsA="data", coordsB="data",
        #                       axesA=axins, axesB=ax_temp)
        # axins.add_artist(con)
    ##====================================##
    plt.savefig(figurepath + 'DR_thres_{0}.pdf'.format(1 - threshold), bbox_inches='tight', pad_inches=0.0)
    plt.show()

if __name__=='__main__':
    # run_test()
    # read_source_info()
    CR_all=np.array([3e-4,4e-4,5e-4,6e-4,7e-4,8e-4,9e-4,1e-3])
    # period_all=np.array([1800,3600,7200])
    period_all = np.array([3600,5400,7200,18000])
    # fDR = read_one_ep_noqpo(CR_all, ep=2, threshold=0.0027, label='REJ1034_R5')
    # print(fDR)
    compute_expect_Det(CR_all, [18000.], k_num=[1,2,3,4], threshold=0.0027)
    # plot_fDR_confirm([5e-4,6e-4,7e-4,8e-4,9e-4,1e-3],ep=4)
    # plot_all_ep(CR_all, period_all, k_num=[1,2,3,4], threshold=0.0027)
    # plot_all_ep_two_model(CR_all, period_all, k_num=[1,2,3,4], threshold=0.0027)
    # plot_sim_one_ep(ep=2, CR_all=np.array([8e-4]),period_all=np.array([18000]),threshold=0.0027,label='A002_R15')
    # plot_sim_one_ep(ep=2,CR_all=np.array([3e-4,4e-4,5e-4,6e-4,7e-4,8e-4,9e-4,1e-3]),period_all=np.array([1800]),threshold=0.0027,label='REJ1034')
    # plot_sim_one_ep(2,threshold=0.0027)
    # plot_sim_one_ep(4, CR_all=CR_all,period_all=period_all,threshold=0.0027)
    # plot_sim_one_ep(4, threshold=0.0027)

