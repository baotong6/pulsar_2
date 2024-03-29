#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import string
from datetime import datetime
from scipy.interpolate import lagrange
from scipy import interpolate
from scipy.optimize import curve_fit
import pandas as pd
from astropy.stats import poisson_conf_interval
import scipy
import hawkeye as hawk
import NGC104.timing_comb as ldata
gcname = ['Tuc', 'terzan5', 'M28', 'omg', 'NGC6397', 'NGC6752', 'NGC6266']

def read_GCsrc(gcname,simid=None):
    path = '/Users/baotong/Desktop/period_' + gcname + '/'+'M31ACIS_txt/'

    if simid:
        id=simid
    else:
        src_info = np.loadtxt(path + 'src_info.txt')
        id = src_info[:, 0]
        counts_all = src_info[:, 3]
        id = id[counts_all > 50];
        counts_all = counts_all[counts_all > 50]
    for k in range(len(id)):
        (src_evt_use, epoch_info_use) = ldata.load_data(dataname=int(id[k]), ecf=90, ifpath=path+'txt_all_obs_p90/',
                                                  ifobsID=[14197])
        obsID = epoch_info_use[:, 2]
        expT = epoch_info_use[:, 3]
        cts = []
        for j in range(len(obsID)):
            cts.append(len(np.where(src_evt_use[:, 2] == obsID[j])[0]))
        cts = np.array(cts)
        CR = cts / expT
        simN=100
        for item in range(simN):
            simt,obsidout=hawk.get_epoch_time_series(cts_rate=CR,period=None,amp=0, model='const',epoch_info=epoch_info_use)
            energy=np.zeros(len(simt))+2000. ##useless
            simevt=np.column_stack((simt,energy,obsidout))
            if os.path.exists(path+'sim_const/'):
                np.savetxt(path+f'sim_const/sim_{int(id[k])}_f{item}.txt',simevt,fmt="%.7f  %5.3f  %d")
            else:
                os.system(f'mkdir {path}sim_const/')
def GLres():
    for i in range(len(gcname)):
        path = '/Users/baotong/Desktop/period_' + gcname[i] + '/'
        src_info = np.loadtxt(path + 'src_info.txt')
        id=src_info[:,0]
        counts_all = src_info[:, 3]
        id=id[counts_all>50];counts_all=counts_all[counts_all>50]
        w_range = 2 * np.pi * np.arange(1 / 50000., 1 / 10000, 1e-8)
        for j in range(len(id)):
            print(int(id[j]))
            srcevt=np.loadtxt(path+f'sim_const/sim_{int(id[j])}.txt')
            epoch_info=np.loadtxt(path+f'txt_all_obs_p90/epoch_src_{int(id[j])}.txt')
            res=hawk.GL_algorithm_single.write_result(srcevt,epoch_info,w_range,
                                                      dataname=str(int(id[j])),if_filter=False,
                                                      pathout='/Users/baotong/Desktop/period_' + gcname[i] + '/sim_GL/')
def read_GLres(srcid):
    path='/Users/baotong/Desktop/period_M31XRB/M31HRC_txt/'
    for j in range(len(srcid)):
        prob_GL = []
        for i in range(100):
            if os.path.isfile(path + 'sim_GL/' + 'result_3h_{0}_f{1}.txt'.format(str(srcid[j]), i)):
                res = np.loadtxt(path + 'sim_GL/' + 'result_3h_{0}_f{1}.txt'.format(str(srcid[j]), i))
            elif os.path.isfile(path + 'sim_GL/' + 'result_10h_{0}_f{1}.txt'.format(str(srcid[j]), i)):
                res = np.loadtxt(path + 'sim_GL/' + 'result_10h_{0}_f{1}.txt'.format(str(srcid[j]), i))
            else:
                print('None')
                continue
            prob_GL.append(res[2])
        prob_GL = np.array(prob_GL)
    # print(prob_GL)
    fD = len(prob_GL[prob_GL > 0.95])
    print(fD)
    print(np.where(prob_GL>0.9)[0])
def plot_GLres(srcid,pathin=None,labellist=None):
    for j in range(len(srcid)):
        prob_GL = []
        for i in range(100):
            if os.path.isfile(pathin + 'sim_GL/' + 'result_3h_{0}_f{1}.txt'.format(str(srcid[j]), i)):
                res = np.loadtxt(pathin + 'sim_GL/' + 'result_3h_{0}_f{1}.txt'.format(str(srcid[j]), i))
            elif os.path.isfile(pathin + 'sim_GL/' + 'result_10h_{0}_f{1}.txt'.format(str(srcid[j]), i)):
                res = np.loadtxt(pathin + 'sim_GL/' + 'result_10h_{0}_f{1}.txt'.format(str(srcid[j]), i))
            else:
                print('None')
                continue
            prob_GL.append(res[2])
        prob_GL = np.array(prob_GL)
        bins=np.linspace(0,1,21)
        if labellist:
            plt.hist(prob_GL, bins=18, lw=2, histtype='step', label=f'{labellist[j]}')
        else:
            plt.hist(prob_GL,bins=18,lw=2,histtype='step',label=f'{srcid[j]}')
    # plt.xlabel(r'$P_{\rm GL}$', hawk.font1)
    # plt.ylabel('Number of trials', hawk.font1)
    # plt.plot([0.95, 0.95], [0, 100], '--', color='grey')
    # plt.tick_params(labelsize=16)
    # plt.semilogy()
    # plt.legend(loc='upper center', ncol=5, handletextpad=0.5, columnspacing=0.5, facecolor=None,
    #                      edgecolor='grey', handlelength=0.5)
    # plt.show()
if __name__=='__main__':
    # read_GCsrc('M31XRB',simid=[2])
    # GLres()
    # read_GLres([6])
    path1='/Users/baotong/Desktop/period_M31XRB/M31HRC_txt/'
    path2='/Users/baotong/Desktop/period_M31XRB/M31ACIS_txt/'
    plot_GLres([38,37,6,83,167],path1,labellist=['1','2','4','5','7'])
    plot_GLres([2,199,21], path2,labellist=['3','5','8'])
    plt.xlabel(r'$P_{\rm GL}$', hawk.font1)
    plt.ylabel('Number of trials', hawk.font1)
    plt.plot([0.95, 0.95], [0, 100], '--', color='grey')
    plt.tick_params(labelsize=16)
    plt.semilogy()
    # plt.semilogx()
    plt.legend(loc='upper center', ncol=5, handletextpad=0.5, columnspacing=0.5, facecolor=None,
                         edgecolor='grey', handlelength=1.5)
    plt.savefig('/Users/baotong/Desktop/aas/M31XRB/fig/'+'hist_sim.pdf',bbox_inches='tight', pad_inches=0.05)
    plt.show()
    # M31 199=18