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
import linecache
from astropy.timeseries import LombScargle

font1 = {'family': 'Normal',
         'weight': 'normal',
         'size': 18, }
def get_line_context(file_path, line_number):
    return linecache.getline(file_path, line_number).strip()



def plot_lc_ap(mode):
    #path_out='/Users/baotong/Desktop/aas/MN_pCV/figure/'+mode+'/'
    def pfold(time,P,flux,):
        def trans(t, p_test, shift=0.5):
            ti = t
            v = 1.0 / p_test
            turns = v * ti
            turns += shift
            # 初始相位
            for i in range(len(turns)):
                turns[i] = turns[i] - int(turns[i])
            return turns
        turns=trans(time,P)
        plt.figure(1,(8,6))
        plt.xlabel('phase')
        plt.ylabel('photon flux')
        #print(useobs[np.where((turns<0.5)&(turns>0.2))])
        plt.errorbar(turns,flux[0],yerr=[flux[0]-flux[1],flux[2]-flux[0]],
                     fmt = 'o', capsize = 3, elinewidth = 1, color='red',ecolor = 'red')
        plt.errorbar(turns+1., flux[0], yerr = [flux[0] - flux[1], flux[2] - flux[0]],
                     fmt = 'o', capsize = 3, elinewidth = 1, color = 'red', ecolor = 'red')
        #plt.savefig('/Volumes/pulsar/WR/1671/pfold_2_4.eps')
        plt.show()
    if mode=='ND':
        src_ID=data.ID_ND
        obs_ID=np.loadtxt('ACIS-I_epoch.txt')[:,2]
        obs_time=(np.loadtxt('ACIS-I_epoch.txt')[:,0]+np.loadtxt('ACIS-I_epoch.txt')[:,1])/2
        time=obs_time/86400+2449352.5-2400000.5
        os.chdir(path_ND)
    elif mode=='LW':
        src_ID=data.ID_LW
        #src_ID=['202']
        obs_ID=np.loadtxt('LW_epoch.txt')[:,2]
        obs_time=(np.loadtxt('LW_epoch.txt')[:,0]+np.loadtxt('LW_epoch.txt')[:,1])/2
        os.chdir(path_LW)
        time = obs_time / 86400 + 2449352.5 - 2400000.5
    elif mode == '47Tuc':
        src_ID = ['245', '453', '407', '182', '304', '224', '206', '223', '211',
                  '462', '402', '294', '485', '549', '258', '5', '508', '520',
                  '292', '350', '283', '314', '284', '62', '367', '345', '364']  # 47Tuc
        # src_ID = ['245']
        path_Tuc = '/Volumes/pulsar/47Tuc/merge_data/spectra/aprates/'
        EPOCH = np.loadtxt(path_Tuc + '47Tuc_epoch.txt')
        obs_ID = EPOCH[:, 2]
        obs_time = (EPOCH[:, 0] + EPOCH[:, 1]) / 2
        os.chdir(path_Tuc)
        time = obs_time / 86400 + 2449352.5 - 2400000.5

    elif mode=='terzan5':
        src_ID = ['98', '283', '282', '8', '377', '299', '33', '91', '74', '60', '117', '84', '37', '85',
                  '128', '89', '78', '372', '374', '66', '43', '83', '375', '14', '292', '352', '55', '294']  # terzan5
        #src_ID=['71']
        path_Tuc = '/Volumes/pulsar/terzan5/merge_data/spectra/aprates/'
        EPOCH = np.loadtxt(path_Tuc + 'terzan5_epoch.txt')
        obs_ID = EPOCH[:, 2]
        obs_time = (EPOCH[:, 0] + EPOCH[:, 1]) / 2
        os.chdir(path_Tuc)
        time = obs_time / 86400 + 2449352.5 - 2400000.5
    elif mode == 'M28':
        src_ID = ['22','92','131','321','320','38','366','134','91','13','85','82','251','326']  #
        #src_ID = ['245']
        path_Tuc = '/Volumes/pulsar/M28/merge_data/spectra/aprates/'
        EPOCH = np.loadtxt(path_Tuc + 'M28_epoch.txt')
        obs_ID = EPOCH[:, 2]
        obs_time = (EPOCH[:, 0] + EPOCH[:, 1]) / 2
        os.chdir(path_Tuc)
        time = obs_time / 86400 + 2449352.5 - 2400000.5
    elif mode == 'omg_cen':
        src_ID = ['69','52','36','87','25']  #
        path_Tuc = '/Volumes/pulsar/omega_cen/merge_data/spectra/aprates/'
        EPOCH = np.loadtxt(path_Tuc + 'omg_cen_epoch.txt')
        obs_ID = EPOCH[:, 2]
        obs_time = (EPOCH[:, 0] + EPOCH[:, 1]) / 2
        os.chdir(path_Tuc)
        time = obs_time / 86400 + 2449352.5 - 2400000.5
    for k in range(len(src_ID)):
        ID = int(src_ID[k])
        ID_note = int(src_ID[k])
        if str(int(ID))[-3:] == '001' or str(int(ID))[-3:] == '002':
            ID = str(ID)[:-3]
        use_obs = []
        use_obs_time = []
        cts_rate_low = []
        cts_rate_high = []
        cts_rate = []
        for i in range(len(obs_ID)):
            outname = str(ID) + '_' + str(int(obs_ID[i])) + '_out.par'
            if os.path.exists(outname):
                use_obs.append(obs_ID[i])
                use_obs_time.append(time[i])
                a = get_line_context(outname, 15)[18:-3]
                a_low = get_line_context(outname, 16)[25:-3]
                a_high = get_line_context(outname, 17)[25:-3]
                if a == 'INDEF':
                    a = 0
                if a_low == 'INDEF':
                    a_low = 0
                if a_high == 0:
                    continue

                cts_rate.append(float(a))
                cts_rate_low.append(float(a_low))
                cts_rate_high.append(float(a_high))

            else:
                print('continue')
                continue

        def get_VI(cts_rate_m, cts_rate_low, cts_rate_high):
            cts_rate = cts_rate_m
            low = cts_rate_low
            high = cts_rate_high
            i = 0
            if len(cts_rate) != len(low) != len(high):
                print('error')
                return None
            else:
                while i < len(low):
                    if low[i] == 0:
                        low = np.delete(low, i)
                        high = np.delete(high, i)
                        cts_rate = np.delete(cts_rate, i)
                    else:
                        i += 1
            if len(cts_rate)==0:
                return 0
            VI = np.max(cts_rate) / np.min(cts_rate)
            return VI

        VI = get_VI(cts_rate, cts_rate_low, cts_rate_high)
        cts_rate=np.array(cts_rate)*1e7
        cts_rate_low=np.array(cts_rate_low)*1e7
        cts_rate_high=np.array(cts_rate_high)*1e7

        path_aprates='/Users/baotong/Desktop/period_terzan5/aprates/'
        out_info = np.column_stack((cts_rate, cts_rate_low, cts_rate_high, use_obs))
        np.savetxt(path_aprates + str(ID) + '.txt', out_info, fmt='%10.3f %10.3f %10.3f %10d')

        # plt.figure(1,(10,7.5))
        # #fig, (ax, ax2) = plt.subplots(1, 2,figsize=(10,7.5),sharey = True)
        # plt.semilogy()
        # plt.title('{0} #{1}, VI={2}'.format(mode, k + 1, round(VI, 2)), fontsize = 20)
        # #ax.set_title('{0} #{1}, VI={2}'.format(mode, k + 1, round(VI, 2)), fontsize = 20)
        # #ax.set_title('{0} {1}, VI={2}'.format(mode, 'F1', round(VI, 2)), fontsize = 20)
        # plt.xlabel('MJD', font1)
        # plt.ylabel('photon flux ( '+r'$\rm 10^{-7} ph \,s^{-1} \, cm^{-2}$'+' )', font1)
        #
        # plt.tick_params(labelsize = 18)

        # plt.figure(1,(10,7.5))
        fig, (ax, ax2) = plt.subplots(1, 2, figsize=(10, 7.5), sharey=True)
        plt.semilogy()
        #plt.title('{0} #{1}, VI={2}'.format(mode, k + 1, round(VI, 2)), fontsize=20)
        plt.title('{0}, VI={1}'.format(src_ID[k], round(VI, 2)), fontsize=20)
        # ax.set_title('{0} #{1}, VI={2}'.format(mode, k + 1, round(VI, 2)), fontsize = 20)
        # ax.set_title('{0} {1}, VI={2}'.format(mode, 'F1', round(VI, 2)), fontsize = 20)
        # plt.ylabel('cts rate  '+r'$(erg\cdot s^{-1} \cdot cm^{-2})$', fontsize = 15)
        # plt.xlabel('MJD', fontsize = 15)

        ax2.set_xlabel('MJD', font1)
        ax.set_ylabel('photon flux ( ' + r'$\rm 10^{-7} ph \,s^{-1} \, cm^{-2}$' + ' )', font1)
        ax.tick_params(labelsize=15)
        ax2.tick_params(labelsize=15)
        # plt.ylim(5e-8,3e-5)
        # if np.min(cts_rate_low) > 1e-10:
        #     plt.ylim(ymin = np.min(cts_rate_low) * 0.8, ymax = np.max(cts_rate_high) * 1.2)
        # else:
        #     plt.ylim(ymin = 6e-8, ymax = np.max(cts_rate_high) * 1.2)
        #plt.xlim(np.min(use_obs_time) - 100, np.max(use_obs_time) + 100)
        for i in range(len(use_obs_time)):
            if cts_rate_low[i] == 0:
                ax.plot([use_obs_time[i] - 5, use_obs_time[i] + 5], [cts_rate_high[i], cts_rate_high[i]], color='red')
                ax.annotate("", xy=(use_obs_time[i], cts_rate_high[i] * 0.8),
                            xytext=(use_obs_time[i], cts_rate_high[i]), color="red",
                            weight="bold",
                            arrowprops=dict(arrowstyle="->", connectionstyle="arc3", color="red"))
                ax2.plot([use_obs_time[i] - 4, use_obs_time[i] + 4], [cts_rate_high[i], cts_rate_high[i]],
                         color='red')
                ax2.annotate("", xy=(use_obs_time[i], cts_rate_high[i] * 0.8),
                             xytext=(use_obs_time[i], cts_rate_high[i]), color="red",
                             weight="bold",
                             arrowprops=dict(arrowstyle="->", connectionstyle="arc3", color="red"))
                # plt.annotate("", xy = (use_obs_time[i], cts_rate_high[i] * 0.5),
                #              xytext = (use_obs_time[i], cts_rate_high[i]), color = "red",
                #              weight = "bold",
                #              arrowprops = dict(arrowstyle = "->", connectionstyle = "arc3", color = "red"))

            else:
                ax.errorbar(use_obs_time[i], cts_rate[i],
                            yerr=[[cts_rate[i] - cts_rate_low[i]], [cts_rate_high[i] - cts_rate[i]]],
                            fmt='o', capsize=3, elinewidth=1, ecolor='red')
                ax2.errorbar(use_obs_time[i], cts_rate[i],
                             yerr=[[cts_rate[i] - cts_rate_low[i]], [cts_rate_high[i] - cts_rate[i]]],
                             fmt='o', capsize=3, elinewidth=1, ecolor='red')
                # plt.errorbar(use_obs_time[i], cts_rate[i],
                #         yerr = [[cts_rate[i] - cts_rate_low[i]], [cts_rate_high[i] - cts_rate[i]]],
                #         fmt = 'o', capsize = 3, elinewidth = 1, ecolor = 'red')

        ###########################################################################################
        # For LW
        # ax.set_xlim(52100, 53150)  # most of the data
        # ax2.set_xlim(53200, 53230)  # outliers only

        # For ND
        # ax.set_xlim(50250, 50300)  # most of the data
        # ax2.set_xlim(55000, 56700)  # outliers only

        # For 47Tuc
        ax.set_xlim(50000, 51100)  # most of the data
        ax2.set_xlim(55400, 55630)  # outliers only

        # # For terzan5
        # ax.set_xlim(51000, 52000)  # most of the data
        # ax2.set_xlim(53500, 56500)  # outliers only
        #
        # # For M28
        # ax.set_xlim(50990, 51080)  # most of the data
        # ax2.set_xlim(53100, 56000)  # outliers only
        #
        # # For omg_cen
        # ax.set_xlim(50100, 50110)  # most of the data
        # ax2.set_xlim(54550, 54590)  # outliers only

        # hide the spines between ax and ax2
        ax.spines['right'].set_visible(False)
        ax2.spines['left'].set_visible(False)
        ax.yaxis.tick_left()
        # ax.tick_params(labeltop = 'off')# don't put tick labels at the top
        ax2.yaxis.tick_right()

        # Make the spacing between the two axes a bit smaller
        plt.subplots_adjust(wspace=0.15)

        # This looks pretty good, and was fairly painless, but you can get that
        # cut-out diagonal lines look with just a bit more work. The important
        # thing to know here is that in axes coordinates, which are always
        # between 0-1, spine endpoints are at these locations (0,0), (0,1),
        # (1,0), and (1,1). Thus, we just need to put the diagonals in the
        # appropriate corners of each of our axes, and so long as we use the
        # right transform and disable clipping.

        d = .015  # how big to make the diagonal lines in axes coordinates
        # arguments to pass plot, just so we don't keep repeating them
        kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
        ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-left diagonal
        ax.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal

        kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes

        ax2.plot((-d, d), (-d, +d), **kwargs)  # top-right diagonal
        ax2.plot((-d, d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal
        ###########################################################################################
        path_out='/Users/baotong/Desktop/aas/pCV_GC/figure/47Tuc/'
        #plt.savefig(path_out + str(ID_note) + '_lc.eps')
        #plt.show()
        plt.close()
plot_lc_ap('terzan5')

