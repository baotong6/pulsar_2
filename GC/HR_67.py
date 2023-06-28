#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
from scipy.interpolate import lagrange
from scipy import optimize as op
from scipy.optimize import curve_fit
import pandas as pd
from astropy.stats import poisson_conf_interval
import hawkeye as hawk

pos_all={'Tuc':[6.0236250, -72.0812833, 3.17 * 60, 3.17 / 8.8 * 60],
         'terzan5':[267.0202083, -24.7790556, 0.72 * 60, 0.72 / 3.4 * 60],
         'M28':[276.1363750, -24.8702972, 1.97 * 60, 1.97 / 8.2 * 60],
         'omg':[201.69700, -47.47947, 5 * 60, 5 / 2.1 * 60],
         'NGC6397':[265.17539, -53.67433, 2.9 * 60, 2.9 / 58 * 60],
         'NGC6752':[287.71713, -59.98455, 1.91* 60, 1.91 / 11.24 * 60],
         'NGC6266':[255.303333,-30.113722,0.92* 60,0.92/4.2*60],
         'M30':[325.092167,-23.179861,1.03* 60,1.03/17.2*60]}

def f(x,a,b):
    logS=a*x+b  #S~V^a
    return logS
def spectrafit(x,y,error):
    popt, pcov = op.curve_fit(f, np.log(x), np.log(y),absolute_sigma=True,sigma=np.log(error))
    # popt, pcov = op.curve_fit(f, np.log(x), np.log(y))
    perr = np.sqrt(np.diag(pcov))
    logydata=f(np.log(x),popt[0],popt[1])
    ydata=np.exp(logydata)
    return (popt,perr)

def plot_HR67_SH(save=0,show=1):
    gclist=['M28','NGC6397','terzan5','omg','NGC6752','NGC6266']
    color_list = ['grey', 'g', 'b', 'k', 'orange', 'purple']
    linesrc=[84,306,193]
    HR=[];HRl=[];HRu=[];F=[];Fl=[];Fu=[]
    for i in range(len(gclist)):
        gc=gclist[i]
        file = pd.read_csv(f'/Users/baotong/Desktop/period_{gc}/HR/S26_H67_evt/all_HR_sup.csv')
        pars=np.loadtxt(f'/Users/baotong/Desktop/period_{gc}/HR/S26_H67_evt/all_inputpars.txt')
        srcinfo=np.loadtxt(f'/Users/baotong/Desktop/period_{gc}/src_info.txt')
        dist=srcinfo[:,6]
        S=pars[:,1]-pars[:,3]/pars[:,5];H=pars[:,2]-pars[:,4]/pars[:,6]
        seq=pars[:,0]
        S[np.isnan(S)]=0;H[np.isnan(H)]=0
        HR=file['HR'];HRl=file['HRl'];HRu=file['HRu']


        SH=H+S
        SH_1sigma = poisson_conf_interval(SH, interval='frequentist-confidence').T
        SHl = SH_1sigma[:,0];SHu = SH_1sigma[:,1]
        # plt.scatter(HR,SH,label=gc)
        ##==plot HR-S+H diagram==##
        plt.errorbar(x=HR,y=SH,yerr=[SH-SHl,SHu-SH],xerr=[HR-HRl,HRu-HR],fmt = '.', capsize = 3, elinewidth = 1, label=gc,color=color_list[i])
        if i<=2:
            lineindex=np.where(pars[:,0]==linesrc[i])[0]
            plt.text(x=HR[lineindex],y=SH[lineindex]*1.1,s=str(int(pars[:,0][lineindex][0])),fontdict=hawk.font1,color=color_list[i])
        outlierindex=np.where((HR>-0.85)&(SH>100))
        for outindex in outlierindex[0]:
            print(str(int(pars[:, 0][outindex])))
            plt.text(x=HR[outindex],y=SH[outindex]*1.1,s=str(int(pars[:,0][outindex])),fontdict=hawk.font1,color=color_list[i])
    plt.legend()
    plt.ylim(bottom=100)
    plt.xlim(left=-1.05, right=-0.5)
    plt.tick_params(labelsize=18)
    plt.semilogy()
    plt.xlabel('HR', fontdict=hawk.font1)
    plt.ylabel('S+H', fontdict=hawk.font1)
    if save:plt.savefig('HR.pdf', bbox_inches='tight', pad_inches=0.05)
    if show:plt.show()

def plot_HR67_dist(save=0,show=1):
    gclist=['M28','NGC6397','terzan5','omg','NGC6752','NGC6266','Tuc']
    color_list = ['grey', 'g', 'b', 'k', 'orange', 'purple','c']
    linesrc=[84,306,193]
    HR=[];HRl=[];HRu=[];F=[];Fl=[];Fu=[]
    for i in range(len(gclist)):
        gc=gclist[i]
        file = pd.read_csv(f'/Users/baotong/Desktop/period_{gc}/HR/S26_H67_evt/all_HR_sup.csv')
        pars=np.loadtxt(f'/Users/baotong/Desktop/period_{gc}/HR/S26_H67_evt/all_inputpars.txt')
        srcinfo=np.loadtxt(f'/Users/baotong/Desktop/period_{gc}/src_info.txt')
        if gc=='Tuc':dist=srcinfo[:,-1]/pos_all[gc][-2]
        else:dist=srcinfo[:,6]/pos_all[gc][-2]
        S=pars[:,1]-pars[:,3]/pars[:,5];H=pars[:,2]-pars[:,4]/pars[:,6]
        seq=pars[:,0]
        seq=seq.astype('int')
        S[np.isnan(S)]=0;H[np.isnan(H)]=0
        HR=file['HR'];HRl=file['HRl'];HRu=file['HRu']
        dist =dist[seq-1]
        SH=H+S
        SH_1sigma = poisson_conf_interval(SH, interval='frequentist-confidence').T
        SHl = SH_1sigma[:,0];SHu = SH_1sigma[:,1]
        ##==plot HR-S+H diagram==##
        goodindex=np.where(SH>100)[0]
        print(len(goodindex))
        plt.errorbar(x=HR[goodindex],y=dist[goodindex],xerr=[HR[goodindex]-HRl[goodindex],HRu[goodindex]-HR[goodindex]],
                     fmt = '.', capsize = 3, elinewidth = 1, label=gc,color=color_list[i])
        # if i<=2:
        #     lineindex=np.where(pars[:,0]==linesrc[i])[0]
        #     plt.text(x=HR[lineindex],y=dist[lineindex]*1.1,s=str(int(pars[:,0][lineindex][0])),fontdict=hawk.font1,color=color_list[i])
        # outlierindex=np.where((HR>-0.85)&(SH>100))
        # for outindex in outlierindex[0]:
        #     # print(str(int(pars[:, 0][outindex])))
        #     plt.text(x=HR[outindex],y=dist[outindex]*1.1,s=str(int(pars[:,0][outindex])),fontdict=hawk.font1,color=color_list[i])
    plt.legend()
    # plt.xlim(left=-1.05, right=-0.5)
    plt.tick_params(labelsize=18)
    plt.semilogy()
    plt.xlabel('HR', fontdict=hawk.font1)
    plt.ylabel('dist', fontdict=hawk.font1)
    if save:plt.savefig('HR-dist.pdf', bbox_inches='tight', pad_inches=0.05)
    if show:plt.show()

def plot_HR67_hist(save=0,show=1):
    # gclist=['M28','NGC6397','terzan5','omg','NGC6752','NGC6266','Tuc']
    # color_list = ['grey', 'g', 'b', 'k', 'orange', 'purple','c']
    gclist=['terzan5']
    color_list = ['b']
    HR=[];HRl=[];HRu=[];F=[];Fl=[];Fu=[]
    for i in range(len(gclist)):
        gc=gclist[i]
        file = pd.read_csv(f'/Users/baotong/Desktop/period_{gc}/HR/S26_H67_evt/all_HR_sup.csv')
        pars=np.loadtxt(f'/Users/baotong/Desktop/period_{gc}/HR/S26_H67_evt/all_inputpars.txt')
        srcinfo=np.loadtxt(f'/Users/baotong/Desktop/period_{gc}/src_info.txt')
        if gc=='Tuc':dist=srcinfo[:,-1]/pos_all[gc][-2]
        else:dist=srcinfo[:,6]/pos_all[gc][-2]
        ## /half-light radius
        S=pars[:,1]-pars[:,3]/pars[:,5];H=pars[:,2]-pars[:,4]/pars[:,6]
        seq=pars[:,0]
        seq=seq.astype('int')
        S[np.isnan(S)]=0;H[np.isnan(H)]=0
        HR=file['HR'];HRl=file['HRl'];HRu=file['HRu']
        dist =dist[seq-1]
        SH=H+S
        SH_1sigma = poisson_conf_interval(SH, interval='frequentist-confidence').T
        SHl = SH_1sigma[:,0];SHu = SH_1sigma[:,1]
        ##==plot HR-histogram==##
        goodindex=np.where(SH>100)[0]
        outindex=np.intersect1d(goodindex,np.where(HR>-0.85)[0])
        print(seq[outindex])
        plt.figure(1,(9,6))
        plt.hist(HR[goodindex],bins=20,histtype='step',label=gc,linewidth=3)
        plt.xticks(size=18)
        plt.yticks(size=18)
        plt.xlabel('HR', hawk.font1)
        plt.ylabel('Number of bin',hawk.font1)
        # plt.scatter(HR[goodindex],dist[goodindex],marker='+',s=50)
        plt.semilogy()
        plt.show()

def plot_HR67_profile(save=0,show=1):
    gclist=['terzan5']
    color_list = ['b']
    for i in range(len(gclist)):
        gc=gclist[i]
        file = pd.read_csv(f'/Users/baotong/Desktop/period_{gc}/HR/S26_H67_evt/all_HR_sup.csv')
        pars=np.loadtxt(f'/Users/baotong/Desktop/period_{gc}/HR/S26_H67_evt/all_inputpars.txt')
        srcinfo=np.loadtxt(f'/Users/baotong/Desktop/period_{gc}/src_info.txt')
        dist=srcinfo[:,6]/pos_all[gc][-2]      ## /half-light radius
        S=pars[:,1]-pars[:,3]/pars[:,5];H=pars[:,2]-pars[:,4]/pars[:,6]
        seq=pars[:,0]
        seq=seq.astype('int')
        S[np.isnan(S)]=0;H[np.isnan(H)]=0
        HR=file['HR'];HRl=file['HRl'];HRu=file['HRu']
        dist =dist[seq-1]
        SH=H+S
        SH_1sigma = poisson_conf_interval(SH, interval='frequentist-confidence').T
        SHl = SH_1sigma[:,0];SHu = SH_1sigma[:,1]
        ##==plot HR-histogram==##
        goodindex=np.where(SH>100)[0]
        outindex=np.intersect1d(goodindex,np.where(HR>-0.85)[0])
        print(seq[outindex])
        print(dist[outindex])
        bins_dist=np.logspace(np.log10(0.03),np.log10(6),5)
        bins2=np.logspace(np.log10(0.03),np.log10(6),10)
        hist1 = plt.hist(dist[outindex], bins_dist)
        hist2 = plt.hist(dist[goodindex], bins2)
        plt.close()
        y = hist1[0];y2=hist2[0]
        x = [(bins_dist[i] + bins_dist[i + 1]) / 2 for i in range(len(bins_dist) - 1)]
        area = [(bins_dist[i + 1] ** 2 - bins_dist[i] ** 2) * np.pi for i in range(len(bins_dist) - 1)]

        x2 = [(bins2[i] + bins2[i + 1]) / 2 for i in range(len(bins2) - 1)]
        area2 = [(bins2[i + 1] ** 2 - bins2[i] ** 2) * np.pi for i in range(len(bins2) - 1)]

        y_err = np.array(poisson_conf_interval(y, interval='frequentist-confidence'))
        y_err[0] = y - y_err[0];y_err[1] = y_err[1] - y
        y = y / area;y_err = y_err / area

        y2_err = np.array(poisson_conf_interval(y2, interval='frequentist-confidence'))
        y2_err[0] = y2 - y2_err[0];y2_err[1] = y2_err[1] - y2
        y2 = y2 / area2;y2_err = y2_err / area2

        xerr = [(bins_dist[i + 1] - bins_dist[i]) / 2 for i in range(len(bins_dist) - 1)]
        x2err = [(bins2[i + 1] - bins2[i]) / 2 for i in range(len(bins2) - 1)]
        (popt1, perr1) = spectrafit(x, y, y_err[0])
        (popt2, perr2) = spectrafit(x2, y2, y2_err[0])
        print('popt1=', popt1)
        print('perr1=', perr1)
        print('popt2=', popt2)
        print('perr2=', perr2)
        plt.figure(1,(9,6))
        plt.errorbar(x, y*10, xerr=xerr, yerr=y_err*10, fmt='ks', capsize=2, elinewidth=2, ecolor='k', color='k',
                     markersize=8, label=r'$\rm Metal-rich~CVs \times 10$')

        plt.errorbar(x2, y2, xerr=x2err, yerr=y2_err, fmt='ro', capsize=1, elinewidth=1, ecolor='r', color='r',
                     markersize=4, label='Bright group')
        plt.plot(x2, 10*np.exp(f(np.log(x2), popt1[0], popt1[1])), '-.', color='k')
        plt.plot(x2, np.exp(f(np.log(x2), popt2[0], popt2[1])), '-.', color='r')
        plt.loglog()
        plt.legend()
        plt.xticks(size=18)
        plt.yticks(size=18)
        plt.xlabel(r'$\rm R/r_c$ ', hawk.font1)
        plt.ylabel(r'$\rm Number~of~source~per~arcmin^2$', hawk.font1)
        plt.show()
if __name__=='__main__':
    # plot_HR67_SH()
    # plot_HR67_dist()
    plot_HR67_hist()
    # plot_HR67_profile()