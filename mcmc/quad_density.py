#!/bin/bash
# -*- coding: utf-8 -*-
# written by Tong
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
from matplotlib import ticker
from astropy.io import fits
import astropy.coordinates as coord
from astropy.coordinates import SkyCoord
from scipy.integrate import quad
import astropy.units as u
import hawkeye as hawk

c_LW = SkyCoord(ra=267.86375 * u.degree, dec=-29.58475 * u.degree)
c_terzan5= SkyCoord(ra=267.02 * u.degree, dec=-24.779167 * u.degree)
c_terzan5_icrs = SkyCoord('17h48m05s', '-24d46m45s', frame='icrs')
c_M28 = SkyCoord(ra=276.136708 * u.degree, dec=	-24.869778 * u.degree)
c_M28_icrs = SkyCoord('18h24m33s', '-24d52m12s', frame='icrs')
dist_terzan5=c_LW.separation(c_terzan5).degree
dist_M28=c_LW.separation(c_M28).degree
def pnd(d,z):
    if d>220:
        return 300*d**(-10)*np.exp(-z/0.045)
    elif 120<=d<=220:
        return 300*d**(-3.5)*np.exp(-z/0.045)
    elif d<120:
        return 300*d**(-0.1)*np.exp(-z/0.045)
def pho(d):
    z=0.19815
    pbulge=1.09*np.sqrt(d**2+(z/0.6)**2)**(-1.8)*np.exp(-d**2/1.9-(z/0.6)**2/1.9)
    pdisk=2.5*np.exp(-(3/d)**3-d/2.5-np.abs(z/130))
    ptot=pbulge+pdisk+pnd(1000*d,z)
    return ptot*1e9

def func_LW(r):
    h0=8*np.sin(1.4/180*3.14)
    d=np.sqrt((8-r)**2+h0**2)
    A0=15*15*(60*0.03867)**2/1e6
    Mtot=A0*pho(d)*r**2/8**2
    return Mtot

def func_terzan(r):
    h0=8*np.sin(dist_terzan5/180*3.14)
    d=np.sqrt((8-r)**2+h0**2)
    A0=3.14*(60*0.03867)**2/1e6
    Mtot=A0*pho(d)*r**2/8**2
    return Mtot
def func_M28(r):
    h0=8*np.sin(dist_M28/180*3.14)
    d=np.sqrt((8-r)**2+h0**2)
    A0=3.14*(60*0.03867)**2/1e6
    Mtot=A0*pho(d)*r**2/8**2
    return Mtot
# plot_aaa()
def plot_cum_density(save=0,show=1):
    A_LW=15*15
    a = quad(func_LW, 0, 8)
    pho_LW=a[0]/(15*15)*3.14
    b = quad(func_terzan, 0, 8)
    c = quad(func_M28, 0, 8)
    pho_terzan=b[0];pho_M28=c[0]
    r_terzan=np.sqrt(15*15/3.14*pho_LW/pho_terzan)
    r_M28=np.sqrt(15*15/3.14*pho_LW/pho_M28)
    r1 = np.arange(0, r_terzan, 0.1)
    r2 = np.arange(0, r_M28, 0.1)
    # 根据半径计算每个圆的面积
    areas1 = np.pi * r1 ** 2;areas2 = np.pi * r2 ** 2
    normalized_areas1 = areas1 / areas1[-1]*22
    d_areas1 = normalized_areas1[1:] - normalized_areas1[:-1]
    cumulative_areas1 = np.cumsum(d_areas1)
    normalized_areas2 = areas2/ areas2[-1]*22
    d_areas2 = normalized_areas2[1:] - normalized_areas2[:-1]
    cumulative_areas2 = np.cumsum(d_areas2)
    fig, ax = plt.subplots(figsize=(8, 6))
    bins=np.linspace(0,10,1000)
    srcdist_terzan=np.sort(np.array([10.42802,50.94517,8.1156,20.74795,293.964,185.54723,171.48937])/60)
    srcdist_M28=np.sort(np.array([153.76256,491.98014,183.42538,213.151])/60)
    ax.plot(r1[:-1], cumulative_areas1, color='green',label='Bulge/disk contribution for Terzan 5',linestyle='dotted',linewidth=4)
    ax.plot(r2[:-1], cumulative_areas2, color='k',label='Bulge/disk contribution for M 28',linestyle='dotted',linewidth=4)
    # n_terzan, bins_terzan, patches_terzan=ax.hist(srcdist_terzan, bins=bins, histtype='step', lw=2, color='green',
    #                                               cumulative=1, linestyle='solid')
    # n_M28, bins_M28, patches_M28=ax.hist(srcdist_M28, bins=bins, histtype='step', lw=2, color='k',
    #                                      cumulative=1,  linestyle='solid')
    ax.plot(srcdist_terzan,np.arange(1,len(srcdist_terzan)+1,1),'o-',lw=2, color='green',label='Deteced periodic CVs in Terzan 5')
    ax.plot(srcdist_M28,np.arange(1,len(srcdist_M28)+1,1), 'o-',lw=2,color='black',label='Deteced periodic CVs in M 28')
    # 标注数字
    for i, value in enumerate(srcdist_terzan):
        plt.text(value, 1.3*(i + 1), str(i+1),fontsize=15)

    for i, value in enumerate(srcdist_M28):
        plt.text(value, 1.3*(i + 1), str(i+1),fontsize=15)
    # num1 = 0;
    # num2 = 0
    # for i in range(len(bins)-2):
    #     if n_terzan[i+1]>n_terzan[i]:
    #         plt.text(bins[i]-0.001, n_terzan[i]*0.5, num1, ha='center', va='center',fontsize=15,color='green')
    #         num1+=1
    #     if n_M28[i+1]>n_M28[i]:
    #         plt.text(bins[i]-0.001, n_M28[i]*0.5, num2, ha='center', va='center',fontsize=15,color='black')
    #         num2+=1
    # r1_cut=r1[np.argmin(np.abs(cumulative_areas1-1))]
    # r2_cut=r2[np.argmin(np.abs(cumulative_areas2-1))]
    # plt.plot([r1_cut,r1_cut],[0,22],'--',color='gray')
    # plt.plot([r2_cut,r2_cut],[0,22],'--',color='gray')
    plt.legend()
    plt.xlabel('Radius',hawk.font1)
    plt.ylabel('Number of periodic CVs',hawk.font1)
    plt.tick_params(labelsize=16)
    ax.axhline(y=1, color='gray', linestyle='--')
    plt.semilogy()
    plt.semilogx()
    path_out='/Users/baotong/Desktop/aas/GCall/figure/'
    if save:
        plt.savefig(path_out+'CV_bkg.pdf', bbox_inches='tight', pad_inches=0.05)
    if show:
        plt.show()


def plot_aaa():
    z=0.19815
    h=np.linspace(0,8,10000)
    pND=[]
    d=8-h
    pbulge=1.09*np.sqrt(d**2+(z/0.6)**2)**(-1.8)*np.exp(-d**2/1.9-(z/0.6)**2/1.9)
    pdisk=2.5*np.exp(-(3/d)**3-d/2.5-np.abs(z/130))
    for i in range(len(h)):
        pND.append(pnd(d[i]*1000,z))
    plt.plot(h,np.array(pND))
    plt.plot(h,pbulge)
    plt.plot(h,pdisk)
    plt.semilogy()
    plt.ylim(0.01,8.2)
    plt.show()
plot_cum_density(save=1,show=1)