#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import pylab as pl
import string
from datetime import datetime
from scipy.interpolate import lagrange
from scipy import interpolate
from scipy.fftpack import fft,ifft
import scipy.signal as ss
import random
import pandas as pd
import read_csv as data
path_table='/Users/baotong/Desktop/period/table/'
result_NSC=pd.read_excel(path_table+'final_all_del.csv','result_NSC_IG')
result_LW=pd.read_excel(path_table+'final_all_del.csv','result_LW')
result_ND=pd.read_excel(path_table+'final_all_del.csv','result_ND')

ID_LW=result_LW['seq']
P_LW=result_LW['P']
flux_LW=result_LW['flux']
label_LW=result_LW['label']
ra_LW=result_LW['ra']
dec_LW=result_LW['dec']

ID_LW_old = []
P_LW_old = []
flux_LW_old = []

ID_LW_new = []
P_LW_new = []
flux_LW_new = []

for i in range(len(flux_LW)):
    if label_LW[i] == 0:
        ID_LW_new.append(ID_LW[i])
        P_LW_new.append(P_LW[i])
        flux_LW_new.append(flux_LW[i])
    else:
        ID_LW_old.append(ID_LW[i])
        P_LW_old.append(P_LW[i])
        flux_LW_old.append(flux_LW[i])
ID_LW_old =np.array(ID_LW_old)
P_LW_old = np.array(P_LW_old)
flux_LW_old = np.array(flux_LW_old)

ID_LW_new =np.array(ID_LW_new)
P_LW_new = np.array(P_LW_new)
flux_LW_new = np.array(flux_LW_new)

ID_ND=result_ND['seq']
P_ND=result_ND['P']
flux_ND=result_ND['flux']


def getlen(a,b):
    length=((a[0]-b[0])**2+(a[1]-b[1])**2)**0.5
    return length

# path='/Users/baotong/Desktop/li_pr/'
# ACIS_I=np.loadtxt(path+'ACIS-I.txt')
# ACIS_S=np.loadtxt(path+'ACIS-S.txt')

def plot_exptime():
    plt.xlabel('ACIS exposure time (ks)',fontsize=12)
    plt.ylabel('N',fontsize=12,verticalalignment='baseline',horizontalalignment='right',rotation='horizontal')
    plt.ylim(0,22)
    #plt.semilogx()
    plt.hist(ACIS_I,bins=20,histtype='step',color='red')
    plt.hist(ACIS_S,bins=20,histtype='step',color='green')
    plt.hist(LW_I,bins=10,histtype='step',color='blue')
    plt.hist(box_I,bins=10,histtype='step',color='orange')
    plt.legend(['NSC_ACIS-I','NSC_XVP','LW_ACIS-I','N_DISK_ACIS-I'])
    plt.show()
#plot_exptime()

def plot_pop_period():
    P_NSC=data.P_NSC_IG
    ID_NSC=data.ID_NSC_IG
    if len(P_NSC)!=len(ID_NSC) or len(ID_LW)!=len(P_LW):
        print('error')

    plt.semilogx()
    plt.hist(P_NSC,bins=20,histtype='step',color='blue')
    plt.hist(P_LW,bins=10,histtype='step',color='green')

    plt.legend('NSC','LW')
    plt.show()
#plot_pop_period()

def compare_counterpart():
    path='/Users/baotong/Desktop/period/'
    cat1=fits.open(path+'Dong2017t3.fit')
    cat2=fits.open(path+'zhu18_3.fits')

    offset=0.5*1./3600.

    x1=cat1[1].data['raj2000']
    y1=cat1[1].data['dej2000']
    n1=cat1[1].data['ID']

    x2=cat2[1].data['raj2000']
    y2=cat2[1].data['dej2000']
    n2=cat2[1].data['seq']


    ID_NSC = np.array([214, 116, 1502, 1624, 1266, 1206, 2532, 3067, 2525, 3120, 2199, 1133,
             2672, 2157, 2187, 214, 2422, 973, 2841, 2344, 1634, 1677, 2730, 3596])
    x2=x2[ID_NSC-1]
    y2=y2[ID_NSC-1]
    n2=n2[ID_NSC-1]
    match=[[] for i in range(len(n2))]
    for i in range(len(n2)):
        for j in range(len(n1)):
            if getlen([x2[i],y2[i]],[x1[j],y1[j]])<offset:
                match[i].append(n1[j])

    print(match)

#compare_counterpart()

def plot_P_flux():
    P_NSC=result_NSC['P']
    L_NSC=result_NSC['L']
    P_LW=result_LW['P']
    L_LW=result_LW['L31']
    P_NSC=np.array(P_NSC)
    L_NSC=np.array(L_NSC)
    P_LW = np.array(P_LW)
    L_LW=np.array(L_LW)
    # flux_LW = np.array(flux_LW)
    # P_ND=np.array(P_ND)
    # flux_ND=np.array(flux_ND)


    P_min=4944.0/3600.
    P_gap=[7740.0/3600.,11448.0/3600.]

    x=P_gap
    y=[0,1e34]
    plt.figure(1,(10,7))
    plt.text(7900/3600., 1e33, 'period gap')
    plt.fill_between(x,y[1],facecolor='yellow',alpha=0.2)

    # plt.scatter(P_NSC,flux_NSC,color='red')
    # plt.scatter(P_LW,flux_LW,color='green')
    plt.semilogy()
    #plt.semilogx()
    plt.scatter(P_NSC/3600., 1e31*L_NSC, color='red')
    plt.scatter(P_LW/3600., 1e31*L_LW, color='green')
    # plt.scatter(P_NSC_old/3600., L_NSC_old, color='',marker = 'o',edgecolors='r')
    # plt.scatter(P_LW_old/3600., L_LW_old, color='',marker = 'o',edgecolors='g')
    plt.legend(['period_gap','NSC', 'LW'])
    plt.plot([P_min,P_min],[0,1e34],'--')
    plt.text(P_min-0.5,1e34,'period minum')
    plt.xlabel('period(hr)')
    plt.ylabel('Luminosity (erg/s)')

    #plt.savefig(path_table+'P_L.eps')
    plt.show()

plot_P_flux()
#
# plt.hist(P_NSC,histtype = 'step')
# plt.hist(P_LW,histtype = 'step')
# plt.show()

def spatial_P():
    P_NSC = result_NSC['P']
    ra_NSC=result_NSC['ra']
    dec_NSC = result_NSC['dec']
    sgra = [266.4168166, -29.0078250]
    d2sgra_NSC=((ra_NSC - sgra[0]) ** 2 + (dec_NSC - sgra[1])**2)**0.5*3600
    #d2sgra_LW=((ra_LW - sgra[0]) ** 2 + (dec_LW - sgra[1])**2)**0.5*3600
    plt.scatter(dec_NSC,P_NSC)
    #plt.scatter(P_LW,d2sgra_LW)
    plt.show()

#spatial_P()