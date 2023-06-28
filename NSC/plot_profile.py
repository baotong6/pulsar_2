#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import random
import os
import pandas as pd
from astropy.stats import poisson_conf_interval
from astropy.coordinates import SkyCoord
from astropy import units as u
import scipy
from stingray.lightcurve import Lightcurve
import hawkeye as hawk
import rocket as rocket
from scipy import optimize as op
from timing_comb import load_data,get_lc_frombkgimg
from plot_result import read_csv as data
font1 = {'family': 'Normal',
         'weight': 'normal',
         'size': 18, }
plt.rc('legend',fontsize=14 )
def adapt_bin(dist,cxb):
    # criteria net>10 and SN>3
    bin_step=0.01
    bin_lf=0;bin_rt=bin_lf+bin_step
    bin_rt_list=[]
    while bin_rt<np.max(dist):
        temp_cts=len(np.where((dist<bin_rt) & (dist>bin_lf))[0])
        temp_area=(bin_rt**2-bin_lf**2)*np.pi
        net_cts=temp_cts-cxb*temp_area
        if net_cts<0 or temp_cts==0:
            bin_rt+=bin_step
            continue
        else:
            SN=net_cts/np.sqrt(temp_cts)
            if net_cts>5 or SN >3:
                bin_rt_list.append(bin_rt)
                bin_lf=bin_rt;bin_rt+=bin_step
                continue
            else:
                bin_rt+=bin_step
                continue

    return bin_rt_list

def f(x,a,b):
    logS=a*x+b  #S~V^a
    return logS
def spectrafit(x,y,error):
    # popt, pcov = op.curve_fit(f, np.log(x), np.log(y),absolute_sigma=True,sigma=np.log(error))
    popt, pcov = op.curve_fit(f, np.log(x), np.log(y))
    perr = np.sqrt(np.diag(pcov))
    logydata=f(np.log(x),popt[0],popt[1])
    ydata=np.exp(logydata)
    return (popt,perr)

def plot_surface_density(save=0,show=1):
    cat=fits.open('/Users/baotong/Desktop/period/zhu18_3.fits')
    ra=cat[1].data['RAJ2000'];dec=cat[1].data['DEJ2000'];seq=cat[1].data['Seq'];H28=cat[1].data['H2-8'];note=cat[1].data['n_Seq']
    seq_bright=seq[np.where(H28>3.5)]
    seq_trans=seq[np.where(note=='m,x')]
    ##--------------------------##
    ra_P=data.ra_NSC_IG;dec_P=data.dec_NSC_IG;seq_P=data.ID_NSC_IG;line_P=data.line_NSC_IG;fore=data.fore
    seq_nolineP=np.array([1748,1083,1080,1873,3242,2268,3370,2478,1182,2961,1180,1854,1487,147,3483,3067,2841])
    # seq_nolineP=np.array([1748,3242,2268,3370,2478,147,3483,442,3067,2841])
    seq_lineP=np.setdiff1d(seq_P,seq_nolineP)
    # seq_lineP=np.array([214,116,1502,2532,1266,1206,2672,1624,2422,3120,1133,2730,2157,1634])
    seq_noP=np.setdiff1d(seq_bright,seq_trans)
    seq_nofore=seq_P[np.where(fore==0)[0]]
    seq_lineP=np.intersect1d(seq_lineP,seq_nofore)
    seq_nolineP=np.intersect1d(seq_nolineP,seq_nofore)
    print(len(seq_lineP),len(seq_nolineP))
    [ra_center,dec_center]=[266.4168250,-29.0077972]
    c3=SkyCoord(266.4197500*u.deg,-29.0040833*u.deg,frame='fk5') ## SWIFT J174540.7-290015
    c1 = SkyCoord(ra * u.deg, dec * u.deg, frame='fk5')
    c2 = SkyCoord(ra_center * u.deg, dec_center * u.deg, frame='fk5')
    dist_all = c1.separation(c2)
    dist_all=dist_all.arcsec
    dist_newtrans=c3.separation(c2).arcsec
    print(dist_newtrans)
    # index_faint = np.where(L058 < 5e30);index_bright = np.where(L058 > 5e30)
    # bin_rt_list=np.logspace(np.log10(2),np.log10(500),15)
    bin_rt_list = np.linspace(3, 500, 14)

    seq_db=seq[np.where(c1.galactic.b>c2.galactic.b)]  ##db>0
    seq_noP=np.intersect1d(seq_noP,seq_db)
    print(len(seq_noP))
    # bin_rt_list=np.concatenate(([0],bin_rt_list))
    # print('bin_rt_list',bin_rt_list)
    x1=[(bin_rt_list[i]+bin_rt_list[i+1])/2 for i in range(len(bin_rt_list)-1)]
    x1err = [(bin_rt_list[i + 1] - bin_rt_list[i]) / 2 for i in range(len(bin_rt_list) - 1)]
    area1 = [(bin_rt_list[i + 1] ** 2 - bin_rt_list[i] ** 2) * np.pi for i in range(len(bin_rt_list) - 1)]
    hist_all=plt.hist(dist_all[seq_noP-1],bin_rt_list)
    plt.close()
    y1=hist_all[0]
    y1_err = np.array(poisson_conf_interval(y1, interval='frequentist-confidence'))
    y1_err[0] = y1 - y1_err[0];y1_err[1] = y1_err[1] - y1

    bin_trans=np.logspace(np.log10(2),np.log10(500),4)
    x2=[(bin_trans[i]+bin_trans[i+1])/2 for i in range(len(bin_trans)-1)]
    x2err = [(bin_trans[i + 1] - bin_trans[i]) / 2 for i in range(len(bin_trans) - 1)]
    area2 = [(bin_trans[i + 1] ** 2 - bin_trans[i] ** 2) * np.pi for i in range(len(bin_trans) - 1)]
    hist_trans=plt.hist(np.concatenate((dist_all[seq_trans-1],[dist_newtrans])),bin_trans)
    plt.close()
    y2 = hist_trans[0]
    y2_err = np.array(poisson_conf_interval(y2, interval='frequentist-confidence'))
    y2*=10;y2_err*=10
    y2_err[0] = y2 - y2_err[0];y2_err[1] = y2_err[1] - y2

    # bin_p=np.logspace(np.log10(20),np.log10(500),5)
    # bin_p=bin_p[1:]
    bin_p = np.logspace(np.log10(10),np.log10(500),5)
    x3=[(bin_p[i]+bin_p[i+1])/2 for i in range(len(bin_p)-1)]
    x3err = [(bin_p[i + 1] - bin_p[i]) / 2 for i in range(len(bin_p) - 1)]
    area3 = [(bin_p[i + 1] ** 2 - bin_p[i] ** 2) * np.pi for i in range(len(bin_p) - 1)]
    hist_pCV=plt.hist(dist_all[seq_lineP-1],bin_p)
    plt.close()
    y3 = hist_pCV[0]
    print(np.sum(y3))
    y3_err = np.array(poisson_conf_interval(y3, interval='frequentist-confidence'))
    y3*=15;y3_err*=15
    y3_err[0] = y3 - y3_err[0];y3_err[1] = y3_err[1] - y3

    bin_p=bin_p[1:]
    x4=[(bin_p[i]+bin_p[i+1])/2 for i in range(len(bin_p)-1)]
    x4err = [(bin_p[i + 1] - bin_p[i]) / 2 for i in range(len(bin_p) - 1)]
    area4 = [(bin_p[i + 1] ** 2 - bin_p[i] ** 2) * np.pi for i in range(len(bin_p) - 1)]
    hist_nolinep=plt.hist(dist_all[seq_nolineP-1],bin_p)
    plt.close()
    y4 = hist_nolinep[0]
    print(np.sum(y4))
    y4_err = np.array(poisson_conf_interval(y4, interval='frequentist-confidence'))
    print(np.sum(y4))
    y4*=10;y4_err*=10
    y4_err[0] = y4 - y4_err[0];y4_err[1] = y4_err[1] - y4


    y1=y1/area1;y2=y2/area2;y3=y3/area3;y4=y4/area4;
    y1_err=y1_err/area1;y2_err=y2_err/area2;y3_err=y3_err/area3;y4_err=y4_err/area4

    plt.figure(1,(10,8))
    plt.errorbar(x1,y1,xerr=x1err,yerr=y1_err,fmt='ro', capsize=1, elinewidth=1, ecolor='r', color='r',markersize=4,label='Total')
    plt.errorbar(x2,y2,xerr=x2err,yerr=y2_err,fmt='o', capsize=1, elinewidth=4, ecolor='orange', color='orange',markersize=1,label='Transient x 10')
    plt.errorbar(x3,y3,xerr=x3err,yerr=y3_err,fmt='o', capsize=3, elinewidth=2, ecolor='c', color='c',markersize=5,label='Periodic CV x 15')
    plt.errorbar(x4,y4,xerr=x4err,yerr=y4_err,fmt='o', capsize=2, elinewidth=3, ecolor='grey', color='grey',markersize=5,label='Other periodic sources x 10')
    (popt, perr) = spectrafit(x1, y1, y1_err[0])
    # x1 = np.concatenate(([2], x1))
    # plt.plot(x1, np.exp(f(np.log(x1), popt[0], popt[1])),'--',linewidth=1,color='red')

    (popt, perr) = spectrafit(x2, y2, y2_err[0])
    print(popt,perr)
    x2 = np.concatenate(([2], x2))
    plt.plot(x2, np.exp(f(np.log(x2), popt[0], popt[1])),'-',linewidth=2,color='orange')
    plt.text(3, 1e-1, r'$\Gamma={0:.1f}\pm{1:.1f}$'.format(popt[0],perr[0]), color='orange',fontdict=font1)
    (popt, perr) = spectrafit(x3, y3, y3_err[0])
    print(popt, perr)
    # x3= np.concatenate(([10], x3))
    plt.plot(x3, np.exp(f(np.log(x3), popt[0], popt[1])),'-.',linewidth=1,color='c')
    plt.text(10,5e-2,r'$\Gamma={0:.1f}\pm{1:.1f}$'.format(popt[0],perr[0]),color='c',fontdict=font1)
    (popt, perr) = spectrafit(x4, y4, y4_err[0])
    print(popt,perr)
    x4 = np.concatenate(([20], x4))
    plt.plot(x4, np.exp(f(np.log(x4), popt[0], popt[1])),'-.',linewidth=2,color='grey')
    plt.text(30,1e-2,r'$\Gamma={0:.1f}\pm{1:.1f}$'.format(popt[0],perr[0]),color='grey',fontdict=font1)
    plt.legend()
    plt.loglog()
    plt.tick_params(labelsize=18)
    plt.xlabel('R (arcsec)',font1)
    plt.ylabel(r'$\rm Surface~Number~Density~(arcsec^2)$',font1)
    plt.savefig('SD.pdf', bbox_inches='tight', pad_inches=0.05)
    plt.show()

# def select_src_noline():
#     cat=fits.open('/Users/baotong/Desktop/period/zhu18_3.fits')
#     ra=cat[1].data['RAJ2000'];dec=cat[1].data['DEJ2000'];seq=cat[1].data['Seq'];H28=cat[1].data['H2-8'];note=cat[1].data['n_Seq']
#     [ra_center,dec_center]=[266.4168250,-29.0077972]
#     c1 = SkyCoord(ra * u.deg, dec * u.deg, frame='fk5')
#     c2 = SkyCoord(ra_center * u.deg, dec_center * u.deg, frame='fk5')
#     dist_all = c1.separation(c2)
#     dist_all=dist_all.arcsec
#
#     seq_nolineP = np.array([1748, 1083, 1080, 1873, 3242, 2268, 3370, 2478, 1182, 2961, 1180, 1854, 1487, 147, 3483])
#     distance_nolineP=dist_all[seq_nolineP-1]
#     flux28range=H28[seq_nolineP-1]
#     seq_inrange=seq[np.where((H28<=np.max(flux28range))&(H28>=np.min(flux28range)))]
#     print(len(seq_inrange))
#     chosen_seq=random.sample(list(seq_inrange),len(seq_nolineP))
#     chosen_seq=np.array(chosen_seq)
#     plt.hist(distance_nolineP, bins=np.logspace(np.log10(0.3),np.log10(500),30), histtype='step')
#     plt.hist(dist_all,bins=np.logspace(np.log10(0.3),np.log10(500),30),histtype='step')
#     plt.loglog()
#     plt.show()
#     # plt.hist(H28, bins=np.logspace(np.log10(0.3),np.log10(500),30), histtype='step')
#     # plt.hist(H28[chosen_seq-1],bins=6,histtype='step')
#     # plt.loglog()
#     # plt.show()
def plot_M_dist():
    P_min = 1.373333333
    P_gap = [7740.0 / 3600., 11448.0 / 3600.]

    cat=fits.open('/Users/baotong/Desktop/period/zhu18_3.fits')
    ra=cat[1].data['RAJ2000'];dec=cat[1].data['DEJ2000'];seq=cat[1].data['Seq'];H28=cat[1].data['H2-8'];note=cat[1].data['n_Seq']
    [ra_center,dec_center]=[266.4168250,-29.0077972]
    c1 = SkyCoord(ra * u.deg, dec * u.deg, frame='fk5')
    c2 = SkyCoord(ra_center * u.deg, dec_center * u.deg, frame='fk5')
    dist_all = c1.separation(c2)
    dist_all=dist_all.arcsec
    seq_nolineP=np.array([1748,1083,1080,1873,3242,2268,3370,2478,1182,2961,1180,1854,1487,147,3483,3067,2841])
    period_nolineP=np.array([592.8268,790.47634,813.75534,1377.72675,1535.0765,2416.97898,4669.6241,
                             4819.50937,4881.85901,4967.46312,5691.19572,6793.53256,13482.54011,19367.0837,
                             30860.32616,5729.01747,6189.86842])

    seq_lineP=np.array([214,2560,116,1502,2380,2532,1206,2672,1624,2508,3357,
                        1219,2422,1853,3120,1133,2730,1084,2157,1634,3401])
    period_lineP=np.array([320.91053,605.57126,666.16334,971.88026,1092.25952,5317.45188,5611.67228,6460.0323,9905.89401,
                           10079.52747,11301.56338,12815.25528,13190.69792,13581.78954,13809.67506,16734.44115,
                           16781.33915,19237.05827,33125.74533,43880.81969,44060.62742])
    L_nolineP=np.array([3.459321519,1.477715839,0.539996479,3.057672226,2.054901118,6.025469999,3.372026107,
                        2.629607492,4.364062166,6.743586436,6.922945751,1.142854463,1.054122124,
                        3.493742786,1.999359674,6.50415569,28.64777707])

    L_lineP=np.array([44.99776346,96.0706941,104.4698397,134.119192,71.10349172,
                      77.53377479,62.05702167,7.194339503,44.69829201,17.32171658,12.83189888,
                      49.3731605,212.0757094,42.84409054,18.49655842,38.86788673,39.26367339,35.65257716,30.02554703,8.012922244,202.57737])
    distance_nolineP=dist_all[seq_nolineP-1]
    distance_lineP=dist_all[seq_lineP-1]
    bins_p=np.logspace(np.log10(0.1),np.log10(10),16)

    # plt.hist(period_lineP/3600,bins=bins_p,histtype='step',color='c',lw=2, linestyle='-',facecolor='c',
    #      hatch='/', edgecolor='k',fill=True,label='Periodic CVs')
    # plt.hist(period_nolineP/3600,bins=bins_p,histtype='step',color='grey',lw=2, linestyle='-',label='Other periodic sources')
    # plt.scatter(period_lineP/3600,distance_lineP,marker='^',color='c',s=50,label='Periodic CVs')
    # plt.scatter(period_nolineP/3600,distance_nolineP,marker='+',color='grey',s=50,label='Other periodic sources')
    plt.scatter(period_lineP/3600,L_lineP*1e31,marker='^',color='c',s=50,label='Periodic CVs')
    plt.scatter(period_nolineP/3600,L_nolineP*1e31,marker='+',color='grey',s=50,label='Other periodic sources')
    plt.xlabel('Period (hour)',font1)
    # plt.ylabel(r'$\rm R (arcsec)$',font1)
    # plt.ylabel(r'$\rm Number~of~sources~per~bin$',font1)
    plt.ylabel(r'$\rm X-ray Luminosity (erg~s^{-1})$',font1)
    plt.tick_params(labelsize=18)
    plt.semilogx()
    plt.semilogy()
    plt.legend()
    plt.show()

if __name__=="__main__":
    # plot_surface_density(save=0, show=1)
    # select_src_noline()
    plot_M_dist()