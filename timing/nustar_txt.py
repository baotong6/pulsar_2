#!/bin/bash
# -*- coding: utf-8 -*-
# written by Tong
# extract the arrival time of photons in nustar data
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import string
#import correct as correct
from scipy.interpolate import lagrange
from scipy import interpolate
from scipy.fftpack import fft,ifft
import scipy.signal as ss
import random
from astropy.timeseries import LombScargle
import pandas as pd
#path='/Volumes/pulsar/xmm_obs_16/'
# path='/Volumes/pulsar/period/'
path='/Volumes/pulsar/WR/40012018002/'
src_file='nu40012018002A01_sr.lc'
srclc=fits.open(path+'products/'+src_file)
time=srclc[1].data.field(0)
rate=srclc[1].data.field(1)

def get_txt():
    evt='nu40012018002A01_cl_bary.evt'
    hdul_evt = fits.open(path +  'out/' + 'nu40012018002A01_cl_bary.evt')

    x=hdul_evt[1].data.field(12)
    y=hdul_evt[1].data.field(13)
    PI=hdul_evt[1].data.field(7)
    time=hdul_evt[1].data.field(0)
    energy=(PI * 0.04 + 1.6)*1000

    # tstart=hdul_evt[1].header['TSTART']
    # tstop=hdul_evt[1].header['TSTOP']

    t_start=hdul_evt[2].data.field(0)
    t_stop=hdul_evt[2].data.field(1)
    print(t_start)

    reg=[583.9681,659.83706,13.539279]

    # #--------------------------------------------------
    def where_region(x,y,reg):
        #输入所有光子x,y，输出region之中光子的index
        r=np.array((x-reg[0],y-reg[1]))
        len_r=np.sqrt(r[0]**2+r[1]**2)
        temp=len_r-reg[2]
        return np.where(temp<=0)



    src_index=where_region(x,y,reg)
    src_x=x[src_index]
    src_y=y[src_index]
    src_t=time[src_index]
    src_E=energy[src_index]
    #src_ID=obs_ID[src_index]
    #src_***就是该源的所有光子的信息

    def delete_photon_ID(time,energy):
        i=0
        while i < len(energy):
            if energy[i]>50000 or energy[i]<3000:
                energy=np.delete(energy,i)
                time=np.delete(time,i)
                i=i-1
            i=i+1
        return [time,energy]

    [src_t, src_E] = delete_photon_ID(src_t, src_E)

    src_txt = np.column_stack((src_t, src_E))
    src_txt = src_txt[src_txt[:, 0].argsort()]
    epoch = np.column_stack((t_start, t_stop))
    # print src_txt
    os.chdir(path)
    os.system('mkdir txt')
    os.system('rm ./txt/*.txt')
    np.savetxt(path + '/txt/' +  'WR1' + '.txt', src_txt, fmt = "%.7f  %.3f ")
    print(epoch)
    np.savetxt(path + '/txt/' + 'epoch_'  + 'WR1' + '.txt', epoch, fmt = "%.2f  %.2f ")

#get_txt()

def get_fig_LS():
    plt.scatter(time,rate)
    plt.show()
    border = 10000
    vary = np.array([i for i in range(0, border)])
    freq = 1 / 50000. + vary * 1.e-8
    x=time
    y=rate
    LS=LombScargle(x,y,dy=1,normalization='standard',fit_mean=True,
                   center_data=True).power(freq,method='cython')
    FP_99 = LombScargle(x, y, dy = 1, normalization = 'standard',
                        fit_mean = True, center_data = True).false_alarm_level(
        0.01,
        minimum_frequency = freq[0], maximum_frequency = freq[-1])
    FP_90 = LombScargle(x, y, dy = 1, normalization = 'standard',
                        fit_mean = True, center_data = True).false_alarm_level(
        0.1,
        minimum_frequency = freq[0], maximum_frequency = freq[-1])
    FP_68 = LombScargle(x, y, dy = 1, normalization = 'standard',
                        fit_mean = True, center_data = True).false_alarm_level(
        0.32,
        minimum_frequency = freq[0], maximum_frequency = freq[-1])

    plt.semilogx()
    plt.step(freq, LS)
    print(max(LS))
    print()
    plt.plot([freq[0], freq[-1]], [FP_99, FP_99], '--')
    plt.plot([freq[0], freq[-1]], [FP_90, FP_90], '--')
    plt.plot([freq[0], freq[-1]], [FP_68, FP_68], '--')
    plt.show()

#get_fig_LS()
def pfold(time, P, flux, bin):
    def trans(t, p_test, shift = 0.0):
        ti = t
        v = 1.0 / p_test
        turns = v * ti
        turns += shift
        # 初始相位
        for i in range(len(turns)):
            turns[i] = turns[i] - int(turns[i])
        return turns


    time=np.loadtxt(path + '/txt/' +  'WR1' + '.txt')[:,0]

    turns = trans(time, P)
    loc=np.zeros(bin)
    for index in turns:
        loc[int(index*bin)] += 1

    x = np.array([(i / bin + 0.5 / bin) for i in range(bin)])
    x2 = np.concatenate((x, x + 1))
    y2 = np.concatenate((loc, loc))
    plt.step(x2, y2, color = 'red')

    # plt.errorbar(turns, flux[0], yerr = [flux[0] - flux[1], flux[2] - flux[0]],
    #              fmt = 'o', capsize = 3, elinewidth = 1, color = 'red', ecolor = 'red')
    # plt.errorbar(turns + 1., flux[0], yerr = [flux[0] - flux[1], flux[2] - flux[0]],
    #              fmt = 'o', capsize = 3, elinewidth = 1, color = 'red', ecolor = 'red')
    # plt.step(turns,flux)
    plt.show()

#get_txt()
pfold(time, 15266.94, rate,bin=20)