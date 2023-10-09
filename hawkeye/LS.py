#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.timeseries import LombScargle
from astropy.io import fits
import stingray as sr
from stingray.events import EventList
from stingray.lightcurve import Lightcurve
from stingray import Lightcurve, Crossspectrum, sampledata,Powerspectrum,AveragedPowerspectrum
from stingray.simulator import simulator, models
import rednoise
from mcmc.gptLS import lomb_scargle
from rednoise import useful_functions as func

font1 = {'family': 'Normal',
         'weight': 'normal',
         'size': 18, }
font2 = {'family': 'Normal',
         'weight': 'normal',
         'size': 16, }


def get_LS(time, flux,freq,outpath=None,outname=None,save=False,show=True):
    x = time
    y = flux
    LS = LombScargle(x, y,normalization = 'standard')
    power = LS.power(freq)
    max_NormLSP=np.max(power)
    period_peak=1./freq[np.where(power==np.max(power))][0]
    FP=LS.false_alarm_probability(power.max(),minimum_frequency = freq[0],maximum_frequency = freq[-1],method='baluev')
    Np=1000
    FAP_N=1-(1-FP)**Np
    print(FAP_N)
    FP_99 = LS.false_alarm_level(0.0027,minimum_frequency = freq[0], maximum_frequency = freq[-1],method='baluev')
    FP_95 = LS.false_alarm_level(0.05, minimum_frequency=freq[0],
                                 maximum_frequency=freq[-1], method='baluev')
    FP_68 = LS.false_alarm_level(0.32,minimum_frequency=freq[0],
                                 maximum_frequency=freq[-1], method='baluev')
    plt.figure(1, (6, 6))
    # plt.title('Period={0:.2f}'.format(period_peak), font1)
    plt.text(freq[np.where(power==np.max(power))][0]*1.3,max_NormLSP*0.95,'P={0:.2f}s'.format(period_peak),fontsize=18,fontweight='semibold')
    plt.plot(freq, power)
    plt.semilogx()
    # plt.semilogy()
    out_period=1./freq[np.where(power==np.max(power))][0]
    plt.plot([freq[0], freq[-1]], [FP_99, FP_99], '--')
    plt.plot([freq[0], freq[-1]], [FP_95, FP_95], '--')
    plt.plot([freq[0], freq[-1]], [FP_68, FP_68], '--')
    plt.text(freq[0], FP_99, '1-FAP 99.73%',font1)
    plt.text(freq[0], FP_95, '95%',font1)
    plt.text(freq[0], FP_68, '68%',font1)
    plt.xlabel('Frequency (Hz)',font1)
    plt.ylabel('Normalized LS Periodogram',font1)
    plt.tick_params(labelsize=16)
    if save:
        plt.savefig(outpath + outname + '_LS.pdf',bbox_inches='tight', pad_inches=0.01)
    if show:plt.show()
    else:plt.close()
    return [FP,out_period,max_NormLSP]

def get_gptLS(time, flux,freq,outpath=None,outname=None,save=False,show=True):
    from numba import jit
    @jit(nopython=True)
    def lombscargle_jit(x, y, f, block_size):
        n = len(x)
        m = len(f)
        power = np.zeros(m)

        for i in range(0, n, block_size):
            x_block = x[i:i + block_size]
            y_block = y[i:i + block_size]
            sin_omega_x = np.sin(np.outer(2 * np.pi * f, x_block))
            cos_omega_x = np.cos(np.outer(2 * np.pi * f, x_block))
            S = np.sum(sin_omega_x ** 2, axis=1)
            C = np.sum(cos_omega_x ** 2, axis=1)
            Sy = y_block @ sin_omega_x
            Cy = y_block @ cos_omega_x
            block_power = 0.5 * ((Sy ** 2 / C) + (Cy ** 2 / S))
            power += block_power

        return power / n

    def lombscargle(x, y, f, block_size=100000):
        # 优化输入数据类型来提高性能
        x = np.asarray(x, dtype=np.float64)
        y = np.asarray(y, dtype=np.float64)
        f = np.asarray(f, dtype=np.float64)
        # 通过广播来扩展 f 的维度，以保证 sin_omega_x 和 cos_omega_x 的正确形状
        # 分块计算
        power = lombscargle_jit(x, y, f, block_size=block_size)
        return power
    power=lombscargle(time, flux, freq)
    max_NormLSP=np.max(power)
    period_peak=1./freq[np.where(power==np.max(power))][0]
    # FP=LS.false_alarm_probability(power.max(),minimum_frequency = freq[0],maximum_frequency = freq[-1],method='baluev')
    # FP_99 = LS.false_alarm_level(0.0027,minimum_frequency = freq[0], maximum_frequency = freq[-1],method='baluev')
    # FP_95 = LS.false_alarm_level(0.05, minimum_frequency=freq[0],
    #                              maximum_frequency=freq[-1], method='baluev')
    # FP_68 = LS.false_alarm_level(0.32,minimum_frequency=freq[0],
    #                              maximum_frequency=freq[-1], method='baluev')
    plt.figure(1, (6, 6))
    # plt.title('Period={0:.2f}'.format(period_peak), font1)
    plt.text(freq[np.where(power==np.max(power))][0]*1.3,max_NormLSP*0.95,'P={0:.2f}s'.format(period_peak),fontsize=18,fontweight='semibold')
    plt.plot(freq, power)
    plt.semilogx()
    # plt.semilogy()
    out_period=1./freq[np.where(power==np.max(power))][0]
    # plt.plot([freq[0], freq[-1]], [FP_99, FP_99], '--')
    # plt.plot([freq[0], freq[-1]], [FP_95, FP_95], '--')
    # plt.plot([freq[0], freq[-1]], [FP_68, FP_68], '--')
    # plt.text(freq[0], FP_99, '1-FAP 99.73%',font1)
    # plt.text(freq[0], FP_95, '95%',font1)
    # plt.text(freq[0], FP_68, '68%',font1)
    plt.xlabel('Frequency (Hz)',font1)
    plt.ylabel('Normalized LS Periodogram',font1)
    plt.tick_params(labelsize=16)
    if save:
        plt.savefig(outpath + outname + '_LS.pdf',bbox_inches='tight', pad_inches=0.01)
    if show:plt.show()
    else:plt.close()
    return [power,out_period,max_NormLSP]

def FAPn():
    F=0.0001;Np=np.arange(100,1000,1)
    FAPn=1-(1-F)**Np
    # plt.plot(Np,FAPn,color='red',linestyle='--')
    # plt.show()
    # print(FAPn)
if __name__=='__main__':
    # path = '/Users/baotong/Downloads/hzq2/'
    # time=np.load(path+'time_h.npy')
    # flux=np.load(path+'flux_h.npy')
    # freq=np.linspace(12,1e4,10000)
    # lc=Lightcurve(time,flux)
    # psd=rednoise.func.plot_psd(lc)
    # get_LS(time,flux,freq)
    path='/Users/baotong/Desktop/XMMcentral/all_lc/'
    file_name='0886081301_268.80_26.05_2976_pnsrc_2_10keV.lc'
    a = fits.open(path + file_name)
    rate = a[1].data['RATE']
    time = a[1].data['TIME']
    T = time[-1] - time[0]
    lc = Lightcurve(time=time, counts=rate * (time[1] - time[0]))
    freq = np.arange(1 / T, 1 / 1, 1e-5)
    # hawk.get_LS(lc.time,lc.counts,freq,show=1)
    CR = np.mean(lc.counts) / lc.dt
    epoch = np.array([time[0], time[-1], '11111', T])
    get_LS(time, rate, freq=np.arange(1 / T, 1 / 1, 1e-6), outpath=None, outname=None, save=0, show=1)
    # FAPn()