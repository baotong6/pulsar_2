#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.timeseries import LombScargle

font1 = {'family': 'Normal',
         'weight': 'normal',
         'size': 18, }
font2 = {'family': 'Normal',
         'weight': 'normal',
         'size': 16, }

def get_LS(time, flux,freq,outpath,outname,save=False,show=True):
    x = time
    y = flux
    LS = LombScargle(x, y,normalization = 'standard')
    power = LS.power(freq)
    max_NormLSP=np.max(power)
    FP=LS.false_alarm_probability(power.max(),minimum_frequency = freq[0],maximum_frequency = freq[-1],method='baluev')
    FP_99 = LS.false_alarm_level(0.0027,minimum_frequency = freq[0], maximum_frequency = freq[-1],method='baluev')
    FP_95 = LS.false_alarm_level(0.05, minimum_frequency=freq[0],
                                 maximum_frequency=freq[-1], method='baluev')
    FP_68 = LS.false_alarm_level(0.32,minimum_frequency=freq[0],
                                 maximum_frequency=freq[-1], method='baluev')

    plt.plot(freq, power)
    plt.semilogx()
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
    # plt.show()
    if save:
        plt.savefig(outpath + outname + '.eps')
    if show:plt.show()
    else:plt.close()
    return [FP,out_period,max_NormLSP]