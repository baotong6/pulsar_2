#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.io import fits
import sys,os
path='/Users/baotong/Desktop/period/'
catalog=fits.open(path+'zhu18_3.fits')
flux=catalog[1].data['F2-8']
L=flux*8.04523361e+30
HR1=catalog[1].data['HR0']
filter_index=np.where(HR1>0.4)
plt.hist(HR1,bins=30)
plt.show()
fig = plt.figure(2, figsize=(9, 6))
ax1 = fig.add_subplot(211)

bins=np.logspace(np.log10(1e31),np.log10(2e33),100)
allnum = ax1.hist(L, bins=bins, histtype='step', linewidth=2, color='grey')
cvnum = ax1.hist(L[filter_index], bins=bins, histtype='step', linewidth=2, facecolor='r',
                 hatch='/', edgecolor='k', fill=True)
ax1.set_xscale('log')
ax1.set_xlim(1e31, 2e33)
const_frac = len(L[filter_index]) / (len(L))
ax2 = fig.add_subplot(212)
ax2.set_xscale('log')
ax2.plot(bins[:-1], np.zeros(len(bins[:-1])) + const_frac, '-.', color='c')
numofcv = cvnum[0]
numofall = allnum[0]
ax2.step(bins, np.concatenate(([0], numofcv / (numofall + 1e-2))), color='k')
ax2.plot([4e31, 4e31], [0, 1], '--', color='r')
ax2.set_xlim(1e31, 2e33)
plt.show()