#!/bin/bash
# -*- coding: utf-8 -*-
# written by Tong
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
from matplotlib import ticker
from astropy.io import fits
from astropy.stats import poisson_conf_interval
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy import units as u
from scipy import optimize as op
import os
import hawkeye as hawk
import NGC104.plot_pXS as plot_pXS
import GC.localCVs as localCVs
import NGC104.CV_model as CV_model

gcname = ['Tuc', 'terzan5', 'M28', 'omg', 'NGC6397', 'NGC6752', 'NGC6266','M30','NGC6121','NGC6304','NGC6656']
# catname = ['xray_properties-592.fits', 'cheng2019_terzan.fit', 'cheng2020_M28.fit',
#            'cheng2020_omg.fit', 'ngc6397_catalog.fits', 'ngc6752_catalog.fits',
#            'NGC6266_p50_i5_src_1_2_4_8.fits']
# gcname = ['M30','NGC6121','NGC6304','NGC6656']

def srcinfo_plt():
    # color_list = ['r', 'g', 'b', 'k', 'orange', 'purple', 'magenta', 'cyan']
    for i in range(len(gcname)):
        path = '/Users/baotong/Desktop/period_' + gcname[i] + '/'
        src_info = np.loadtxt(path + 'src_info.txt')
        counts_all = src_info[:, 3];
        exptime_all = src_info[:, 4];
        print(gcname[i])
        print('src=',len(src_info))
        print('bright_src=', len(np.where(counts_all > 100)[0]))
        # bins_counts = np.logspace(0.5, 4, 20)
        plt.hist(counts_all, bins=20, histtype='step', lw=2, linestyle='-',
                 label=gcname[i])
        plt.semilogx()
        plt.semilogy()

srcinfo_plt()
plt.show()
