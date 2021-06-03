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
from astropy.stats import poisson_conf_interval

path='/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8/'
star_ID=np.array([991,950,925,910,872,780,518,457,398,319,153,64])
for i in range(len(star_ID)):
    srcevt = np.loadtxt(path + '{0}.txt'.format(star_ID[i]))
    bkgevt = np.loadtxt(path + '{0}_bkg.txt'.format(star_ID[i]))
    epoch=np.loadtxt(path+'epoch_src_{0}.txt'.format(star_ID[i]))
    t_start = epoch[:, 0]
    t_end = epoch[:, 1]
    obsID = epoch[:, 2]
    expT = epoch[:, 3]



    cts = [];bkg_cts=[]
    for k in range(len(obsID)):
        cts.append(len(np.where(srcevt[:, 2] == obsID[k])[0]))
        bkg_cts.append(len(np.where(bkgevt[:,2] == obsID[k])[0]))
    cts = np.array(cts);bkg_cts=np.array(bkg_cts)
    bkg_cts=bkg_cts/12.
    net_cts=cts-bkg_cts
    net_cts[np.where(net_cts<0)]=0
    net_cts_all=np.sum(net_cts)

    y2_err=np.array(poisson_conf_interval(net_cts,interval='frequentist-confidence'))
    y2_err[0]=net_cts-y2_err[0]
    y2_err[1]=y2_err[1]-net_cts

    CR = net_cts / expT
    err_CR=y2_err/expT
    VI=np.max(CR)/np.min(CR[np.where(CR>0)])

    # plt.title('Star {0};VI {1}'.format(star_ID[i],VI))
    plt.title('Star {0}; net_counts={1}'.format(star_ID[i],int(net_cts_all)))
    # plt.plot(t_start, CR, marker='+')
    plt.ylim(1e-7,5e-3)
    plt.errorbar(t_start, CR, yerr=err_CR, fmt='.', capsize=1, elinewidth=1, ecolor='red')
    plt.semilogy()
    plt.show()


