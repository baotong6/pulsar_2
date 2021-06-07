import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import pandas as pd
import sys
import os
from scipy import interpolate
from scipy.optimize import curve_fit
import astropy.units as u
import astropy.constants as c
from scipy import interpolate
import stingray as sr
import useful_functions as func


def plot_scatter_single(filename,qpo_P):
    temp=np.loadtxt(filename)
    FP=temp[:,0];period=temp[:,1];
    id=np.where((FP<0.05)&(np.abs(period-qpo_P)<1/10*qpo_P))[0]
    print(len(id))
    plt.scatter(period,FP,marker='x',color='green')
    plt.loglog()
    plt.show()

def run_test():
    epoch1 = '/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep1/CDFS_epoch_ep1.txt'
    epoch2 = '/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep2/CDFS_epoch_ep2.txt'
    epoch3 = '/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep3/CDFS_epoch_ep3.txt'
    epoch4 = '/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep4/CDFS_epoch_ep4.txt'
    epoch_all=[epoch1,epoch2,epoch3,epoch4]
    CR_all=[5e-4,6e-4,7e-4,8e-4,9e-4,1e-3]
    period_all=[1800,3600,7200]
    ep=2;CR=CR_all[5];period=period_all[2]
    path = '/Users/baotong/Desktop/CDFS/simulation/EP{0}/'.format(int(ep + 1))
    # filename=path+'CR_{0}_P_{1}_REJ1034.txt'.format("%.0e"%CR,str(int(period)))

    filename = path + 'CR_{0}_P_{1}_REJ1034_simpleway.txt'.format("%.0e" % CR, str(int(period)))
    # filename = path + 'CR_1e-03_P_1800_REJ1034_three_way.txt'.format("%.0e" % CR, str(int(period)))
    plot_scatter_single(filename,qpo_P=period)

if __name__=='__main__':
    run_test()

