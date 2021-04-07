import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import pandas as pd
import sys
import os
from scipy import interpolate
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from tkinter import _flatten
from astropy.stats import poisson_conf_interval

path='/Users/baotong/Desktop/CDFS/'
os.chdir(path)
epoch_file=np.loadtxt('CDFS_epoch.txt')

epoch1=[55e6,100e6]
epoch2=[300e6,320e6]
epoch3=[380e6,400e6]
epoch4=[500e6,580e6]
def plot_epoch_obs():
    tstart=epoch_file[:,0]
    tstop=epoch_file[:,1]
    obsid=epoch_file[:,2]
    print(np.sum(epoch_file[:,3]))

    x=(tstart+tstop)/2
    err_x=(tstop-tstart)/2
    y=np.zeros(len(x))+1
    plt.errorbar(x,err_x,y,marker='o',fmt='--')
    # plt.loglog()
    plt.show()
# plot_epoch_obs()

def plot_exp_epoch(k):
    path = '/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep{0}/'.format(k)
    epoch = np.loadtxt(path + 'CDFS_epoch_ep{0}.txt'.format(k))
    ID=epoch[:,2]
    exptime=np.sum(epoch[:,3])
    ID=ID.astype('int')
    print(exptime)
    # print(ID)

plot_exp_epoch(1)
plot_exp_epoch(2)
plot_exp_epoch(3)
plot_exp_epoch(4)