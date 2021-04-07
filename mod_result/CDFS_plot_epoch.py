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

    x=(tstart+tstop)/2
    err_x=(tstop-tstart)/2
    y=np.zeros(len(x))+1
    plt.errorbar(x,err_x,y,marker='o',fmt='--')
    plt.loglog()
    plt.show()
plot_epoch_obs()