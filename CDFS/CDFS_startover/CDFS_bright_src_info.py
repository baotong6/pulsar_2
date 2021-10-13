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
from CDFS.CDFS_startover import useful_functions as func
from CDFS.CDFS_startover import sim_psd as sim

bright_source_id=np.array([495,175,716,479,949,855,291,730,711,938,567,208,220,242,
                    200,19,856,168,979,557,996,89,13,785,485,941,306,643,100,
                    943,28,876,102,988,796,947,507,98,81])

path = '/Users/baotong/Desktop/CDFS/'
info = pd.read_excel(path + 'source_infomation.xlsx')
print(len(info))
catalog=fits.open(path+'7Ms_catalog.fit')
Seq=catalog[1].data['Seq']
Gamma=catalog[1].data['Gamma']
NH=catalog[1].data['NH']
def plot_CR_NH_bright():
    gamma_bright_src=Gamma[bright_source_id-1]
    NH_bright_src=NH[bright_source_id-1]/1e22
    label = 'CR'
    CR = info[label]
    CR_bright_src=CR[0:len(bright_source_id)]
    plt.scatter(CR_bright_src*1e4,NH_bright_src,marker='X')
    # plt.semilogx()
    plt.tick_params(labelsize=16)
    plt.xlabel('Count rate ($10^{-4}$ counts$~$$s^{-1}$ )', func.font1)
    plt.ylabel('NH ($10^{22}  cm^{-2}$)',func.font1)
    plt.show()

def plot_CR_NH_all():
    info_sort=info.sort_values(by='ID')
    label = 'CR'
    CR=info[label]
    print(len(CR),len(NH))
    plt.scatter(CR,NH/1e22,marker='X')
    plt.show()


if __name__=='__main__':
    plot_CR_NH_bright()
    plot_CR_NH_all()



