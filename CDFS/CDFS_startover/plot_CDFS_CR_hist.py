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
import useful_functions as func

path='/Users/baotong/Desktop/CDFS/'
source_info=pd.read_excel(path+'source_infomation.xlsx')
net_CR=source_info['CR']
bkg_CR=source_info['bkg_CR']
plt.figure(1,(7,5))

net_CR=np.array(net_CR)
bkg_CR=np.array(bkg_CR)
# net_CR=net_CR[np.where(net_CR>0)]
def plot_DIS_NET():
    DIS_NET=plt.hist(net_CR*10000,bins=100,linewidth=3,histtype='step',color='green')
    DIS_BKG = plt.hist(bkg_CR * 10000 / 12., bins=100, histtype='step', color='grey')
    print(len(bkg_CR))
    plt.legend(['Net count rate', 'Background count rate'])
    plt.plot([2.8,2.8],[0,200],'--')
    # plt.xlim(-1,30)
    plt.semilogy()
    plt.semilogx()
    bin_bright=DIS_NET[1][np.where(DIS_NET[1]>3)]
    # DIS_bright=plt.hist(net_CR*10000,bin_bright,histtype='step',color='red',facecolor='r',hatch='/',edgecolor='k',fill=True)
    DIS_BKG = plt.hist(bkg_CR * 10000/12., bins=100,histtype='step',color='grey')
    plt.xlabel('Count rate ($10^{-4}$ counts$~$$s^{-1}$ )', func.font1)
    plt.ylabel('Number of sources',func.font1)
    plt.tick_params(labelsize=16)
    plt.savefig(func.figurepath+'src_hist.pdf',bbox_inches='tight',pad_inches=0.0)
    plt.show()
def plot_DIS_bkg():
    plt.semilogy()
    plt.show()

def plot_bright_netp():
    bright_net_CR=net_CR[np.where(net_CR>5e-4)]
    netp=(bkg_CR/12.)/(net_CR+bkg_CR/12.)
    bright_netp=netp[np.where(net_CR>5e-4)]
    print(bright_netp)
plot_DIS_NET()