import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
import pandas as pd
import sys
import os
from tkinter import _flatten
import funcs_fits2txt as funcs
from scipy import interpolate
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.stats import poisson_conf_interval

def get_net_info(obsid,ra,dec,srcIDlist):
    ecf = 75
    expmap_file='expmap_02_5_{0}.fits'.format(obsid)
    image = fits.open(path + expmap_file)[0].data
    w = WCS(path + expmap_file)
    (xbin, ybin) = image.shape
    src_x, src_y = w.all_world2pix(ra, dec, 1)
    src_x=src_x.astype('int');src_y=src_y.astype('int')

    not_include_index=[]
    for i in range(len(ra)):
        if src_x[i]<=0 or src_x[i]>=xbin-1 or src_y[i]<=0 or src_y[i]>=xbin-1:
            src_x[i]=0;src_y[i]=0
            not_include_index.append(i)

    expmap=fits.open(path+expmap_file)[0].data.T

    path_txt=path+'txt/txt_psf{0}_{1}/'.format(ecf,obsid)
    src_info=np.loadtxt(path_txt+'src_info.txt')
    src_cts=src_info[:,0];bkg_cts=src_info[:,1]
    net_cts=src_cts-bkg_cts

    expvalue=expmap[src_x,src_y]
    return (net_cts,expvalue)

if __name__=='__main__':
    path = '/Users/baotong/eSASS/data/raw_data/47_Tuc/'
    (ra, dec, srcIDlist) = funcs.read_erosita_cat(filename='/Users/baotong/Desktop/period_Tuc/erosita_cat_coord.xlsx')
    obsIDlist = [700011, 700163, 700013, 700014, 700173, 700174, 700175]
    for i in range(len(obsIDlist)):
        (net_cts,expvalue)=get_net_info(obsIDlist[i],ra,dec,srcIDlist)
        net_CR=net_cts/expvalue/1500
        net_CR=net_CR[np.where(net_CR>0)]
        net_CR=net_CR*1e4
        plt.hist(net_CR,bins=np.logspace(-2,2,30),histtype='step')
        plt.loglog()
        plt.title('obsID: {0}'.format(obsIDlist[i]))
        plt.xlabel('Net count rate')
        plt.ylabel('Number')
        plt.show()
        # net_cts_err = poisson_conf_interval(net_cts, interval='frequentist-confidence').T
        # print(net_cts_err)

        print(len(np.where(net_cts<0)[0]))
    # print(expvalue)
