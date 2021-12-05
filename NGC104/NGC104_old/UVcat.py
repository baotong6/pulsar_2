import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import pandas as pd
import sys
import os
from scipy import interpolate
from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import SkyCoord
path='/Users/baotong/Desktop/period_Tuc/'

def trans_coord():
    ##useless
    path='/Users/baotong/Desktop/period_Tuc/ngc0104/'
    filename='hlsp_hugs_hst_wfc3-uvis_ngc0104_f336w_v1_stack-0225s.fits'
    img_file=fits.open(path+filename)
    img_data=img_file[0].data
    w = WCS(path + filename)
    # src_x, src_y = w.all_world2pix(ra, dec, 1)
    ra,dec=w.all_pix2world(5000,5000,1)

def read_uv_cat(filename):
    cat_dat=np.loadtxt(filename)
    # cat_dat=pd.read_csv(filename)
    ra=cat_dat[:,33];dec=cat_dat[:,34]

    return (ra,dec)

def read_chandra_cat(filename):
    cat_file=fits.open(filename)
    cat_data=cat_file[1].data
    ra=cat_data['RAdeg']
    dec=cat_data['DEdeg']

    return (ra,dec)

def read_erosita_cat(filename):
    cat=pd.read_excel(filename,header=0)
    ra_hms=cat['RA']
    dec_hms=cat['DEC']
    ra=[];dec=[]
    for i in range(len(dec_hms)):
        skycoord=ra_hms[i]+dec_hms[i]
        c = SkyCoord(skycoord, unit=(u.hourangle, u.deg))
        ra.append(c.ra.value)
        dec.append(c.dec.value)

    return (ra,dec)
def plot_three_cat():
    (ra1,dec1)=read_chandra_cat(path+'xray_properties-592.fits')
    (ra2,dec2)=read_erosita_cat(path+'erosita_cat_coord.xlsx')
    (ra3,dec3)=read_uv_cat(path+'ngc0104/ngc104_meth1.txt')
    plt.xlabel('RA (deg)')
    plt.ylabel('DEC (deg)')

    plt.scatter(ra1, dec1, marker='o',color='red')
    plt.scatter(ra2, dec2, marker='.', color='green')
    plt.scatter(ra3, dec3, marker='.', color='blue')
    plt.legend(['Chandra','eROSITA','HST_UV'])
    plt.show()

if __name__=='__main__':
    plot_three_cat()