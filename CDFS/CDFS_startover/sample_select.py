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
import warnings
import useful_functions as func


def read_catalog(fitsname):
    ra_center=53.1166667
    dec_center=-27.8083333
    catalog=fits.open(fitsname)
    ra=catalog[1].data['RAJ2000']
    dec=catalog[1].data['DEJ2000']
    type=catalog[1].data['OType']
    off_axis=np.sqrt((ra-ra_center)**2+(dec-dec_center)**2)*60
    gamma=catalog[1].data['Gamma']
    # circle = plt.Circle(xy=(ra_center, dec_center),
    #                     radius=8/60, alpha=1,color='black',linestyle='--', fill=False)
    # plt.gcf().gca().add_artist(circle)
    # plt.scatter(ra,dec)
    # plt.show()
    plt.scatter(off_axis,gamma)
    plt.show()

    return off_axis

if __name__=='__main__':
    fitsname='/Users/baotong/Desktop/CDFS/7Ms_catalog.fit'
    print(read_catalog(fitsname))
