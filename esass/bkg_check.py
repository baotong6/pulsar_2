import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
import pandas as pd
import sys
import os
from tkinter import _flatten
import esass.funcs_fits2txt as funcs
from esass.funcs_fits2txt import Circle
from scipy import interpolate
from astropy.coordinates import SkyCoord
from astropy import units as u

obsIDlist = [700011, 700163, 700013, 700014, 700173, 700174, 700175]
expT = [25808.86, 25268.79, 25208.83, 25208.82, 8809.16, 8389.71, 8330.68]
def make_srfinfo_from_bkgimg(ecf):
    path='/Users/baotong/eSASS/data/raw_data/47_Tuc/'
    os.chdir(path)
    srcid = np.arange(1, 889, 1)
    for obsid in obsIDlist:
        bkgimgfile=fits.open('bkg_map_02_5_{0}.fits'.format(obsid))
        psffile=fits.open('psfmap_{0}.fits'.format(obsid))
        psfmap_all = psffile[0].data
        ecflist = np.array([0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95])
        indexpsf = np.where(ecflist == ecf)[0][0]
        psfmap_ecf = psfmap_all[indexpsf]
        psfmap_ecf = psfmap_ecf.T
        area_file=np.loadtxt('./txt/txt_psf{0}_{1}/src_area.txt'.format(int(ecf*100),obsid))
        src_area=area_file[:,1]
        src_area/=80*80  ##转回img 坐标


