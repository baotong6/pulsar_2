import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
import pandas as pd
import sys
import os
from tkinter import _flatten
import funcs_fits2txt as funcs
from funcs_fits2txt import Circle
from scipy import interpolate
from astropy.coordinates import SkyCoord
from astropy import units as u

ra_center=6.022318
dec_center=-72.081443
inter_radius=20  ##角分

path = '/Users/baotong/eSASS/data/raw_data/47_Tuc/'
ID = [700011, 700013, 700014, 700163, 700173, 700174, 700175]
catalog_file = '/Users/baotong/Desktop/period_Tuc/erosita_cat_coord.xlsx'
(ra_eR,dec_eR,srcIDlist)=funcs.read_erosita_cat(catalog_file)

c1 = SkyCoord(ra_eR * u.deg, dec_eR * u.deg, frame='fk5')
c2 = SkyCoord(ra_center * u.deg, dec_center * u.deg, frame='fk5')
dist=c1.separation(c2)
dist=dist.arcmin
inter_srcID=srcIDlist[np.where(dist<inter_radius)[0]]
np.savetxt(path+'txt/inter_srcID.txt',inter_srcID,fmt='%10d')
print(inter_srcID)