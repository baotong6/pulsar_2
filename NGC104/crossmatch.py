import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import pandas as pd
import sys
import os
from scipy import interpolate
from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import SkyCoord,match_coordinates_sky

def load_data():
    ## for LW x-ray and IR
    path='/Users/baotong/Desktop/period_LW/'
    xray_info=pd.read_excel(path+'final_all_del_add.xlsx','result_LW')
    ra_x = np.array(xray_info['ra'])
    dec_x = np.array(xray_info['dec'])
    seq_x = np.array(xray_info['seq'])

    IR_info=pd.read_excel(path+'IR_match_kumiko.xlsx')
    ra_IR=np.array(IR_info['RA_IR'])
    dec_IR=np.array(IR_info['DEC_IR'])

    return (ra_x,dec_x,ra_IR,dec_IR)

(ra_x,dec_x,ra_IR,dec_IR)=load_data()
c1 = SkyCoord(ra=ra_x * u.degree, dec=dec_x * u.degree)
c2 = SkyCoord(ra=ra_IR * u.degree, dec=dec_IR * u.degree)
plt.scatter(ra_x,dec_x)
plt.scatter(ra_IR,dec_IR)

idx, d2d, d3d = match_coordinates_sky(matchcoord=c1, catalogcoord=c2, nthneighbor=1)
max_sep = 5.0 * u.arcsec
sep_constraint = d2d < max_sep
c_matches = c1[sep_constraint]
catalog_matches = c2[idx[sep_constraint]]
d2d_matches = d2d[sep_constraint]
print(idx[sep_constraint])

plt.scatter(ra_IR[idx[sep_constraint]],dec_IR[idx[sep_constraint]],marker='+',s=20)
print(c_matches)
print(catalog_matches)
print(d2d_matches.arcsec)

plt.show()