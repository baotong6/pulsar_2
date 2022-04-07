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
    path='/Users/baotong/Desktop/period_Tuc/'
    xray_info=pd.read_excel(path+'result_0.5_8_all.xlsx','47Tuc')
    ra_x = np.array(xray_info['RA'])
    dec_x = np.array(xray_info['DEC'])
    seq_x = np.array(xray_info['Seq'])

    IR_info=fits.open(path+'Cheng2019.fit')
    ra_IR=IR_info[1].data['RAJ2000']
    dec_IR=IR_info[1].data['DEJ2000']

    return (ra_x,dec_x,ra_IR,dec_IR)

def match_twofits():
    path='/Users/baotong/Desktop/period_Tuc/'
    file1=fits.open(path+'xray_properties-592.fits')
    ra_1=file1[1].data['RAdeg']
    dec_1=file1[1].data['DEdeg']

    file2=fits.open(path+'Cheng2019.fit')
    ra_2=file2[1].data['RAJ2000']
    dec_2=file2[1].data['DEJ2000']

    c1 = SkyCoord(ra=ra_2 * u.degree, dec=dec_2 * u.degree)
    c2 = SkyCoord(ra=ra_1 * u.degree, dec=dec_1 * u.degree)
    idx, d2d, d3d = match_coordinates_sky(matchcoord=c1, catalogcoord=c2, nthneighbor=1)

    max_sep = 0.1 * u.arcsec
    sep_constraint = d2d < max_sep
    c_matches = c1[sep_constraint]
    catalog_matches = c2[idx[sep_constraint]]
    d2d_matches = d2d[sep_constraint]
    print(idx)
    judge=np.array(sep_constraint).astype('int')
    print(idx[sep_constraint])
    match_result=np.column_stack((np.arange(1,len(ra_2)+1,1),idx+1,judge))

    np.savetxt(path+'match_old_newfits.txt',match_result,fmt='%10d %10d %10d')
    return idx[sep_constraint]

if __name__=='__main__':
    match_twofits()
    # (ra_x,dec_x,ra_IR,dec_IR)=load_data()
    # c1 = SkyCoord(ra=ra_x * u.degree, dec=dec_x * u.degree)
    # c2 = SkyCoord(ra=ra_IR * u.degree, dec=dec_IR * u.degree)
    # plt.scatter(ra_x,dec_x)
    # plt.scatter(ra_IR,dec_IR)
    #
    # idx, d2d, d3d = match_coordinates_sky(matchcoord=c1, catalogcoord=c2, nthneighbor=1)
    # max_sep = 0.1 * u.arcsec
    # sep_constraint = d2d < max_sep
    # c_matches = c1[sep_constraint]
    # catalog_matches = c2[idx[sep_constraint]]
    # d2d_matches = d2d[sep_constraint]
    # print(idx[sep_constraint])
    #
    # plt.scatter(ra_IR[idx[sep_constraint]],dec_IR[idx[sep_constraint]],marker='+',s=20)
    # print(c_matches)
    # print(catalog_matches)
    # print(d2d_matches.arcsec)

    plt.show()