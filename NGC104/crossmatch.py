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
import esass.funcs_fits2txt as funcs
import latex.poisson_error as pos_err
ra_center=6.022318
dec_center=-72.081443

def load_data():
    ## for LW x-ray and IR
    path='/Users/baotong/Desktop/period_Tuc/'
    xray_info=pd.read_csv(path+'radio_cat/ngc104_5sig.csv')
    ra_x = np.array(xray_info['RA'])
    dec_x = np.array(xray_info['DEC'])
    seq_x = np.arange(1,len(ra_x)+1,1)

    IR_info=fits.open(path+'Cheng2019.fit')
    ra_IR=IR_info[1].data['RAJ2000']
    dec_IR=IR_info[1].data['DEJ2000']
    # catname1=path+'erosita_cat_coord.xlsx'
    # (ra_IR, dec_IR, srcIDlist) = funcs.read_erosita_cat(catname1)

    return (ra_x,dec_x,ra_IR,dec_IR)

def match_twofits():
    path='/Users/baotong/Desktop/period_Tuc/'


    file1=fits.open(path+'cheng2019.fit')
    ra_1=file1[1].data['RAJ2000']
    dec_1=file1[1].data['DEJ2000']

    file2=fits.open(path+'2019spectra.fit')
    ra_2=file2[1].data['RAJ2000']
    dec_2=file2[1].data['DEJ2000']

    c1 = SkyCoord(ra=ra_2 * u.degree, dec=dec_2 * u.degree)
    c2 = SkyCoord(ra=ra_1 * u.degree, dec=dec_1 * u.degree)
    idx, d2d, d3d = match_coordinates_sky(matchcoord=c1, catalogcoord=c2, nthneighbor=1)

    max_sep =1* u.arcsec
    sep_constraint = d2d < max_sep
    c_matches = c1[sep_constraint]
    catalog_matches = c2[idx[sep_constraint]]
    d2d_matches = d2d[sep_constraint]
    print(idx)
    judge=np.array(sep_constraint).astype('int')
    print(idx[sep_constraint])
    match_result=np.column_stack((np.arange(1,len(ra_2)+1,1),idx+1,judge))

    np.savetxt(path+'match_spectra.txt',match_result,fmt='%10d %10d %10d')
    return idx[sep_constraint]

if __name__=='__main__':
    (ra_x,dec_x,ra_IR,dec_IR)=load_data()
    print(len(ra_x))
    c1 = SkyCoord(ra=ra_x * u.degree, dec=dec_x * u.degree)
    c2 = SkyCoord(ra=ra_IR * u.degree, dec=dec_IR * u.degree)
    cc=SkyCoord(ra_center*u.deg,dec_center*u.deg,frame='fk5')
    dist1 = c1.separation(cc)
    dist=dist1.arcsec
    # plt.scatter(ra_x,dec_x)
    # plt.scatter(ra_IR,dec_IR)

    idx, d2d, d3d = match_coordinates_sky(matchcoord=c1, catalogcoord=c2, nthneighbor=1)
    print(idx)
    max_sep =1* u.arcsec
    sep_constraint = d2d < max_sep
    c_matches = c1[sep_constraint]
    catalog_matches = c2[idx[sep_constraint]]
    d2d_matches = d2d[sep_constraint]
    print(idx[sep_constraint])
    ra_IR=np.array(ra_IR);dec_IR=np.array(dec_IR)
    # plt.scatter(ra_IR[idx[sep_constraint]],dec_IR[idx[sep_constraint]],marker='+',s=20)
    print(c_matches)
    print(dist)

    # bins=np.linspace(0,500,8)
    # hist1=plt.hist(dist,bins=bins,histtype='step')
    # plt.close()
    # number=hist1[0]
    # area=[]
    # for i in range(len(hist1[1])-1):
    #     area.append(np.pi*(hist1[1][i+1]**2-hist1[1][i]**2))
    # area=np.array(area)
    # print(number)
    # x=bins[1:]-35.7
    # x_err=35.7
    # y=number/area
    # y_err=np.sqrt(number)/area
    # plt.errorbar(x,y,xerr=x_err,yerr=y_err, fmt='co', capsize=4, elinewidth=2, ecolor='red', color='green')
    # plt.xlabel('Arcsec')
    # plt.ylabel('Number of sources per unit arcsec^2')
    # print(catalog_matches)
    # print(d2d_matches.arcsec)
    #
    # plt.show()
    match_twofits()