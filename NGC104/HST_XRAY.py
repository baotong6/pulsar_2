import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord,match_coordinates_sky
import pandas as pd
import sys
import os
from scipy import interpolate

path='/Users/baotong/Desktop/period_Tuc/'
res = fits.open(path + 'xray_properties-592.fits')[1].data
ra_xray = res['RAdeg'];dec_xray = res['DEdeg'];ID_xray = res['Seq']

res2=fits.open(path+'HST2001.fit')[1].data
ra_HST=res2['_RAJ2000'];dec_HST = res2['_DEJ2000'];Name = res2['VName']
period=res2['Per']

res3=fits.open(path+'Grindlay2001.fit')[1].data
ra_Grind = res3['_RAJ2000'];dec_Grind = res3['_DEJ2000'];Name_Grind=res3['__GHE2001_']

def get_chandra_HST2001():
    c_chandra= SkyCoord(ra=ra_xray*u.degree, dec=dec_xray*u.degree)
    c_HST=SkyCoord(ra=ra_HST*u.degree, dec=dec_HST*u.degree)

    idx, d2d, d3d=match_coordinates_sky(matchcoord=c_HST, catalogcoord=c_chandra, nthneighbor=1)
    max_sep = 5.0 * u.arcsec
    sep_constraint = d2d < max_sep
    c_matches = c_HST[sep_constraint]
    catalog_matches = c_chandra[idx[sep_constraint]]
    d2d_matches=d2d[sep_constraint]


    print(Name[sep_constraint])
    print(ID_xray[idx[sep_constraint]])
    print(period[sep_constraint]*86400)
    # print(d2d_matches.arcsec)
    a=[245,453,245,182,224,206,453,223,211,462,402,294,364,304,485,182,258,283,345]
    print(np.intersect1d(a,ID_xray[idx[sep_constraint]]))
    print(period[sep_constraint][np.where(ID_xray[idx[sep_constraint]]==345)]*86400)
    print(d2d_matches.arcsec[np.where(ID_xray[idx[sep_constraint]]==345)])
    print(Name[sep_constraint][np.where(ID_xray[idx[sep_constraint]]==345)])

    return None
# get_chandra_HST2001()

def get_chandra_2003_now():
    c_chandra= SkyCoord(ra=ra_xray*u.degree, dec=dec_xray*u.degree)
    c_Grind=SkyCoord(ra=ra_Grind*u.degree, dec=dec_Grind*u.degree)

    idx, d2d, d3d=match_coordinates_sky(matchcoord=c_Grind, catalogcoord=c_chandra, nthneighbor=1)

    a=[245,453,245,182,224,206,453,223,211,462,402,294,364,304,485,182,258,283,345]

    a_opt03=np.array([12,14,15,18,26,27,29,30,34,36,38,41,42,44,47,72,73,75,1,2,3,8,11,23,
             33,45,69,70,76,92,66,9,21,22,25,59,68,94])
    a_opt03_extend=[167,163,182,120,121,197,184]
    # print(ID_xray[idx])
    print(d2d.arcsec)
    print(np.intersect1d(a,ID_xray[idx]))
    print(ID_xray[idx[a_opt03-1]])

    # max_sep = 3.0 * u.arcsec
    # sep_constraint = d2d < max_sep
    # c_matches = c_Grind[sep_constraint]
    # catalog_matches = c_chandra[idx[sep_constraint]]



    return None

get_chandra_2003_now()
