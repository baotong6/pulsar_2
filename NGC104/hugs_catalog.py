import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import pandas as pd
import sys
import os
from scipy import interpolate
from astropy.wcs import WCS
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord,match_coordinates_sky
import pandas as pd

def read_uv_cat(filename):
    cat_dat=np.loadtxt(filename)
    ra=cat_dat[:,33];dec=cat_dat[:,34]
    F275W=cat_dat[:,2];F336W=cat_dat[:,8]
    F435W=cat_dat[:,14];F606W=cat_dat[:,20]
    F814W=cat_dat[:,26]
    return (ra,dec,F275W,F336W,F435W,F606W,F814W)

def make_region_dot(ra,dec,outpath,outname,srcid=None):
    if not srcid: srcid = np.arange(1, len(ra) + 1, 1)
    for i in range(len(srcid)):
        # with open(outpath + '{0}.reg'.format(srcid[i]), 'w+') as f1:
        #     f1.writelines('fk5' + '\n')
        #     reg = 'circle(' + str(ra) + ',' + str(dec) + ',' + str(psfradii) + '"' + ')'
        #     f1.writelines(reg)
        with open(outpath + '{0}.reg'.format(outname), 'a+') as f2:
            if i==0: f2.writelines('fk5'+'\n')
            # reg = 'point(' + str(ra[i]) + ',' + str(dec[i]) + ')' + '# point=cross  text=' + '{' + str(
            #     srcid[i]) + '}'
            reg = 'point(' + str(ra[i]) + ',' + str(dec[i]) + ')' + '# point=cross'
            f2.writelines(reg + '\n')
def func_ngc6752():
    path='/Users/baotong/Desktop/period_NGC6752/hugs_ngc6752/'
    (ra1, dec1, F275W1, F336W1, F435W1, F606W1, F814W1)=read_uv_cat(path+'hlsp_hugs_hst_wfc3-uvis-acs-wfc_ngc6752_multi_v1_catalog-meth1.txt')
    (ra2, dec2, F275W2, F336W2, F435W2, F606W2, F814W2)=read_uv_cat(path+'hlsp_hugs_hst_wfc3-uvis-acs-wfc_ngc6752_multi_v1_catalog-meth2.txt')
    (ra3, dec3, F275W3, F336W3, F435W3, F606W3, F814W3)=read_uv_cat(path+'hlsp_hugs_hst_wfc3-uvis-acs-wfc_ngc6752_multi_v1_catalog-meth3.txt')

    ra4=[287.71497,287.71425,287.66826,287.76237]
    dec4=[-59.98379,-59.98476,-59.97816,-59.99503]
    outpath='/Users/baotong/Desktop/period_NGC6752/'
    make_region_dot(ra1,dec1,outpath,outname='meth1')
    make_region_dot(ra2,dec2,outpath,outname='meth2')
    make_region_dot(ra3,dec3,outpath,outname='meth3')
    # plt.scatter(ra1,dec1,marker='o',label='meth1')
    # plt.scatter(ra2,dec2,marker='*',label='meth2')
    # plt.scatter(ra3,dec3,marker='^',label='meth3')
    # plt.scatter(ra4,dec4,marker='o',label='pCV')
    # plt.legend()
    # plt.show()

def func_m62():
    path='/Users/baotong/Desktop/period_NGC6266/'
    catname='ngc6266_2002.xlsx'
    catinfo=pd.read_excel(path+catname)
    print(catinfo)

func_m62()