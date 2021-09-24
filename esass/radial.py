import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
import pandas as pd
import sys
import os
from tkinter import _flatten
import funcs_timing as funcs
from scipy import interpolate
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.stats import poisson_conf_interval

ra_center=6.022318
dec_center=-72.081443
def read_erosita_cat(filename):
    cat=pd.read_excel(filename,header=0)
    srcid=cat['NAME']
    srcid=np.array(srcid)
    ra_hms=cat['RA']
    dec_hms=cat['DEC']
    ra=[];dec=[]
    for i in range(len(dec_hms)):
        skycoord=ra_hms[i]+dec_hms[i]
        c = SkyCoord(skycoord, unit=(u.hourangle, u.deg))
        ra.append(c.ra.value)
        dec.append(c.dec.value)

    return (ra,dec,srcid)

def dist_center_twocat(catname1,catname2):
    (ra_eR,dec_eR,srcIDlist)=read_erosita_cat(catname1)
    cat2=fits.open(catname2)
    ra_chand=cat2[1].data['RAdeg']
    dec_chand=cat2[1].data['DEdeg']
    srcIDlist_chand=np.arange(1,len(ra_chand)+1,1)

    c1=SkyCoord(ra_eR*u.deg,dec_eR*u.deg,frame='fk5')
    c2=SkyCoord(ra_chand*u.deg,dec_chand*u.deg,frame='fk5')
    cc=SkyCoord(ra_center*u.deg,dec_center*u.deg,frame='fk5')
    dist1=c1.separation(cc)
    dist2=c2.separation(cc)
    return (dist1,dist2)
if __name__=='__main__':
    # run_many_sim()
    path='/Users/baotong/Desktop/period_Tuc/'
    catname1=path+'erosita_cat_coord.xlsx'
    catname2=path+'xray_properties-592.fits'
    (dist1,dist2)=dist_center_twocat(catname1,catname2)
    dist1=dist1.arcmin;dist2=dist2.arcmin
    bins1=np.arange(int(np.min(dist1)),int(np.max(dist1))+2,1)
    bins2=np.arange(int(np.min(dist2)),int(np.max(dist2))+2,0.5)

    num1=plt.hist(dist1,bins=bins1,histtype='step')[0]
    num2=plt.hist(dist2,bins=bins2,histtype='step')[0]
    plt.close()

    num1_err=np.array(poisson_conf_interval(num1,interval='frequentist-confidence'))
    num1_err[0]=num1-num1_err[0]
    num1_err[1]=num1_err[1]-num1

    num2_err=np.array(poisson_conf_interval(num2,interval='frequentist-confidence'))
    num2_err[0]=num2-num2_err[0]
    num2_err[1]=num2_err[1]-num2
    print(num2_err[0])

    for i in range(len(num1)):
        num1[i]/=np.pi*(bins1[i+1]**2-bins1[i]**2)
        num1_err[0][i] /= np.pi * (bins1[i + 1] ** 2 - bins1[i] ** 2)
        num1_err[1][i] /= np.pi * (bins1[i + 1] ** 2 - bins1[i] ** 2)

    for i in range(len(num2)):
        num2[i]/=np.pi*(bins2[i+1]**2-bins2[i]**2)
        num2_err[0][i] /= np.pi * (bins2[i + 1] ** 2 - bins2[i] ** 2)
        num2_err[1][i] /= np.pi * (bins2[i + 1] ** 2 - bins2[i] ** 2)


    plt.errorbar(x=bins1[:-1]+0.5,y=num1,yerr=num1_err,xerr=np.zeros(len(num1))+0.5,marker='.',linestyle='')
    plt.errorbar(x=bins2[:-1]+0.25,y=num2,yerr=num2_err,xerr=np.zeros(len(num2))+0.25,marker='.',linestyle='')
    plt.xlabel('R (arcmin)',funcs.font2)
    plt.ylabel('Number of source per arcmin^2',funcs.font2)
    plt.tick_params(labelsize=16)
    plt.semilogy()
    plt.show()