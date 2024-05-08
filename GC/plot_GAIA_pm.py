'''
Author: baotong && baotong@smail.nju.edu.cn
Date: 2023-10-18 12:12:31
LastEditors: baotong && baotong@smail.nju.edu.cn
LastEditTime: 2023-10-25 16:57:30
FilePath: /pulsar/GC/plot_GAIA_pm.py
Description: 

Copyright (c) 2023 by baotong, All Rights Reserved. 
'''
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import pandas as pd
from astropy import units as u
from astropy import constants as c
from astropy.coordinates import SkyCoord,match_coordinates_sky
from astroquery.gaia import Gaia
from astropy.coordinates import Angle
from matplotlib.patches import Circle

font1 = {'family': 'Normal',
         'weight': 'normal',
         'size': 18, }

def masyr2kms(masyr,dist):
    ## masyr in the unit of mas/year
    ## dist in the unit of kpc
    a=masyr*dist*u.au/u.year
    b=a.to(u.km/u.s)
    return b.value

def plot_oneGC(GCname,GCpmra,GCpmdec,save=0,show=1):
    path='/Users/baotong/Desktop/period_terzan5/CSC/'
    cat=pd.read_csv(path+f'GAIA_{GCname}.csv')
    pmra=np.array(cat['pmra'])
    pmdec=np.array(cat['pmdec'])
    pmra = pmra[np.logical_not(np.isnan(pmra))]
    pmdec = pmdec[np.logical_not(np.isnan(pmdec))]
    dist=5.50;vesc=49.5
    xall=pmra-GCpmra;yall=pmdec-GCpmdec
    x1=1.04355343-GCpmra;y1=-6.896789407-GCpmdec
    x2=-5.300021308-GCpmra;y2=-8.205846514-GCpmdec
    x3=0.120670036-GCpmra;y3=-10.16596735-GCpmdec
    x4=-19.3736737-GCpmra;y4=-26.86913636-GCpmdec
    xall,yall=masyr2kms(xall,dist),masyr2kms(yall,dist)
    x1,y1=masyr2kms(x1,dist),masyr2kms(y1,dist)
    x2,y2=masyr2kms(x2,dist),masyr2kms(y2,dist)
    x3,y3=masyr2kms(x3,dist),masyr2kms(y3,dist)
    x4,y4=masyr2kms(x4,dist),masyr2kms(y4,dist)
    fig, ax = plt.subplots(figsize=(8,8))
    circle1 = Circle((0., 0.), vesc, fill=False, color='blue',linewidth=2, linestyle='--',label=r'$\rm V_{esc}$')
    circle2 = Circle((0., 0.), vesc*2, fill=False, color='cyan',linewidth=2, linestyle='--',label=r'$\rm 2\times V_{esc}$')
    plt.title(f'{GCname} (Method 1)',font1)
    plt.scatter(xall,yall,s=5,color='red',marker='.')
    plt.scatter([x1,x2,x3,x4],[y1,y2,y3,y4],s=100,color='green',marker='o')
    plt.text(x1,y1,'Seq.2',font1)
    plt.text(x2,y2,'Seq.3',font1)
    plt.text(x3,y3,'Seq.6',font1)
    plt.text(x4,y4,'Seq.7',font1)
    ax.add_patch(circle1);ax.add_patch(circle2)
    plt.xlim(-500,500)
    plt.ylim(-500,500)
    plt.xlabel('Proper motion in RA (km/s)',font1)
    plt.ylabel('Proper motion in DEC (km/s)',font1)
    plt.tick_params(labelsize=16)
    plt.legend()
    if save:
        # plt.savefig(path+f'{GCname}_GAIA_pm.png',bbox_inches='tight', pad_inches=0.05)
        # plt.savefig(path+f'{GCname}_GAIA_pm.pdf',bbox_inches='tight', pad_inches=0.05)
        plt.savefig(path+f'{GCname}_GAIA_pm_GCpos_3as.png',bbox_inches='tight', pad_inches=0.05)
        plt.savefig(path+f'{GCname}_GAIA_pm_GCpos_3as.pdf',bbox_inches='tight', pad_inches=0.05)
    if show:plt.show()

if __name__=='__main__':
    plot_oneGC(GCname='M28',GCpmra=-0.301,GCpmdec=-8.913,show=1,save=1)


