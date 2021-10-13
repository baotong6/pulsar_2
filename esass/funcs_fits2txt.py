import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import pandas as pd
import sys
import os
from tkinter import _flatten
import math
from dataclasses import dataclass
from typing import Tuple
from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import SkyCoord

@dataclass
class Circle:
    x: float
    y: float
    r: float

    @property
    def coord(self):
        return self.x, self.y

def read_region(regname):
    reg_file=[]
    with open(regname, 'r') as file_to_read:
        while True:
            lines = file_to_read.readline()  # 整行读取数据
            reg_file.append(lines)
            if not lines:
                break
                pass
    region = reg_file[-2][7:-2]
    reg_x, reg_y, reg_r = [float(i) for i in region.split(',')]
    return [reg_x, reg_y, reg_r]

def where_region(x,y,reg):
    dist=np.sqrt((x-reg[0])**2+(y-reg[1])**2)
    temp=dist-reg[2]
    return np.where(temp<=0)

def where_region_annulus(x, y, reg):
    an=reg
    dist = np.sqrt((x-reg[0])**2+(y-reg[1])**2)
    temp_in= dist - an[2]
    temp_out=dist-an[3]
    return np.where((temp_in >= 0)&(temp_out <= 0))

def delete_photon_ID(time,energy,ID,emin=500,emax=8000):
    i=0
    while i < len(energy):
        if energy[i]>emax or energy[i]<emin:
        #if energy[i] > 8000 or energy[i] < 2000:
            energy=np.delete(energy,i)
            time=np.delete(time,i)
            ID=np.delete(ID,i)
            i=i-1
        i=i+1
    return [time,energy,ID]

def get_evt_srcreg(evtname,obsid,src_reg,emin,emax):
    evtall = fits.open(evtname)[1]
    time_all = evtall.data['TIME']
    energy_all = evtall.data['PI']
    RA_all = evtall.data['RA']
    DEC_all = evtall.data['DEC']
    X_all=evtall.data['X']
    Y_all=evtall.data['Y']

    index_out = where_region(X_all, Y_all, src_reg)
    time = time_all[index_out];
    energy = energy_all[index_out]
    obsID = np.array([obsid for i in range(len(time))])
    [src_t, src_E, src_ID] = delete_photon_ID(time, energy, obsID, emin=emin, emax=emax)
    return [src_t, src_E, src_ID]

def trans_radec2xy(imagename,ra,dec):
    image_file=fits.open(imagename)[0]
    w = WCS(image_file)
    src_x, src_y = w.all_world2pix(ra, dec, 1)
    phy_x=src_x*80-41000;phy_y=src_y*80-41000
    return (phy_x,phy_y)

def get_evt_bkgreg(evtname,imgname,obsid,bkg_reg,ra_list,dec_list,radius_list,emin,emax):
    ##reg must be annulus##
    evtall = fits.open(evtname)[1]
    time_all = evtall.data['TIME']
    energy_all = evtall.data['PI']
    RA_all = evtall.data['RA']
    DEC_all = evtall.data['DEC']
    X_all=evtall.data['X']
    Y_all=evtall.data['Y']
    index_out = where_region_annulus(X_all, Y_all, bkg_reg)

    X_evt=X_all[index_out]
    Y_evt=Y_all[index_out]
    time = time_all[index_out];
    energy = energy_all[index_out]
    ##去掉其他源的overlap##
    del_index=[];overlap_area=0
    (phy_x,phy_y)=trans_radec2xy(imgname,ra_list,dec_list)

    cir1=Circle(bkg_reg[0],bkg_reg[1],bkg_reg[2])
    cir2=Circle(bkg_reg[0],bkg_reg[1],bkg_reg[3])

    for i in range(len(radius_list)):
        dist=np.sqrt((X_evt-phy_x[i])**2+(Y_evt-phy_y[i])**2)

        del_index_single = np.where(dist < 2 * radius_list[i])[0]
        del_index = np.union1d(del_index, del_index_single)
        if len(del_index_single)>0:
            cir3=Circle(phy_x[i],phy_y[i],radius_list[i])
            overlap_area_single=find_intersection(cir2,cir3)-find_intersection(cir1,cir3)
        else:
            overlap_area_single=0
        overlap_area+=overlap_area_single

    del_index = del_index.astype('int64')
    X_evt = np.delete(X_evt, del_index);
    Y_evt = np.delete(Y_evt, del_index)
    time=np.delete(time,del_index)
    energy=np.delete(energy,del_index)
    obsID = np.array([obsid for i in range(len(time))])

    bkg_area=np.pi*(bkg_reg[3]**2-bkg_reg[2]**2)-overlap_area

    [bkg_t, bkg_E, bkg_ID] = delete_photon_ID(time, energy, obsID, emin=emin, emax=emax)

    return [bkg_t, bkg_E, bkg_ID,bkg_area]

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

def find_intersection(c1: Circle, c2: Circle)-> float:
    """Finds intersection area of two circles.

    Returns intersection area of two circles otherwise 0
    """
    d = math.dist(c1.coord, c2.coord)
    rad1sqr = c1.r ** 2
    rad2sqr = c2.r ** 2

    if d == 0:
        # the circle centers are the same
        return math.pi * min(c1.r, c2.r) ** 2

    angle1 = (rad1sqr + d ** 2 - rad2sqr) / (2 * c1.r * d)
    angle2 = (rad2sqr + d ** 2 - rad1sqr) / (2 * c2.r * d)

    # check if the circles are overlapping
    if (-1 <= angle1 < 1) or (-1 <= angle2 < 1):
        theta1 = math.acos(angle1) * 2
        theta2 = math.acos(angle2) * 2

        area1 = (0.5 * theta2 * rad2sqr) - (0.5 * rad2sqr * math.sin(theta2))
        area2 = (0.5 * theta1 * rad1sqr) - (0.5 * rad1sqr * math.sin(theta1))

        return area1 + area2
    elif angle1 < -1 or angle2 < -1:
        # Smaller circle is completely inside the largest circle.
        # Intersection area will be area of smaller circle
        # return area(c1_r), area(c2_r)
        return math.pi * min(c1.r, c2.r) ** 2
    return 0

def find_last_negative(arr):
    if (np.abs(arr)+arr).any()==0:
        print('All negative!')
        return len(arr)-1
    if (np.abs(arr) - arr).any() == 0:
        print('All positive!')
        return 0
    if arr.all()==np.sort(arr).all():
        i=0
        while i <len(arr):
            if arr[i]<=0:
                i+=1
            else:
                return i-1

    else:
        print('Not sorted!')
        return None