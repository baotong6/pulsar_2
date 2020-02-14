import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import string
import math
from astropy.coordinates import FK5
from astropy.time import Time

# time=Time(2457061.5,format='jd',scale='utc')
# print
# print time.tt.delta_tdb_tt
path='/Users/baotong/chandra/input/CDFS/2239/repro/'
srclist='source_list.fits'
hdul_src=fits.open(path+srclist)
#print hdul.info()
x=hdul_src[1].data.field(4)
y=hdul_src[1].data.field(5)

ctsimg="acis_I_test_broad_thresh.img"
hdul_cts=fits.open(path+ctsimg)
pixel_data=hdul_cts[0].data
#index=np.argwhere(pixel_data>0)
def convert_phy_imgCod(x,y):
    min_X=hdul_src[1].header['TLMIN5']
    max_X=hdul_src[1].header['TLMAX5']
    min_Y=hdul_src[1].header['TLMIN6']
    max_Y=hdul_src[1].header['TLMAX6']

    min_pixel_x=0
    min_pixel_y=0
    max_pixel_x=len(pixel_data[0])
    max_pixel_y=len(pixel_data)

    img_x=(x-min_X)/4.0
    img_y=(y-min_Y)/4.0

    return (img_x,img_y)

(img_x,img_y)=convert_phy_imgCod(x,y)
img_x=np.rint(img_x)
img_y=np.rint(img_y)
mat_data=np.mat(pixel_data)
mat_data=np.transpose(mat_data)
counts=[]
for i in range(len(img_x)):
    counts.append(mat_data[int(img_x[i])-2:int(img_x[i])+2,int(img_y[i])-2:int(img_y[i])+2])

print(counts)