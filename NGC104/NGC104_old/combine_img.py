#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import string
from scipy.interpolate import lagrange
from scipy import interpolate
from scipy.optimize import curve_fit
import pandas as pd

def combine_image_fits(imgfile):
    img_out = fits.open(imgfile[0])[0].data
    for i in range(1,len(imgfile)):
        img_out+=fits.open(imgfile[i])[0].data

    return img_out
path1='/Users/baotong/xmm/0651790101/cal/'
path2='/Users/baotong/xmm/0802580201/cal/'
imgfilelist=[path1+'mos1_filt.img',path1+'mos2_filt.img',path1+'pn_filt.img',
             path2+'mos1_filt.img',path2+'mos2_filt.img',path2+'pn_filt.img']
img_out=combine_image_fits(imgfilelist)
print(img_out)

header = fits.getheader(imgfilelist[0], 0)
grey = fits.PrimaryHDU(img_out)
greyHDU = fits.HDUList([grey])
print(header)
fits.writeto(path1+'combine_img.fits', img_out, header, overwrite=True)
