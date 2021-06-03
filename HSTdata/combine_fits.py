#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import string
#import correct as correct
from datetime import datetime
import pandas as pd
from astropy.stats import poisson_conf_interval
import scipy
from scipy.interpolate import lagrange
from scipy import interpolate
from scipy.optimize import curve_fit
from scipy.special import comb, perm
from astropy.table import Table, hstack

def combine_drc():
    path='/Users/baotong/Desktop/HST/anonymous47355/rawdata/07/'
    os.chdir(path)
    f1=fits.open('F547N_drc_ctx.fits')
    f2=fits.open('F547N_drc_sci.fits')
    f3=fits.open('F547N_drc_wht.fits')
    f_org=fits.open('ibir07020_drc.fits')
    hdul1=fits.ImageHDU(data=f2[0].data,header=f2[0].header,name='SCI')
    hdul2 = fits.ImageHDU(data=f3[0].data,header=f3[0].header,name='WHT')
    hdul3=fits.ImageHDU(data=f1[0].data,header=f1[0].header,name='CTX')
    hdul4=f2[1]
    hdul0=f_org[0]

    combine=fits.HDUList([hdul0,hdul1,hdul2,hdul3,hdul4])
    combine.writeto('combine_drc.fits',overwrite=True)

combine_drc()

def delete_rows():
    path='/Volumes/pulsar/xmmCDFS/merge_evt/'
    os.chdir(path)
    f1=fits.open('ep5_temp21_pn_filt_time.fits')
