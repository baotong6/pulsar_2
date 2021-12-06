#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import string

#path='/Volumes/halo/GC/merge_data/xdata/'
path='/Volumes/pulsar/47Tuc/merge_data/xdata/'
# path='/Volumes/pulsar/CDFS/merge_data/xdata/'
band=['4']
#band=['2','3','4']
for k in band:
    psfmap='expweighted_75_mean_i'+k+'.psfmap'
    hdul_psf=fits.open(path+psfmap)
    psf90=hdul_psf[0].data
    psf90=psf90.T
    for i in range(0,len(psf90)):
        for j in range(0,len(psf90[i])):
            if psf90[i][j]<100:
                kkk=0
                #print('sss')
            else:
                kkk=1
            if kkk==1:
                psf90[i][j]=0

    if os.path.exists(path + 'mod_'+psfmap):
        os.remove(path + 'mod_'+psfmap)
    psf90 = psf90.T
    grey = fits.PrimaryHDU(psf90)
    greyHDU = fits.HDUList([grey])
    greyHDU.writeto(path + 'mod_'+psfmap)
