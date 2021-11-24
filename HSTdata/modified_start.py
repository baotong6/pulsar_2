#!/usr/bin/env python
# coding: utf-8
import numpy as np
import astropy

from astroquery.mast import Observations
from astropy.table import QTable
import matplotlib.pyplot as plt
import astropy
from stsci.tools import teal 
from wfc3tools import calwf3
import wfc3tools
from glob import glob
import os


path='/Users/baotong/Desktop/HST/mastDownload/HSt/calwww/'
os.chdir(path)

wfc3tools.wf3rej(path+'ibir07meq_raw.fits')
# teal.unlearn('tweakreg')
# tweakreg.TweakReg('ibir01*flc.fits',updatehdr=False,updatewcs=True,conv_width=3.5,threshold=200,peakmax=50000)
#
# teal.unlearn('astrodrizzle')

# astrodrizzle.AstroDrizzle('ibir01*flc.fits',output='F658N',driz_sep_bits='64,32',driz_cr_corr=True,final_bits='64,32')


# In[ ]:




