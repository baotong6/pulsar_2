#!/usr/bin/env python
# coding: utf-8
import numpy as np
import astropy

from astropy import coordinates

from astroquery.mast import Observations
from astropy.table import QTable
import matplotlib.pyplot as plt
import astropy
import pandas as pd

from stsci.tools import teal 
from wfc3tools import calwf3
import wfc3tools
from glob import glob
import os

path='/Users/baotong/Desktop/HST/anonymous47355/rawdata/01'
os.chdir(path)

# wfc3tools.wf3cte('')
# wfc3tools.wf3ccd('')
# wfc3tools.wf32d('')
# wfc3tools.calwf3('ibir07mkq_raw.fits')

# wfc3tools.calwf3('ibir03020_asn.fits')

import drizzlepac
from drizzlepac import astrodrizzle
from drizzlepac import tweakreg

teal.unlearn('tweakreg')
tweakreg.TweakReg('ibir01*flc.fits',updatehdr=False,updatewcs=True,conv_width=3.5,threshold=200,peakmax=50000)

teal.unlearn('astrodrizzle')

# astrodrizzle.AstroDrizzle('ibir01*flc.fits',output='F658N',driz_sep_bits='64,32',driz_cr_corr=True,final_bits='64,32')


# In[ ]:




