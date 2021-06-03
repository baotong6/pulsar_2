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
from scipy.interpolate import lagrange
from scipy import interpolate
from scipy.optimize import curve_fit
import pandas as pd
from astropy.stats import poisson_conf_interval
import scipy
from astropy.io import fits
from astropy.wcs import WCS
def two_fits_min(fits1,fits2):
    err=fits1[1].data-fits2[1].data
    plt.imshow(err)
    print(err)
    plt.show()

# path='/Users/baotong/Desktop/HST/anonymous47355/'
# fitsname='ibir02nfq_flc.fits'
# fits1=fits.open(path+fitsname)
# fits2=fits.open(path+'rawdata/02/'+fitsname)
# two_fits_min(fits1,fits2)

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import string
#import correct as correct
from datetime import datetime
from scipy.interpolate import lagrange
from scipy import interpolate
from scipy.optimize import curve_fit
import pandas as pd
from astropy.stats import poisson_conf_interval
import scipy
from astropy.io import fits
from astropy.wcs import WCS
imgfile='/Users/baotong/Downloads/F373N_drc.fits'
f=fits.open(imgfile)
w = WCS(f[1].header)

ra=10.6847447
dec=41.2631402

img_x, img_y = w.all_world2pix(ra, dec, 1)

x=2554.7211;y=2240.2517
ra,dec=w.all_pix2world(x,y)
print(ra)