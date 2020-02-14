import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import string

path='/Volumes/pulsar/SgrA/3549/repro/'
filename='fracarea_test.fits'
hdul=fits.open(path+filename)
time=hdul[1].data.field(0)
time=time-time[0]
frac=hdul[1].data.field(4)
plt.plot(time,frac)
plt.show()
