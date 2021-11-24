import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

path='/Users/baotong/Desktop/period_Tuc/'
a=fits.open(path+'omg_cen_p50_i5_src_1_2_4_8.fits')
try:
    b=a[1].header['fuck']
except:
    print('No fuck')
    b=a[1].header['TTYPE1']
    print(b)
