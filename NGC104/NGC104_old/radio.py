import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord,match_coordinates_sky
import pandas as pd
import sys
import os

path='/Users/baotong/Desktop/radio/'
wcen_list=pd.read_pickle(path+'47Tuc.ATCA.combine.5500.pickle')
print(wcen_list.columns)
print(wcen_list)
ra=wcen_list['ra'];dec=wcen_list['dec']
src_radius=np.zeros(len(ra))+5

print(wcen_list['peak_flux'])
for i in range(len(ra)):
    reg = 'circle(' + str(ra[i]) + ',' + str(dec[i]) + ',' + str(src_radius[i]) +'"'+ ')'
    with open(path+'47Tuc_radio.reg', 'a+') as f2:
        f2.writelines(reg + '\n')