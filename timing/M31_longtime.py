import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import lagrange
from scipy import interpolate
from astropy.io import fits
from scipy import optimize as op
import hawkeye as hawk


path='/Users/baotong/Desktop/command/zjc/'
catalog=fits.open(path+'M31acis_catalog_final.fits')
data=catalog[1].data
BCOUNT_TOT=np.array(data['BCOUNT_TOT'][0])
SCOUNT_TOT=np.array(data['SCOUNT_TOT'][0])
HCOUNT_TOT=np.array(data['HCOUNT_TOT'][0])
VCOUNT_TOT=np.array(data['VCOUNT_TOT'][0])

id=np.arange(1,len(BCOUNT_TOT)+1,1)
print(id)
plt.plot(id,BCOUNT_TOT)
plt.show()