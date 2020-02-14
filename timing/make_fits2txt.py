import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import string

path='/Users/baotong/Desktop/period_LW/'
fitsname='LimWin_p50_2000_8000_src.fits'
catalog_LW=fits.open(path+fitsname)[1].data
seq=np.linspace(1,847,847)
ra=catalog_LW.field(0)
dec=catalog_LW.field(1)
phy_x=catalog_LW.field(4)
phy_y=catalog_LW.field(5)
net_cts=catalog_LW.field(9)
net_cts_err=catalog_LW.field(10)
bkg_cts=catalog_LW.field(11)
bkg_cts_err=catalog_LW.field(12)
net_rate=catalog_LW.field(13)
net_rate_err=catalog_LW.field(14)
bkg_rate=catalog_LW.field(15)
bkg_rate_err=catalog_LW.field(16)


catalog_info=np.column_stack((seq,ra,dec,phy_x,phy_y,net_cts,net_cts_err,bkg_cts,bkg_cts_err,net_rate,net_rate_err,bkg_rate,bkg_rate_err))
np.savetxt(path+'catalog_LW.txt',catalog_info,
           fmt='%5d %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %e %e %e  %e')