import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import string

obsID = ['5934', '6362', '6365', '9500', '9501', '9502',
         '9503', '9504', '9505', '9854', '9855', '9892', '9893']
srcid = '528'

for i in range(len(obsID)):
    id=obsID[i]
    path='/Volumes/pulsar/LimWin_damage/'
    regionfile=path+f'merge_data/timing/region_{id}/region_90/{srcid}.reg'
    asolname = 'aspect_1_bcc.fits'
    os.chdir(path + id + '/cal')
    cmd='dither_region infile={0} region="region({1})" maskfile=maskfile.fits ' \
            'wcsfile=evt2file_bcc.fits outfile=fracarea_vs_time_{2}_src{3}.fits clobber=yes verbose=5'.format(asolname,regionfile,id,srcid)
    print(cmd)
    os.system(cmd)
