import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import pandas as pd
import sys
import os

path='/Volumes/pulsar/47Tuc/merge_data/timing/'
obsid=[78,953,954,955,956,2735,3384,2736,3385,2737,3386,2738,3387,16527,15747,16529,17420,15748,16528]
for id in obsid:
    evtfile=fits.open(path+'all_bcc_{0}_reproj_evt.fits'.format(id))
    detname=evtfile[1].header['DETNAM']
    print(id,detname)