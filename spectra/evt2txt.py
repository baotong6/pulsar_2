#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import string
from datetime import datetime
from scipy.interpolate import lagrange
from scipy import interpolate
from scipy.optimize import curve_fit
import pandas as pd

path='/Volumes/pulsar/LimWin_damage/merge_data/spectra/'
bkg_evt='513_bkg_evt.fits'
hdul_evt=fits.open(path+bkg_evt)
x=hdul_evt[1].data.field(10)
y=hdul_evt[1].data.field(11)
energy=hdul_evt[1].data.field(14)
time=hdul_evt[1].data.field(0)
obs_ID=[111 for i in range(len(x))]

src_txt=np.column_stack((time,energy,obs_ID))
np.savetxt('513_bkg.txt',src_txt,fmt="%.7f  %5.3f  %d")