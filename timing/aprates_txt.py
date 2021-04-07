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


expmap_name = 'reproj_evt2_sou_' + str(int(obs)) + '_i4.fits'
reg = read_region(str(int(ID)) + '.reg')
print(src_cts, bkg_cts)
expmap = fits.open(expmap_name)
exptime_all = expmap[0].data
exptime_all = exptime_all.T
img_x = reg[0]
img_y = reg[1]
##特别注意，这里的expmap，physical坐标与image坐标一致(忘了是为啥了)
r = reg[2] ** 0.707
##为了简便，取个圆的内接正方形对exposure map平均吧
exp_src = exptime_all[np.arange(int(img_x - r), int(img_x + r))]
exp_src = exp_src[:, np.arange(int(img_y - r), int(img_y + r))]
exp_s = np.mean(exp_src)
aprates_text = 'aprates n={0} m={1} A_s={2} A_b={3} alpha=0.9 beta=0.02 T_s=1 ' \
               'E_s={4} eng_s=1 flux_s=1 T_b=1 E_b={5} eng_b=1 flux_b=1 clobber=yes ' \
               'outfile={6} conf=0.9973'.format(src_cts, bkg_cts, backscale, b_backscale, exp_s, exp_s,
                                                str(int(ID)) + '_' + str(int(obs)) + '_out.par')
with open('aprates/' + 'run_' + str(int(ID)) + '_' + str(int(obs)) + '.e', 'w+') as f:
    f.writelines(aprates_text)