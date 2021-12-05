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
from scipy.fftpack import fft,ifft
import scipy.signal as ss
import random
from astropy.wcs import WCS
def read_region(regname):
    reg_file = []
    with open(regname, 'r') as file_to_read:
        while True:
            lines = file_to_read.readline()  # 整行读取数据
            reg_file.append(lines)
            if not lines:
                break
                pass
    region = reg_file[-2][7:-2]
    reg_x, reg_y, reg_r = [float(i) for i in region.split(',')]
    return [reg_x, reg_y, reg_r]
path='/Volumes/pulsar/47Tuc/merge_data/timing/region/'
srcid=np.arange(1,593,1)
print(srcid)
for i in range(len(srcid)):
    [reg_x, reg_y, reg_r]=read_region(path+str(srcid[i])+'.reg')
    with open('/Volumes/pulsar/47Tuc/merge_data/timing/region_fk5/{0}.reg'.format(srcid[i]),'w+') as f1:
        f1.writelines('fk5'+'\n')
        reg = 'circle(' + str(reg_x) + ',' + str(reg_y) + ',' + str(reg_r) +'"'+')'
        f1.writelines(reg)
    with open('/Volumes/pulsar/47Tuc/merge_data/timing/region_fk5/all.reg','a+') as f2:
        reg = 'circle(' + str(reg_x) + ',' + str(reg_y) + ',' + str(reg_r) + '"' + ')'+'# text='+'{'+str(srcid[i])+'}'
        f2.writelines(reg+'\n')
