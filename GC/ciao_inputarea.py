#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 19:33:00 2023
@author: baotong
"""
import numpy as np
import os
import glob
import subprocess
import csv
import pandas as pd
import math
from astropy.io import fits

def input_area(path,image,regname):
    os.chdir(path)
    for i in range(1,377):
        print(i)
        cmd="dmstat "+ '"'+image+'[sky=region(region_startover/region_7460/region_90/'+f'{i}.reg)]" '+ "centroid = no"
        os.system(cmd)
        cmd='good=`pget dmstat out_good`'
        os.system(cmd)
        srcarea = os.popen('pget dmstat out_good').read()
        cmd="dmstat "+ '"'+image+'[sky=region(region_startover/region_7460/region_90/'+f'{i}_bkg.reg)]" '+ "centroid = no"
        os.system(cmd)
        cmd='good=`pget dmstat out_good`'
        os.system(cmd)
        bkgarea = os.popen('pget dmstat out_good').read()
        str_tmp="{0:10d} {1:10.5f} {2:10.5f}".format(int(i),float(srcarea)*0.492**2,float(bkgarea)*0.492**2)
        with open(path+'backscale/' + 'all_area.txt', 'a+') as f3:
            f3.writelines(str_tmp+'\n')

    return None

if __name__=='__main__':
    image='NGC6397_5_exp.fits'
    path='/Users/baotong/Desktop/period_NGC6397/spectra_startover/'
    input_area(path='/Users/baotong/Desktop/period_NGC6397/',image=image,regname=None)