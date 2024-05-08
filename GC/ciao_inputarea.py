'''
Author: baotong && baotong@smail.nju.edu.cn
Date: 2023-05-16 15:18:31
LastEditors: baotong && baotong@smail.nju.edu.cn
LastEditTime: 2024-03-27 18:56:26
FilePath: /pulsar/GC/ciao_inputarea.py
Description: 

Copyright (c) 2024 by baotong, All Rights Reserved. 
'''
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
    for i in range(1,3620):
        print(i)
        cmd="dmstat "+ '"'+image+'[sky=region(region_3392/region_90/'+f'{i}.reg)]" '+ "centroid = no"
        os.system(cmd)
        cmd='good=`pget dmstat out_good`'
        os.system(cmd)
        srcarea = os.popen('pget dmstat out_good').read()
        cmd="dmstat "+ '"'+image+'[sky=region(region_3392/region_90/'+f'{i}_bkg.reg)]" '+ "centroid = no"
        os.system(cmd)
        cmd='good=`pget dmstat out_good`'
        os.system(cmd)
        bkgarea = os.popen('pget dmstat out_good').read()
        str_tmp="{0:10d} {1:10.5f} {2:10.5f}".format(int(i),float(srcarea)*0.492**2,float(bkgarea)*0.492**2)
        with open(path+'all_area.txt', 'a+') as f3:
            f3.writelines(str_tmp+'\n')

    return None

if __name__=='__main__':
    image='SgrA_5.fits'
    path='/Users/baotong/Desktop/period/backscale/'
    input_area(path=path,image=image,regname=None)


    