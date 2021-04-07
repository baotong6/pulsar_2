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
import read_csv as data
import pandas as pd
path='/Users/baotong/Desktop/period_LW/'

ra=data.ra_LW
dec=data.dec_LW
label=data.label_LW
with open(path+'all_pCV.reg','a+') as f2:
    f2.writelines('fk5'+'\n')
for i in range(len(ra)):
    if i==13 or i==21:
        continue
    elif label[i]==0:
        with open(path+'all_pCV.reg','a+') as f2:
            reg = 'circle(' + str(ra[i]) + ',' + str(dec[i]) + ',' + str('10"') + ')'+" # color=red width=2 text={"+str(i+1)+"}"
            f2.writelines(reg+'\n')

    elif label[i]!=0:
        with open(path+'all_pCV.reg','a+') as f2:
            reg = 'circle(' + str(ra[i]) + ',' + str(dec[i]) + ',' + str('10"') + ')'+" # color=magenta width=2 text={"+str(i+1)+"}"
            f2.writelines(reg+'\n')