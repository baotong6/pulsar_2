#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import pylab as pl
import string
from datetime import datetime
from scipy.interpolate import lagrange
from scipy import interpolate
from scipy.fftpack import fft,ifft
import scipy.signal as ss
import random
import pandas as pd
import read_csv as data

def getlen(a,b):
    length=((a[0]-b[0])**2+(a[1]-b[1])**2)**0.5
    return length
def compare_counterpart():
    path='/Users/baotong/Desktop/period_LW/'
    cat1=fits.open(path+'Maureen.fit')

    x2=data.ra_LW
    y2=data.dec_LW
    n2=data.ID_LW

    offset=5*1./3600.

    x1=cat1[1].data['_RAJ2000']
    y1=cat1[1].data['_DEJ2000']
    n1=cat1[1].data['LW']
    # ID_IR=np.array([311,350,421,461,505,538,599,632,
    #                 694,782,974,997,1000,1039,1053,1151,1277,1339,1399,1495,1493,1577,1667])
    #
    # x1=x1[ID_IR-1]
    # y1=y1[ID_IR-1]
    # n1=n1[ID_IR-1]

    match=[[] for i in range(len(n2))]
    for i in range(len(n2)):
        for j in range(len(n1)):
            if getlen([x2[i],y2[i]],[x1[j],y1[j]])<offset:
                match[i].append(n1[j])

    # print(x1)
    # print(y2)
    print(match)
    plt.scatter(x1,y1,color='red')
    plt.scatter(x2,y2,color='green')
    plt.show()
    plt.legend(['IR','X-ray'])
    # print(match)
compare_counterpart()