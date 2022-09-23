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
src_ID = ['1180', '1182', '1628', '1538', '2961', '1487', '1514']
# src_ID = ['2560', '1502', '2532', '1206', '3067', '1624', '2508',
#           '3357', '2841', '1219', '2672', '2422', '1853', '3120', '1133', '2730',
#           '1084', '2525', '2157', '2187', '2344', '2199', '1677', '1634', '973']
#src_ID=['1671','2148','2338']
#path="/Volumes/pulsar/SgrA/merge_data/spectra/spectra_p/"
path="/Volumes/pulsar/SgrAGRT/merge_data/spectra/"
os.chdir(path)
for ID in src_ID:
    cmd="ls "+str(ID)+"*pi " +"> "+str(ID)+"_spectra.txt"
    print(cmd)
    os.system(cmd)
    cmd="combine_spectra @"+str(ID)+"_spectra.txt " +str(ID)+"_stack method=sum verbose=2 clobber=yes"
    print(cmd)
    os.system(cmd)
    cmd="cp "+str(ID)+"*stack* "+"./spectra_stack_src"
    print(cmd)
    os.system(cmd)

