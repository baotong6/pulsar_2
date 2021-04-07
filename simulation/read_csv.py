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
from scipy.optimize import curve_fit
import pandas as pd

path_table = '/Users/baotong/Desktop/period/table/'
result_NSC = pd.read_excel(path_table + 'final_all_del.csv', 'result_NSC')
result_LW=pd.read_excel(path_table+'final_all_del.csv','result_LW')
result_ND=pd.read_excel(path_table+'final_all_del.csv','result_ND')
result_NSC_IG=pd.read_excel(path_table + 'final_all_del.csv', 'result_NSC_IG')

ID_NSC_IG=result_NSC_IG['seq']
P_NSC_IG=result_NSC_IG['P']
ra_NSC_IG=result_NSC_IG['ra']
dec_NSC_IG=result_NSC_IG['dec']
net_percent_NSC_IG=result_NSC_IG['net_percent']

ID_NSC=result_NSC['seq']
P_NSC=result_NSC['P']
flux_NSC=result_NSC['flux']
label_NSC=result_NSC['label']
ra_NSC=result_NSC['ra']
dec_NSC=result_NSC['dec']
net_percent_NSC=result_NSC['net_percent']

ID_LW=result_LW['seq']
P_LW=result_LW['P']
flux_LW=result_LW['flux']
label_LW=result_LW['label']
ra_LW=result_LW['ra']
dec_LW=result_LW['dec']
net_percent_LW=result_LW['net_percent']

ID_ND=result_ND['seq']
P_ND=result_ND['P']
flux_ND=result_ND['flux']
ra_ND=result_ND['ra']
dec_ND=result_ND['dec']
net_percent_ND=result_ND['net_percent']