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
import read_csv as data
type=['NSC','ND','LW']
mode='NSC'
if mode=='ND':
    path='/Volumes/pulsar/GC/spectra/aprates'
    obs_ID=np.loadtxt('ACIS-I_epoch.txt')[:,2]
    os.chdir(path)
    src_ID=data.ID_ND
    for ID in src_ID:
        for obs in obs_ID:
            os.system('./run_{0}_{1}.e'.format(str(int(ID)),str(int(obs))))
if mode=='LW':
    path = '/Volumes/pulsar/LimWin_damage/merge_data/spectra/aprates'
    obs_ID=np.loadtxt('LW_epoch.txt')[:,2]
    os.chdir(path)
    src_ID=data.ID_LW
    for ID in src_ID:
        if str(int(ID))[-3:]=='001' or str(int(ID))[-3:]=='002':
            ID=str(ID)[:-3]
        for obs in obs_ID:
            os.system('./run_{0}_{1}.e'.format(str(int(ID)),str(int(obs))))
if mode=='NSC':
    path = '/Volumes/pulsar/SgrA/merge_data/spectra/spectra_p/aprates'
    obs_ID=np.loadtxt('SgrA_I_epoch.txt')[:,2]
    os.chdir(path)
    src_ID=data.ID_NSC
    for ID in src_ID:
        if str(int(ID))[-3:]=='001' or str(int(ID))[-3:]=='002':
            ID=str(ID)[:-3]
        for obs in obs_ID:
            os.system('./run_{0}_{1}.e'.format(str(int(ID)),str(int(obs))))



