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
import random
import pandas as pd
from tkinter import _flatten
from astropy.wcs import WCS
import warnings

obs_id=[581,1431,441,582,2406,2405,2312,1672,2409,2313,2239,8591,9593,9718,8593,8597,8595,8592,8596,
9575,9578,8594,9596,12043,12123,12044,12128,12045,12129,12135,12046,12047,12137,12138,12055,12213,12048,
12049,12050,12222,12219,12051,12218,12223,12052,12220,12053,12054,12230,12231,12227,12233,12232,12234,16183,
16180,16456,16641,16457,16644,16463,17417,17416,16454,16176,16175,16178,16177,16620,16462,17535,17542,16184,
16182,16181,17546,16186,16187,16188,16450,16190,16189,17556,16179,17573,17633,17634,16453,16451,16461,16191,
16460,16459,17552,16455,16458,17677,18709,18719,16452,18730,16185]
path_out='/Users/baotong/Desktop/CDFS/'
TSTART=[]
TSTOP=[]
for id in obs_id:
    evtfits='/Volumes/pulsar/CDFS/merge_data/xdata/all_bcc_{0}_reproj_evt.fits'.format(str(id))
    hdul=fits.open(evtfits)
    TSTART.append(hdul[1].header['TSTART'])
    TSTOP.append(hdul[1].header['TSTOP'])
TSTART=np.array(TSTART)
TSTOP=np.array(TSTOP)
obs_id=np.array(obs_id)
T_exp=TSTOP-TSTART

epoch=np.column_stack((TSTART,TSTOP,obs_id,T_exp))
epoch = epoch[epoch[:,0].argsort()]
np.savetxt(path_out+'CDFS_epoch.txt',epoch,fmt='%15.2f %15.2f %10d %20.2f')