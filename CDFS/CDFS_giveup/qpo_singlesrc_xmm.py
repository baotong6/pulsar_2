#!/bin/bash
# -*- coding: utf-8 -*-
# written by Tong
# plot the phase-folded light curve from txt file (version for xmm)
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
from astropy.timeseries import LombScargle
import stingray
from stingray import Lightcurve, Powerspectrum, AveragedPowerspectrum
font1 = {'family': 'Normal',
         'weight': 'normal',
         'size': 16, }
