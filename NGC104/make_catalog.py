#!/bin/bash
# -*- coding: utf-8 -*-
# written by Tong
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
from astropy.timeseries import LombScargle
import stingray
from stingray import Lightcurve, Powerspectrum, AveragedPowerspectrum
from stingray.events import EventList

gcname=['47Tuc','terzan5','omg_cen','M28','NGC6397','NGC6752','NGC6266','M30']
i=0
path='/Users/baotong/Desktop/period_'+gcname[i]
catalogname=[]