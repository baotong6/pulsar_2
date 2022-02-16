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
from astropy.stats import poisson_conf_interval
import scipy
import hawkeye as hawk

def load_4longobs(dataname,ecf):
    for obsid in [700011,700163,700013,700014]:
        path=f'/Users/baotong/eSASS/data/raw_data/47_Tuc/txt/txt_psf{ecf}_{obsid}/'

