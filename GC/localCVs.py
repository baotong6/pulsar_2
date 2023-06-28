#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import csv
from datetime import datetime
from scipy.interpolate import lagrange
from scipy import interpolate
from scipy.optimize import curve_fit
import pandas as pd
from astropy.stats import poisson_conf_interval
import hawkeye as hawk

mCVs={'CC Scl':[1.402,31.4],
         'EX Hya':[1.638,31.59],
         'AO Psc':[3.591,32.2],
         'DO Dra':[3.969,31.4],
         'TV Col':[5.486,32.22],
         'EI UMa':[6.434,33.2],
         'CV Hyi':[1.297,32.3],
         'V4738 Sgr':[1.300,31.8],
         'EV UMa':[1.328,33.3],
         'GG Leo':[1.331,31.5],
         'FH UMa':[1.336,32.0],
         'EF Eri':[1.350,32.0],
         'IW Eri':[1.452,31.6],
         'EU UMa':[1.502,32.4],
         'EQ Cet':[1.547,31.4],
      'V393 Pav': [1.647,32.4],
      'EG Lyn':[1.656,32.4],
      'RS Cae':[1.699,33.0],
      'CD Ind':[1.848,31.7],
      'BL Hyi':[1.894,31.7],
      'EK UMa':[1.909,32.6],
      'AN UMa':[1.914,31.9],
      'V1007 Her':[1.999,32.3],
      'HU Aqr': [2.084,32.2],
      'UW Pic': [2.223,31.9],
      'RX J0859': [2.397,32.6],
      'CW Hyi': [3.030,32.4],
      'RX J1610': [3.176,32.1],
      '1RXS J231603': [3.491,32.8],
      'AI Tri': [4.602,32.5]
      }
non_mCVs={'SW UMa':[1.364,30.6],
         'CC Scl':[1.41,31.4],
         'T Leo':[1.412,30.6],
         'BZ UMa':[1.632,31.1],
         'RX J1715':[1.64,30.4],
         'VW Hyi':[1.783,30.5],
         'WX Hyi':[1.795,31.4],
         'SU UMa':[1.832,31.9],
         'SDSS J1730':[1.84,31.2],
         'TY PsA':[2.018,31.3],
         'TT Ari':[3.301,31.7],
         'EF Tuc':[3.5,31.4],
         'RX J1831':[4.01,31.5],
         'WW Cet':[4.22,31.2],
      'V405 Peg': [4.264,30.5],
      'EX Dra':[5.039,29.7]
      # 'TW Pic':[1.699,33.0],
      # 'RBS490':[1.848,31.7],
      # 'RBS1411':[1.894,31.7],
      # 'IQ Eri':[1.909,32.6]
          }