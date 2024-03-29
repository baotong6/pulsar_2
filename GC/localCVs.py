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

GDNe={'BZ Uma':[1.63,31.17],
      'HT Cas':[1.77,30.79],
      'SS Aur':[4.39,30.98],
      'SW UMma':[1.36,30.69],
      'U Gem':[4.25,30.92],
      'T Leo':[1.42,30.81],
      'V893 Sco':[1.82,31.7],
      'SS Cyg':[6.603,31.888],
      'Z Cam':[6.98,31.79],
      'WZ Sge':[1.36,29.85],
      'VY Aqr':[1.51,30.11],
      'ASAS J0025':[1.37,30.2],
      'GW Lib': [1.28, 28.7]}

Asypolar={'BYCam':[3.354,32.813],
          'V1432Aql':[3.366,33.021]}

IPlow={'EXHya':[1.6376,31.3979],
     'DODra':[3.98898,32.0414],
     'V1025cen':[1.41024,31.8846],
     'GKPer':[47.9233,32.663],
     'CTCVJ2056-3014':[1.76,30.778],
     'FO Aqr':[4.8508,32.87]}

IPhigh={'TVCol':[5.4864,33.29],
        'V405Aur':[4.1429,33.38],
        'MUCam':[4.7245,33.40],
        'V667Pup':[5.6112,33.30],
        'V24000Oph':[3.408,33.53],
        'BGCMi':[3.234,33.64],
        'V709Cas':[5.33,33.54],
        'V2731Oph':[15.42,33.61],
        'PQGem':[5.193,33.39],
        'V1223Sgr':[3.366,33.86],
        'AOPsc':[3.59,33.308],
        'V1062Tau':[9.98,33.83],
        'FOAqr':[4.8508,33.42],
        'TXCol':[5.7192,32.94],
        'GKPer':[47.9233,34.29],
        'V2306Cyg':[4.35708,33.56],
        'NYLup':[9.864,34.13],
        'V1033Cas':[4.032,33.41],
        'RXJ2133':[7.14,33.99],
        'V418Gem':[4.3704,34.15],
        'V515And':[2.731,33.49],
        'V647Aur':[3.47,33.875],
        'EIUMa':[6.4344,33.57],
        'V2069Cyg':[7.48032,33.47],
        'IGRJ1719':[4.0056,33.37],
        'IGRJ0457':[7.2,34.64],
        'IGRJ1817':[1.5312,34.22],
        'IGRJ0838':[7.92,33.64],
        'IGRJ1509':[5.8892,33.55],
        'IGRJ1649':[3.6168,33.55],
        'IGRJ1654':[3.7152,33.54]}

Polarhigh={'1RXSJ165424':[2.87,31.70],
           'RXJ0154':[1.48,32.23],
           'RXJ1002':[1.67,32.73],
           'RXJ0600':[1.31,32.83],
           'RXJ0859':[2.40,33.09],
           'RXJ0953':[1.73,32.95],
           'CW Hyi': [3.030, 32.4],
           'RX J1610': [3.176, 32.1],
           '1RXS J231603': [3.491, 32.8],
           'AI Tri': [4.602, 32.5]
           }

Polarlow={'1RXSJ165424':[2.87,30.25],
          'RXJ0154':[1.48,31.65],
          'RXJ1002':[1.67,32.53],
          'V379Vir':[1.47,29.36],
          'SDSS1250':[1.44,29.06],
          'SDSS1514':[1.48,29.36]}

periodbouncer={'V379Vir':[1.47,29.36],
          'SDSS1250':[1.44,29.06],
          'SDSS1514':[1.48,29.36],
          'GWLib': [1.28, 28.77],
          'WZSge': [1.36, 29.92],
          'VYAqr': [1.51, 30.18],
          'ASASJ0025': [1.37, 30.27],
          'OVBoo': [1.11, 30.36],
          'SDSSJ1035': [1.37, 28.91]}

nonmCVs_2015={'BZUma':[1.63,31.24],
              'HTCas':[1.77,30.86],
              'SSAur':[4.39,31.054],
              'SWUma':[1.36,30.76],
              'UGem':[4.25,30.99],
              'TLeo':[1.42,30.88],
              'V893Sco':[1.82,31.77],
              'SSCyg':[6.603,31.96],
              'ZCam':[6.98,31.86],
              '1RXHJ0826':[10.4,32.27],
              'EYCyg':[11.02,32.01],
              'CXOGBSJ17':[10.3,32.35],
              'GWLib':[1.28,28.77],
              'WZSge':[1.36,29.92],
              'VYAqr':[1.51,30.18],
              'ASASJ0025':[1.37,30.27],
              'OVBoo':[1.11,30.36],
              'SDSSJ1035':[1.37,28.91],
              'SWUMa':[1.364,31.1],
              'CCScl':[1.41,31.4],
              'TLeo':[1.412,31.3],
              'BZUMa':[1.632,31.3],
              'RXJ1715':[1.64,31.4],
              'VWHyi':[1.783,30.8],
              'WXHyi':[1.795,31.7],
              'SUUMa':[1.832,32.1],
              'SDSSJ1730':[1.84,31.9],
              'TYPsA':[2.018,31.6],
              'TTAri':[3.301,31.9],
              'EFTuc':[3.5,33.0],

              'SW UMa': [1.364, 30.6],
              'CC Scl': [1.41, 31.4],
              'T Leo': [1.412, 30.6],
              'BZ UMa': [1.632, 31.1],
              'SU UMa': [1.832, 31.9],
              'RX J1831': [4.01, 31.5],
              'WW Cet': [4.22, 31.2],
              'V405 Peg': [4.264, 30.5],
              'EX Dra': [5.039, 29.7]

              }

