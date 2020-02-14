# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 19:21:40 2019

@author: H.S.Wang
"""

from scipy.optimize import fsolve
import numpy as np

delta = 0.1
w = 1.0/5000
t = 2500


def f(t1):
    return t1-delta/w*np.cos(w*t1) 


result = fsolve(f,2800)