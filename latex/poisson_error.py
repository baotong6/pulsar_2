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
from astropy.stats import poisson_conf_interval

def xony_poserr(x,y):
    ## x and y must be integers
    u=x/y
    if x<0:x=-x
    x_1sigma = poisson_conf_interval(x, interval='frequentist-confidence').T
    y_1sigma = poisson_conf_interval(y, interval='frequentist-confidence').T
    x_low=x_1sigma[0];x_high=x_1sigma[1];y_low=y_1sigma[0];y_high=y_1sigma[1]
    x_err = min((x_high-x),(x-x_low))
    y_err = min((y_high - y), (y - y_low))
    u_err=np.sqrt((x_err**2*y**2+y_err**2*x**2)/y**4)
    return (u,u_err)
# (u,u_err)=xony_poserr(6,7)
# print(u,u_err)
if __name__=="__main__":
    (u,u_err)=xony_poserr(6,7)
    print(u,u_err)