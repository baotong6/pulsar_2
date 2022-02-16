#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.timeseries import LombScargle
from astropy.stats import poisson_conf_interval
from scipy import optimize as op
def f(x,a,b):
    logS=a*x+b  #S~V^a
    return logS



def spectrafit(x,y,error):
    popt, pcov = op.curve_fit(f, np.log(x), np.log(y),absolute_sigma=True,sigma=np.log(error)/y)
    perr = np.sqrt(np.diag(pcov))
    logydata=f(np.log(x),popt[0],popt[1])
    ydata=np.exp(logydata)
    return (popt,perr)


v=np.array([3.0,6.0,10.0,15.0])
error=3*np.array([6,5,7,6])
A=np.array([1522.3,710.3,395.1,216.3])
(popt,perr)=spectrafit(v,A,error)

plt.scatter(v,A)

x_test=np.linspace(v[0],v[-1],1000)
(popt,perr)=spectrafit(v,A,error)
print(popt)
print(perr)
plt.plot(x_test,np.exp(f(np.log(x_test),popt[0],popt[1])))
plt.show()