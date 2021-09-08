import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import pandas as pd
import sys
import os
from scipy import interpolate
from scipy.optimize import curve_fit
import astropy.units as u
import astropy.constants as c
from scipy import integrate
import stingray as sr
from stingray.events import EventList
from stingray.lightcurve import Lightcurve
from stingray import Lightcurve, Crossspectrum, sampledata,Powerspectrum,AveragedPowerspectrum
from stingray.simulator import simulator, models
import matplotlib.font_manager as font_manager
from astropy.timeseries import LombScargle
import warnings
from functools import reduce
import csv
import useful_functions as func

f=func.standard_lorentzian
qpo_f=2.7e-4
p_lorentz=[qpo_f, qpo_f/ 16, 100, 2]

p_bend_REJ1034=[4.3e-4, 3.4, 0., 2.3e-3]
bend=func.bending_po(2.67e-4,p_bend_REJ1034)
print(bend)


