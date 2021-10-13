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
from scipy import interpolate
import stingray as sr
from stingray.events import EventList
from stingray.lightcurve import Lightcurve
from stingray import Lightcurve, Crossspectrum, sampledata,Powerspectrum,AveragedPowerspectrum
from stingray.simulator import simulator, models
from stingray.modeling import PSDParEst
import matplotlib.font_manager as font_manager
from astropy.timeseries import LombScargle
import warnings
from functools import reduce
from CDFS.CDFS_startover import useful_functions as func
from CDFS.CDFS_startover import sim_psd as sim
from scipy import integrate
from astropy.modeling import models
from astropy.modeling.fitting import _fitter_to_model_params
from stingray.modeling import PSDLogLikelihood

def plot():
    return 1