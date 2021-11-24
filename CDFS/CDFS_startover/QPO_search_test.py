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
import scipy.stats
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
from stingray.modeling import PSDLogLikelihood,PSDPosterior,set_logprior
# define power law component
pl = models.PowerLaw1D()

# fix x_0 of power law component
pl.x_0.fixed = True

# define constant
c = models.Const1D()

# make compound model
plc = pl + c

# parameters for fake data.
alpha = 2.0
amplitude = 5.0
white_noise = 2.0
freq = np.linspace(0.01, 10.0, int(10.0/0.01))
from astropy.modeling.fitting import _fitter_to_model_params
_fitter_to_model_params(plc, [amplitude, alpha, white_noise])
psd_shape = plc(freq)
powers = psd_shape*np.random.chisquare(2, size=psd_shape.shape[0])/2.0

ps = Powerspectrum()
ps.freq = freq
ps.power = powers
ps.df = ps.freq[1] - ps.freq[0]
ps.m = 1
loglike = PSDLogLikelihood(ps.freq, ps.power, plc, m=ps.m)
test_pars = [1, 5, 100]
parest = PSDParEst(ps, fitmethod="L-BFGS-B", max_post=True)

# flat prior for the power law index
p_alpha = lambda alpha: ((-1. <= alpha) & (alpha <= 5.))

# flat prior for the power law amplitude
p_amplitude = lambda amplitude: ((0.01 <= amplitude) & (amplitude <= 10.0))

# normal prior for the white noise parameter
p_whitenoise = lambda white_noise: scipy.stats.norm(2.0, 0.1).pdf(white_noise)

priors = {}
priors["alpha_0"] = p_alpha
priors["amplitude_0"] = p_amplitude
priors["amplitude_1"] = p_whitenoise

starting_pars=[3.0, 1.0, 2.4]
lpost = PSDPosterior(ps.freq, ps.power, plc, priors=priors, m=ps.m)
parest = PSDParEst(ps, fitmethod='BFGS', max_post=True)
res = parest.fit(lpost, starting_pars)

sample = parest.sample(lpost, res.p_opt, cov=res.cov, nwalkers=400,
             niter=100, burnin=200, namestr="psd_modeling_test")

max_power, max_freq, max_ind = parest._compute_highest_outlier(lpost, res)
print(max_power)
pval = parest.calibrate_highest_outlier(lpost, starting_pars, sample=sample,
                                  max_post=True,
                                  nsim=100, niter=100, nwalkers=400,
                                  burnin=200, namestr="test")
print(pval)