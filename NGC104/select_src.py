import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
import pandas as pd
import sys
import os
from tkinter import _flatten
import esass.funcs_fits2txt as funcs
from esass.funcs_fits2txt import Circle
from scipy import interpolate
from astropy.coordinates import SkyCoord
from astropy import units as u

ra_center=6.022318
dec_center=-72.081443
inter_radius=20  ##角分
path='/Volumes/pulsar/47Tuc/merge_data/timing/txt_startover/txt_all_obs_p50/'