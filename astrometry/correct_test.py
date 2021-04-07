import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import string
import math
from jplephem.spk import SPK
from astropy.time import Time
import sofa_test as sofa_test

c=299792.458
#km
G=6.67259*10**(-11)
Msun=1.9891*10**(30)
ra=(5+34./60+31.97232/3600)*15
dec=(22+0./60.+52.0690/3600)*15
distance=6300*9460730472580.8


def jpl(number, jd):
    kernel = SPK.open('/Users/baotong/desktop/pulsar/de430.bsp')
    position, velocity = kernel[0, number].compute_and_differentiate(jd)

    # position=np.array(position)
    # velocity=np.array(velocity)
    velocity = velocity / 86400.0
    return [position, velocity]

# jd = np.array([2457061.5, 2457062.5, 2457063.5, 2457064.5])
# print jpl(3,jd)[0]
def angle(x,y):
    Lx=x.dot(x)
    Ly=y.dot(y)
    cos_angle=x.dot(y)/(Lx*Ly)
    angle=np.arccos(cos_angle)
    return angle

#print angle(np.array([0,1,2]),np.array([3,12,3]))


def UTC_to_TT(t):
    # t must be format "jd"
    time = Time(t, format='jd', scale='utc')
    time_tt = time.tt
    dt = (time.value - time_tt.value) * 86400
    return dt

tobs=time
clock=-UTC_to_TT(time)
dispersion=0
#robs is the vector from the SSB to the observatory
#rearth is the vector from the SSB to the earth geocenter
#rsat is the vector from the earth geocenter to the spacecraft
rsat = sofa_test.r_CIS(orbit,tobs)
#print rsat
rearth=jpl(3,tobs)[0]
vearth=jpl(3,tobs)[1]
rsun=jpl(10,tobs)[0]
n=np.array([np.cos(ra)*np.cos(dec),np.sin(ra)*np.cos(dec),np.sin(dec)])
rpulsar=distance*n
#print rearth
#print rsun
#print rpulsar
th=angle((rearth-rsun),(rpulsar-rsun))
robs = rearth + rsat
geometric = (robs.dot(n)) / c
Einstein_geo = 1.48 * 10 ** (-8)
Einstein = Einstein_geo + (rsat.dot(vearth)) / (c ** 2)
Shapiro = - (2 * G * Msun / ((c * 1000) ** 3)) * np.log(1 + np.cos(th))
# print geometric
# print Einstein
tb = (tobs - 2454466.5) * 86400 + clock - dispersion + geometric + Einstein - Shapiro
# print tb-(tobs-2454466.5)*86400
# print robs
# print geometric
# print clock
# print Einstein
# print Shapiro
