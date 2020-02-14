import numpy as np
import matplotlib.pyplot as plt
from tkinter import _flatten
from astropy.io import fits
import sys
import os
import string

path='/Users/baotong/Desktop/period/simulation/simulation_NSC_10y/result_5.0_0.5/'
therehold=0.9
a=0
for i in range(1,101):
    temp=np.loadtxt(path+'result_sim_{0}.txt'.format(i))
    if temp[2]>therehold:
        a+=1
print(a/100.)


