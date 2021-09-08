import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import pandas as pd
import sys
import os
from tkinter import _flatten

def where_region(x,y,reg):
    r=np.array((x-reg[0],y-reg[1]))
    len_r=np.sqrt(r[0]**2+r[1]**2)
    temp=len_r-reg[2]
    return np.where(temp<=0)

def delete_photon_ID(time,energy,ID,emin=500,emax=8000):
    i=0
    while i < len(energy):
        if energy[i]>emax or energy[i]<emin:
        #if energy[i] > 8000 or energy[i] < 2000:
            energy=np.delete(energy,i)
            time=np.delete(time,i)
            ID=np.delete(ID,i)
            i=i-1
        i=i+1
    return [time,energy,ID]