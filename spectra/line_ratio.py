import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import string
import datetime
from scipy.interpolate import lagrange
from scipy import interpolate
from scipy.fftpack import fft,ifft
import scipy.signal as ss
import scipy.stats as stats
import pandas as pd
import read_csv as data

def get_line(ID,path):
    bins = np.linspace(np.log10(0.01), np.log10(5), 100)  # note that +1 to nbins
    in64=[];in67=[];in70=[]
    in70on67=[];low_70on67=[];high_70on67=[]
    in64on67=[];low_64on67=[];high_64on67=[]
    i=0
    while i <(len(ID)):
        id=ID[i]
        fake_name='fake_{0}.fits'.format(str(id))
        fake_file=fits.open(path+fake_name)
        line_in_64=fake_file[1].data['norm6']
        line_in_67=fake_file[1].data['norm9']
        line_in_70=fake_file[1].data['norm12']
        plt.figure()
        plt.title('{0}'.format(str(id)))
        # print(line_in_64)
        # print(line_in_67)
        a64=plt.hist(line_in_64,bins=100,histtype = 'step')
        a67=plt.hist(line_in_67,bins=100,histtype = 'step')
        a70=plt.hist(line_in_70,bins=100,histtype = 'step')
        in64.append((a64[1][np.where(a64[0]==np.max(a64[0]))[0]]/2.+a64[1][np.where(a64[0]==np.max(a64[0]))[0]+1]/2.)[0])
        in67.append((a67[1][np.where(a67[0]==np.max(a67[0]))[0]]/2.+a67[1][np.where(a67[0]==np.max(a67[0]))[0]+1]/2.)[0])
        in70.append((a70[1][np.where(a70[0]==np.max(a70[0]))[0]]/2.+a70[1][np.where(a70[0]==np.max(a70[0]))[0]+1]/2.)[0])
        plt.show()
        line_64on67 = (6.4 / 6.7) * (line_in_64 / line_in_67)
        line_70on67 = (7.0 / 6.7) * (line_in_70 / line_in_67)
        # line_64on67=np.delete(line_64on67,np.where(line_64on67>2))
        line_64on67_sort = np.sort(line_64on67)
        line_70on67_sort = np.sort(line_70on67)

        low_70on67.append(line_70on67_sort[int(0.32 * len(line_70on67_sort))])
        high_70on67.append(line_70on67_sort[int(0.68 * len(line_70on67_sort))])
        in70on67_mod = plt.hist(line_70on67_sort, bins = bins, histtype = 'step')

        low_64on67.append(line_64on67_sort[int(0.32 * len(line_64on67_sort))])
        high_64on67.append(line_64on67_sort[int(0.68 * len(line_64on67_sort))])
        in64on67_mod = plt.hist(line_64on67_sort, bins = bins, histtype = 'step')

        in64on67.append(((in64on67_mod[1][np.where(in64on67_mod[0] == np.max(in64on67_mod[0]))[0]] +
                          in64on67_mod[1][np.where(in64on67_mod[0] == np.max(in64on67_mod[0]))[0] + 1]) / 2.)[0])

        in70on67.append(((in70on67_mod[1][np.where(in70on67_mod[0] == np.max(in70on67_mod[0]))[0]] +
                          in70on67_mod[1][np.where(in70on67_mod[0] == np.max(in70on67_mod[0]))[0] + 1]) / 2.)[0])

        plt.show()
        i += 1
    # in70on67=np.array(in70on67)
    # low_70on67=np.array(low_70on67)
    # high_70on67=np.array(high_70on67)
    return [in64on67,in70on67,low_64on67,low_70on67,high_64on67,high_70on67]
ID=np.array([194,500,114,153,191,16,118])
path='/Users/baotong/Desktop/period_LW/spectra_stack_src/'
result=get_line(ID,path)
print(result[0])
print(result[1])
print(result[2])
print(result[3])
print(result[4])
print(result[5])


