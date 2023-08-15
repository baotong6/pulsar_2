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
from astropy.modeling import models, fitting


def func_gaosi(x, miu, sigma):
    return 1 / np.sqrt(2 * np.pi) / sigma * np.exp(-(x - miu) ** 2 / 2 / sigma ** 2)


def fit_gaosi(data, bins):
    hx, xedge = np.histogram(data, bins)
    xedge = (xedge[1:] + xedge[:-1]) / 2

    g_init = models.Gaussian1D(amplitude=np.max(hx), mean=np.mean(data), stddev=np.std(data))
    fit_g = fitting.LevMarLSQFitter()
    g = fit_g(g_init, xedge, hx)

    return g.mean.value, g.stddev.value, g

def get_line(ID,path):
    bins = np.linspace(np.log10(0.01), np.log10(5), 100)  # note that +1 to nbins
    in64=[];in67=[];in70=[]
    in70on67=[];low_70on67=[];high_70on67=[]
    in64on67=[];low_64on67=[];high_64on67=[]
    i=0
    # sigma=0.683
    while i <(len(ID)):
        id=ID[i]
        fake_name='fake_{0}_3line.fit'.format(str(id))
        #fake_name='fake_all_but_two.fits'.format(str(id))
        fake_file=fits.open(path+fake_name)
        line_in_64=fake_file[1].data['norm6']
        line_in_67=fake_file[1].data['norm9']
        line_in_70=fake_file[1].data['norm12']

        line_in_64_sort=np.sort(line_in_64)
        line_in_67_sort = np.sort(line_in_67)
        line_in_70_sort = np.sort(line_in_70)
        in67.append(line_in_67_sort[int(0.50 * len(line_in_67_sort))])

        low_64=line_in_64_sort[int(0.05 * len(line_in_64_sort))]
        low_67=line_in_67_sort[int(0.05 * len(line_in_67_sort))]
        low_70 = line_in_70_sort[int(0.05 * len(line_in_70_sort))]
        # if low_64>0 or low_67>0 or low_70>0:
        #     print(ID[i])
        # bins = np.linspace(0, 5, 51)
        # a64=plt.hist(line_in_64,bins=100,histtype = 'step')
        # a67=plt.hist(line_in_67,bins=100,histtype = 'step')
        # a70=plt.hist(line_in_70,bins=100,histtype = 'step')
        # plt.close()

        plt.figure()
        plt.title('{0}'.format(str(id)))
        line_64on67 = (6.4 / 6.7) * (line_in_64 / line_in_67)
        line_70on67 = (7.0 / 6.7) * (line_in_70 / line_in_67)

        line_64on67_sort = np.sort(line_64on67)
        line_70on67_sort = np.sort(line_70on67)
        low_70on67.append(line_70on67_sort[int((0.16) * len(line_70on67_sort))])
        in70on67.append(line_70on67_sort[int((0.5) * len(line_70on67_sort))])
        high_70on67.append(line_70on67_sort[int(0.84 * len(line_70on67_sort))])
        in70on67_mod = plt.hist(line_70on67_sort, bins = bins, histtype = 'step')

        low_64on67.append(line_64on67_sort[int(0.16 * len(line_64on67_sort))])
        in64on67.append(line_64on67_sort[int(0.5 * len(line_64on67_sort))])
        high_64on67.append(line_64on67_sort[int(0.84 * len(line_64on67_sort))])
        in64on67_mod = plt.hist(line_64on67_sort, bins = bins, histtype = 'step')
        plt.plot([line_70on67_sort[int((0.16) * len(line_70on67_sort))],line_70on67_sort[int((0.16) * len(line_70on67_sort))]],[0,200],'--')
        plt.plot([line_70on67_sort[int((0.5) * len(line_70on67_sort))],line_70on67_sort[int((0.5) * len(line_70on67_sort))]],[0,200],'--')
        plt.plot([line_70on67_sort[int((0.84) * len(line_70on67_sort))],line_70on67_sort[int((0.84) * len(line_70on67_sort))]],[0,200],'--')
        plt.legend(['7.0/6.7','6.4/6.7'])
        # in64on67.append(((in64on67_mod[1][np.where(in64on67_mod[0] == np.max(in64on67_mod[0]))[0]] +
        #                   in64on67_mod[1][np.where(in64on67_mod[0] == np.max(in64on67_mod[0]))[0] + 1]) / 2.)[0])
        #
        # in70on67.append(((in70on67_mod[1][np.where(in70on67_mod[0] == np.max(in70on67_mod[0]))[0]] +
        #                   in70on67_mod[1][np.where(in70on67_mod[0] == np.max(in70on67_mod[0]))[0] + 1]) / 2.)[0])
        i += 1
        plt.show()


    in64on67=np.array(in64on67)
    low_64on67=np.array(low_64on67)
    high_64on67=np.array(high_64on67)
    in70on67=np.array(in70on67)
    low_70on67=np.array(low_70on67)
    high_70on67=np.array(high_70on67)
    in67=np.array(in67)
    plt.close()

    return [in64on67,in70on67,low_64on67,low_70on67,high_64on67,high_70on67,in67]


def read_line_70on67(ID,path):
    fake_name = 'fake_{0}_3line.fit'.format(str(ID))
    fake_file = fits.open(path + fake_name)
    line_in_67 = fake_file[1].data['norm9']
    line_in_70 = fake_file[1].data['norm12']
    validindex=np.where((line_in_67>0)&(line_in_70>0))[0]
    print(len(validindex))
    line_70on67 = (7.0 / 6.7) * (line_in_70[validindex] / line_in_67[validindex])

    line_70on67_sort = np.sort(line_70on67)
    bins=np.logspace(np.log10(1e-5),np.log10(2),100)
    plt.hist(line_70on67_sort,bins=bins,histtype='step')
    plt.plot(
        [line_70on67_sort[int((0.16) * len(line_70on67_sort))], line_70on67_sort[int((0.16) * len(line_70on67_sort))]],
        [0, 200], '--')
    plt.plot(
        [line_70on67_sort[int((0.5) * len(line_70on67_sort))], line_70on67_sort[int((0.5) * len(line_70on67_sort))]],
        [0, 200], '--')
    plt.plot(
        [line_70on67_sort[int((0.84) * len(line_70on67_sort))], line_70on67_sort[int((0.84) * len(line_70on67_sort))]],
        [0, 200], '--')
    plt.show()
    line_70on67_mod=line_70on67_sort[np.where((line_70on67_sort>bins[0])&(line_70on67_sort<bins[-1]))[0]]
    print(len(line_70on67_mod))
    low_70on67=line_70on67_mod[int((0.16) * len(line_70on67_mod))]
    in70on67=line_70on67_mod[int((0.5) * len(line_70on67_mod))]
    high_70on67=line_70on67_mod[int(0.84 * len(line_70on67_mod))]

    return [low_70on67,in70on67,high_70on67]

if __name__=='__main__':
    ID='lsrc_in'
    path='/Users/baotong/Desktop/period_terzan5/spectra_startover/'
    result=read_line_70on67(ID,path)
    print('ID={0}'.format(ID))
    print('in70on67={0}'.format(result[1]))
    print('low_err_70on67={0}'.format(result[1] - result[0]))
    print('high_err_70on67={0}'.format(result[2] - result[1]))
# print(result[0])
# print(result[2]-result[0])
# print(result[4]-result[0])
# print(result[1])
# print(result[3]-result[1])
# print(result[5]-result[1])

# for i in range(len(ID)):
#     if (result[2][i]-result[0][i])<0 and (result[4][i]-result[0][i])>0:
#         print(ID[i])
#
#     if (result[3][i]-result[0][i])<0 and (result[5][i]-result[0][i])>0:
#         print(ID[i])

