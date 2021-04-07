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
    # sigma=0.683
    while i <(len(ID)):
        id=ID[i]
        fake_name='fake_{0}.fits'.format(str(id))
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
        bins = np.linspace(0.1, 5, 51)

        # print(line_in_64)
        # print(line_in_67)
        a64=plt.hist(line_in_64,bins=100,histtype = 'step')
        a67=plt.hist(line_in_67,bins=100,histtype = 'step')
        a70=plt.hist(line_in_70,bins=100,histtype = 'step')
        plt.close()
        #plt.legend([6.4,6.7,7.0])
        # in64.append((a64[1][np.where(a64[0]==np.max(a64[0]))[0]]/2.+a64[1][np.where(a64[0]==np.max(a64[0]))[0]+1]/2.)[0])
        # in67.append((a67[1][np.where(a67[0]==np.max(a67[0]))[0]]/2.+a67[1][np.where(a67[0]==np.max(a67[0]))[0]+1]/2.)[0])
        # in70.append((a70[1][np.where(a70[0]==np.max(a70[0]))[0]]/2.+a70[1][np.where(a70[0]==np.max(a70[0]))[0]+1]/2.)[0])
        #plt.show()
        plt.figure()
        plt.title('{0}'.format(str(id)))
        line_64on67 = (6.4 / 6.7) * (line_in_64 / line_in_67)
        line_70on67 = (7.0 / 6.7) * (line_in_70 / line_in_67)

        line_64on67_sort = np.sort(line_64on67)
        line_70on67_sort = np.sort(line_70on67)

        low_70on67.append(line_70on67_sort[int((0.05) * len(line_70on67_sort))])
        in70on67.append(line_70on67_sort[int((0.5) * len(line_70on67_sort))])
        high_70on67.append(line_70on67_sort[int(0.84 * len(line_70on67_sort))])
        in70on67_mod = plt.hist(line_70on67_sort, bins = bins, histtype = 'step')

        low_64on67.append(line_64on67_sort[int(0.16 * len(line_64on67_sort))])
        in64on67.append(line_64on67_sort[int(0.5 * len(line_64on67_sort))])
        high_64on67.append(line_64on67_sort[int(0.84 * len(line_64on67_sort))])
        in64on67_mod = plt.hist(line_64on67_sort, bins = bins, histtype = 'step')

        plt.legend(['7.0/6.7','6.4/6.7'])
        # in64on67.append(((in64on67_mod[1][np.where(in64on67_mod[0] == np.max(in64on67_mod[0]))[0]] +
        #                   in64on67_mod[1][np.where(in64on67_mod[0] == np.max(in64on67_mod[0]))[0] + 1]) / 2.)[0])
        #
        # in70on67.append(((in70on67_mod[1][np.where(in70on67_mod[0] == np.max(in70on67_mod[0]))[0]] +
        #                   in70on67_mod[1][np.where(in70on67_mod[0] == np.max(in70on67_mod[0]))[0] + 1]) / 2.)[0])
        i += 1
        #plt.show()


    in64on67=np.array(in64on67)
    low_64on67=np.array(low_64on67)
    high_64on67=np.array(high_64on67)
    in70on67=np.array(in70on67)
    low_70on67=np.array(low_70on67)
    high_70on67=np.array(high_70on67)
    in67=np.array(in67)
    plt.close()

    return [in64on67,in70on67,low_64on67,low_70on67,high_64on67,high_70on67,in67]
ID=np.array([214,2560,116,1502,2380
,2532,1266,1206,1624,3357,1219
,2422,1853,3120,1133,2730,1084
,2525,1529,2157,1677,3564,2338,3107,1634,3401,3596])
ID=np.array([3596])

path='/Users/baotong/Desktop/period/spectra_stack_src/'
result=get_line(ID,path)

for i in range(len(ID)):
#for i in range(index,index+1):
    # if result[2][i]>0 or result[3][i]>0 or result[6][i]>0 and result[5][i]>0:
    # if result[2][i] > 0 or result[3][i] > 0 :
        print('ID={0}'.format(ID[i]))
        print('in64on67={0}'.format(result[0][i]))
        print('low_err_64on67={0}'.format(result[2][i] - result[0][i]))
        print('high_err_64on67={0}'.format(result[4][i] - result[0][i]))
        print('in70on67={0}'.format(result[1][i]))
        print('low_err_70on67={0}'.format(result[3][i] - result[1][i]))
        print('high_err_70on67={0}'.format(result[5][i] - result[1][i]))
        print('in67={0}'.format(result[6][i]))
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

