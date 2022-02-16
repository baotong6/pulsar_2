#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import hawkeye as hawk

path='/Users/baotong/eSASS/data/raw_data/47_Tuc/txt/txt_psf75_700163/'
dataname='481';obsid=700163
data_file=path+f'{dataname}_{obsid}.txt'
epoch_file=path+f'epoch_47Tuc_{obsid}.txt'
w_range=2*np.pi*np.arange(1/10000.,1/3000.,1e-5)
[srcid,runtime,Prob,wpeak,wmean,mopt,wconf_lo,wconf_hi,counts]=hawk.GL.write_result(data_file,epoch_file,w_range,dataname='1',if_filter=False)
print(2*np.pi/wpeak)

def get_result_fromid(id_range):
    result_srcid = []
    result_runtime = []
    result_Prob = []
    result_wpeak = []
    result_wmean = []
    result_mopt = []
    result_wconf_lo = []
    result_wconf_hi = []
    result_counts=[]

    for i in id_range:
        res = write_result(i)
        result_srcid.append(res[0])
        result_runtime.append(res[1])
        result_Prob.append(res[2])
        result_wpeak.append(res[3])
        result_wmean.append(res[4])
        result_mopt.append(res[5])
        result_wconf_lo.append(res[6])
        result_wconf_hi.append(res[7])
        result_counts.append(res[8])
    result_wpeak = np.array(result_wpeak)
    result_period = 2 * np.pi / result_wpeak
    result = np.column_stack((result_runtime, result_Prob, result_wpeak,result_wmean, result_period, result_mopt,
                              result_wconf_lo, result_wconf_hi,result_counts))
    print(result)
    #np.savetxt('result_1h-3h_{0}.txt'.format(id_range[0]), result, fmt='%10d %10.2f %10.2f %10.5f %10.5f %10d %10.5f %10.5f')
    np.savetxt(path+'result_1h_{0}.txt'.format(id_range[0]), result,
               fmt='%10.2f %10.5f %10.5f %10.5f %10.5f %10d %10.5f %10.5f %10d')
