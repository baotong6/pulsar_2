# -*- coding: utf-8 -*-
"""
Created on Sun Sep 20 18:13:40 2015
@author: Felix Darvas
@author: baotong
compute the Gregory-Laredo algorithm on arival times
This function computes the likelihood of a set of arrival times originating
from a periodic system rather than constant rate (poisson) background noise
based on
Gregory, P. C. and Thomas. J. Loredo, 1992,
"A New Method For The Detection Of A Periodic Signal Of Unknown Shape And Period"
in
The Astrophysical Journal, Astrophysical J., 398, p.146
inputs:
Tlist    -  list of arrival times, numpy int array
m_max    -  max number of bins typically 12-15, we use 12 as default
w_range  -  frequency range to scan numpy float array of frequency values
           default is  w_lo=20*pi/T at delta w = pi/T to w_hi=pi*N/T
           where N=#arrival times, T=observation time
ni       - number of integration steps, default ni=10
parallel - use parallel execution - default is off
outut:
O_period - Odds ratio for a periodic process vs. constant rate process
p_period - probability of a periodic process 0<=p_period<=1
m_opt    - optimal bin size 1<= m_opt <=m_max
S        - The probability spectrum
w        - The frequency range for S
w_peak   - the expected frequency of the process
w_conf   - 95% confidence interval of w_peak
"""
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp
import functools
import datetime
import GL_algorithm_esass as GL
import funcs_sim_timeseries as funcs_sim
import sys

def write_result(dataname,cts_num,amp_num):
    path='/Users/baotong/eSASS/data/raw_data/47_Tuc/txt/'
    epoch_file = path + 'epoch_47Tuc.txt'
    epoch_info=np.loadtxt(epoch_file)
    cts_rate_standard=0.01
    amp_standard=1.0
    time = funcs_sim.get_epoch_time_series(cts_rate = cts_rate_standard*cts_num, period =3600.,
                                           amp =amp_standard*amp_num,epoch=epoch_info,model = 'sin',varydelta=5)
    # time = funcs_sim.get_epoch_time_series(cts_rate = cts_rate_standard*cts_num, period =3600.,
    #                                        amp =amp_standard*amp_num,epoch=epoch_info,model = 'sin')
    print(len(time))

    w_range = 2 * np.pi * np.arange(1. /10000, 1. /3000., 5e-7)
    starttime = datetime.datetime.now()
    GL_R = GL.compute_GL(time, epoch_file,w_range=w_range, m_max=12, parallel=True)
    endtime = datetime.datetime.now()
    srcid = dataname
    runtime = (endtime - starttime).seconds
    Prob = GL_R[1]
    wpeak = GL_R[5]
    mopt = GL_R[2]
    wconf_lo = GL_R[7][0]
    wconf_hi = GL_R[7][1]
    O_per_w=GL_R[9]
    p_per_w=GL_R[10]
    N_cts = len(time)

    return [srcid, runtime, Prob, wpeak, mopt, wconf_lo, wconf_hi,O_per_w,p_per_w,N_cts]

def get_result_fromid(dataname,cts_num,amp_num):
    path_out='/Users/baotong/eSASS/data/raw_data/47_Tuc/simulation/vary_delta_GL/'
    result_srcid=int(dataname)
    res = write_result(dataname = dataname,cts_num=cts_num,amp_num=amp_num)
    #result_srcid=res[0]
    result_runtime=res[1]
    result_Prob=res[2]
    result_wpeak=res[3]
    result_mopt=res[4]
    result_wconf_lo=res[5]
    result_wconf_hi=res[6]
    result_O_per_w=res[7]
    result_p_per_w=res[8]
    result_period = 2 * np.pi / result_wpeak
    result_N_cts=res[9]
    result = np.column_stack((result_srcid, result_runtime, result_Prob, result_wpeak, result_period, result_mopt,
                              result_wconf_lo, result_wconf_hi,result_N_cts))
    #print(result)
    #np.savetxt('result_1h-3h_{0}.txt'.format(id_range[0]), result, fmt='%10d %10.2f %10.5f %10.10f %10.5f %10d %10.10f %10.10f %10.5f %10.5f')
    np.savetxt(path_out+'result_{0}_{1}/result_sim_{2}.txt'.format(str(cts_num),str(amp_num),str(dataname)), result,
               fmt='%10d %10.2f %10.5f %10.10f %10.5f %10d %10.10f %10.10f %10d')

if __name__ == '__main__':
    cts_range = [1.5]
    amp_range = [0.2, 0.3, 0.4, 0.5]
    for x in range(len(cts_range)):
        cts_num=cts_range[x]
        for y in range(len(amp_range)):
            amp_num=amp_range[y]
            for i in range(100):
                get_result_fromid(i+1,cts_num=cts_num,amp_num=amp_num)
#cand_id=np.linspace(1,518,518)
# cand_id=np.loadtxt('cand_id.txt')
# cand_id=cand_id.astype(int)
# def choose_id(a1, a2):
#     a1=int(a1)
#     a2=int(a2)
#     get_result_fromid(cand_id[a1:a2])
#
# import sys
# if __name__ == '__main__':
#     choose_id(sys.argv[1],sys.argv[2])

#get_result_fromid([29])
