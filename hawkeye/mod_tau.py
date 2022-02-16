import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp
import functools
import datetime
import hawkeye as hawk
import math

def get_rate_info(rate_file):
    rate_=np.loadtxt(rate_file)
    rate_=rate_[np.lexsort(rate_[:,::-1].T)]
    #排序
    dx=np.diff(rate_[:,0])
    #diff_x[diff_x>100]=0
    rx=(rate_[1:,0]+rate_[:-1,0])/2
    #梯形中心坐标
    ry=(rate_[1:,1]+rate_[:-1,1])/2
    #梯形中心高度
    rs=dx*ry
    #梯形面积

def get_T_in_mbins(rate_file,w,m,fi):
    rate_info=np.loadtxt(rate_file)
    rate_info=rate_info[np.lexsort(rate_info[:,::-1].T)]

    dx=np.diff(rate_info[:,0])
    #diff_x[diff_x>100]=0
    rx=(rate_info[1:,0]+rate_info[:-1,0])/2
    #梯形中心坐标
    ry=(rate_info[1:,1]+rate_info[:-1,1])/2
    #梯形中心高度
    T_in_perbin = np.zeros(m)
    t=rate_info[:,0]
    rate=rate_info[:,1]
    dx=np.diff(rate_info[:,0])
    T=2*np.pi/w
    T_in_perbin = np.zeros(m)
    # 每个bin的总积分时间
    tbin = T/m
    phase=rx / T + fi / (2 * np.pi)
    phase=np.modf(phase)[0]
    phase*=m
    print(np.sum(dx[dx<100]))
    for i in range(len(phase)-1):
        if dx[i]<100:
            phase_left=int(np.floor(phase[i]))
            T_in_perbin[phase_left]+=dx[i]*ry[i]
        else:
            continue
    print(np.sum(T_in_perbin))
    return T_in_perbin

if __name__ == '__main__':
    dataname=528
    path = '/Users/baotong/Desktop/command/zjc/'
    path_out='/Users/baotong/Desktop/command/zjc/'
    rate_file = path + f'{dataname}_info.txt'
    get_T_in_mbins(rate_file,w=2*np.pi/15000.,m=10,fi=0)
    data_file = path + f'{dataname}.txt'
    epoch_file = path + f'epoch_src_{dataname}.txt'
    src_evt_use=np.loadtxt(data_file)
    epoch_info_use=np.loadtxt(epoch_file)
    # lc=hawk.get_hist(src_evt_use[:,0],len_bin=100)
    # freq = np.arange(1 / 1e7, 0.5 / 100, 1 / (10* 1e6))
    # freq = freq[np.where(freq > 1 / 30000.)]
    # (FP, out_period, max_NormLSP) = hawk.get_LS(lc.time, lc.counts, freq=freq, outpath=None,
    #                                             outname=str(dataname), save=False, show=1)
    print(np.sum(epoch_info_use[:,-1]))
    #
    # hawk.phase_fold(time=src_evt_use[:,0], epoch_info=epoch_info_use, net_percent=0.9, p_test=1000., outpath=None, bin=15,
    #                 label='')