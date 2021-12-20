import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp
import functools
import datetime


def get_T_in_mbins(epoch_info,w,m,fi):
    T=2*np.pi/w
    T_in_perbin = np.zeros(m)
    # 每个bin的总积分时间
    tbin = T/m
    # 每个bin的时间长度
    t_start = epoch_info[:, 0];
    t_end = epoch_info[:, 1]

    N_bin_t_start=t_start/tbin+m*fi/(2*np.pi)
    N_bin_t_end=t_end/tbin+m*fi/(2*np.pi)
    intN_bin_t_start=np.floor(N_bin_t_start).astype(int)+1
    intN_bin_t_end=np.floor(N_bin_t_end).astype(int)

    for i in range(len(N_bin_t_start)):
        if intN_bin_t_end[i]>=intN_bin_t_start[i]:
            T_in_perbin+=int((intN_bin_t_end[i]-intN_bin_t_start[i])/m)*tbin
            #print(intN_bin_t_start[i]-1)
            T_in_perbin[np.mod(intN_bin_t_start[i],m)-1]+=(intN_bin_t_start[i]-N_bin_t_start[i])*tbin
            T_in_perbin[np.mod(intN_bin_t_end[i],m)]+=(N_bin_t_end[i]-intN_bin_t_end[i])*tbin
            rest=np.mod(intN_bin_t_end[i]-intN_bin_t_start[i],m)
            for k in range(rest):
                T_in_perbin[int(np.mod((intN_bin_t_start[i] + k), m))] += tbin
            #print(rest)
        else:
            T_in_perbin[np.mod(intN_bin_t_start[i],m)-1]+=(N_bin_t_end[i]-N_bin_t_start[i])*tbin
    print(np.sum(T_in_perbin))
    print(T_in_perbin)
    return T_in_perbin

def get_T_bin_new(epoch_info,w,m,fi):
    T=2*np.pi/w
    T_in_perbin = np.zeros(m)
    dT = T / m
    for i in range(len(epoch_info)):
        t_start=epoch_info[i][0]
        t_stop=epoch_info[i][1]
        fi_start=np.mod((w*t_start+fi),(2*np.pi))/(2*np.pi)
        pm=np.floor(m*fi_start)+1
        temp = pm * dT - fi_start * T
        T_in_perbin[pm] += temp
        while temp < t_stop:
            if t_stop > temp + dT:
                T_in_perbin[pm] += dT
            else:
                T_in_perbin[pm] += t_stop - temp
if __name__=="__main__":
    path='/Users/baotong/Desktop/period_Tuc/txt_all_obs_p90/'
    epoch_info=np.loadtxt(path+'epoch_src_402.txt')
    print(np.sum(epoch_info[:,-1]))
    T_in_perbin=get_T_in_mbins(np.array([[0,1000],[2000,3000]]), 2*np.pi/5000., 8, 0.0)