import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp
import functools
path='/Users/baotong/Desktop/period/txt_90/'
epoch_file =path + 'SgrA_I_epoch.txt'
T=50000.
w=2*np.pi/T
m=2
fi=0.5
t0=54268974.41
dataname='1502'
time=np.loadtxt(path+ str(dataname)+'.txt')[:,0]

def compute_bin(Tlist, m, w, fi):
    n = np.zeros(m, 'int')
    j = np.floor(m * np.mod(w * Tlist + fi, 2 * np.pi) / (2 * np.pi))
    j.astype(int)
    for u in range(0, m):
        n[u] = np.size(np.extract(j == u, j))
    return n
#print(compute_bin(time,m,w,fi))

def get_T_in_mbins(epoch_file,w,m,fi):
    T=2*np.pi/w
    T_in_perbin = np.zeros(m)
    # 每个bin的总积分时间
    tbin = T/m
    # 每个bin的时间长度
    epoch_info = np.loadtxt(epoch_file)
    t_start = epoch_info[:, 0]
    t_end = epoch_info[:, 1]
    ID = epoch_info[:, 2]

    N_bin_t_start=t_start/tbin+fi/(2*np.pi)
    N_bin_t_end=t_end/tbin+fi/(2*np.pi)
    intN_bin_t_start=np.floor(N_bin_t_start)+1
    intN_bin_t_end=np.floor(N_bin_t_end)
    intN_bin_t_start=intN_bin_t_start.astype(int)
    intN_bin_t_end=intN_bin_t_end.astype(int)
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
    return T_in_perbin

def compute_S(epoch_file,Tlist,w,m,fi):
    n = compute_bin(Tlist, m, w, fi)
    tao=get_T_in_mbins(epoch_file,w,m,fi)
    s_wfi=tao/(sum(tao)/m)
    #print(s_wfi)
    #ln_S_wfi=len(Tlist)-sum(n*s_wfi)
    ln_S_wfi = -sum(n * np.log(s_wfi))
    print(n)
    print(tao)
    print(ln_S_wfi)
    return ln_S_wfi
print(compute_S(epoch_file,time,w,m,fi))
#print(get_T_in_mbins(epoch_file,w,m,fi))
#加上所有整数周期

