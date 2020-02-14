import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp
import functools

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
    N_bin_t_start=t_start/tbin+m*fi/(2*np.pi)
    N_bin_t_end=t_end/tbin+m*fi/(2*np.pi)
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
def compute_bin(Tlist, m, w, fi):
    n = np.zeros(m, 'int')
    j = np.floor(m * np.mod(w * Tlist + fi, 2 * np.pi) / (2 * np.pi))
    j.astype(int)
    for u in range(0, m):
        n[u] = np.size(np.extract(j == u, j))
    return n
path='/Users/baotong/Desktop/period_gc/txt/'
dataname='29'
time=np.loadtxt(path+ str(dataname)+'.txt')[:,0]
epoch_file =path + 'ACIS-I_epoch.txt'

m=2
fi=0.2
#period=26785.71429
#print(compute_bin(time,m,1/period,0.5))
N=len(time)

def compute_S(epoch_file,Tlist,w,m,fi):
    n = compute_bin(Tlist, m, w, fi)
    tao=get_T_in_mbins(epoch_file,w,m,fi)
    s_wfi=tao/(sum(tao)/m)
    #if(n[0]-n[1])*(s_wfi[0]-s_wfi[1])<0:
        # print(n)
        # print(s_wfi)
    #ln_S_wfi=len(Tlist)-sum(n*s_wfi)
    ln_S_wfi = -sum(n * np.log(s_wfi))
    return ln_S_wfi

def precompute_binmult(N):
    #return [lg(n!) for n in range(1,N)], 第一项忽略
    # precompute all potential bin factorials
    fbin = np.zeros(int(N) + 1)
    for i in range(2, int(N) + 1):  # n=0 -> define log(n)=0, n=1, log(n)=0, so no need to compute n=0,1
        fbin[i] = fbin[i - 1] + np.log(i)
    return fbin


def Omwf(time,m,w,fi):
    n=compute_bin(time,m,w,fi)
    fbin=precompute_binmult(N+m)
    w_down=0
    for i in range(len(n)):
        w_down+=fbin[n[i]]
    ln_S_wfi = compute_S(epoch_file, time, w, m, fi)
    #print(ln_S_wfi)
    fx=N*np.log(m)+fbin[m-1]+w_down-fbin[N+m-1]+ln_S_wfi-np.log(m)
    return fx

w_all=2 * np.pi*np.arange(1/1000.,1/100.,1e-6)
y=np.zeros(len(w_all))
for i in range(len(w_all)):
    y[i]= Omwf(time, m, w_all[i], fi)
    # for fi in np.arange(0, 2 * np.pi, 0.1):
    #     y[i]+=Omwf(time,m,w_all[i],fi)
    #     print(y[i])
    # y[i]/=len(np.arange(0, 2 * np.pi, 0.1))


x=np.arange(1/1000.,1/100.,1e-6)
plt.step(x,y)
plt.plot(x,x-x+np.log(100),'--')
plt.show()


# def test_lns(x):
#     N=2000
#     m=2
#     fbin = precompute_binmult(N + m)
#     def f(x):
#         s1 = x
#         s2 = 2 - x
#         er=np.random.random(1)*80
#         n1 = N * x / 2.+er
#         print(er)
#         n2 = N-n1
#         f = -n1 * np.log(s1) - n2 * np.log(s2)
#         return f,n1,n2
#     f,n1,n2=f(x)
#     n1=n1.astype(int)
#     n2=n2.astype(int)
#     all=N*np.log(m)+f-fbin[N]+fbin[m]+fbin[n1]+fbin[n2]
#     return all
# x=np.linspace(0.1,1.9,10000)
# #z=test_lns(np.array([0.4999]))
# y=test_lns(x)
# #print(z)
# plt.scatter(x,y)
# plt.show()

