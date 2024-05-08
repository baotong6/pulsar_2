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

starttime = datetime.datetime.now()
def compute_bin(Tlist, m, w, fi):
    n = np.zeros(m, 'int')
    j = np.floor(m * np.mod(w * Tlist + fi, 2 * np.pi) / (2 * np.pi))
    j.astype(int)
    for u in range(0, m):
        n[u] = np.size(np.extract(j == u, j))
    return n

def compute_factor(N, m, v):
    # compute m^N *(N+m-1)! / (2pi*v*(m-1)!
    # which is used to scale the multiplicity function (see eqn.5.25 of the paper )
    # we return the log of the integral scale here to avoid numerical overflow
    f1 = N * np.log(m)  # log of m^N
    f2 = np.sum(np.log(np.arange(1, N + m)))  # log of (N+m-1)!
    #f3 = np.sum(np.log(np.arange(1, m + 1)))  # log of m!
    f3 = np.sum(np.log(np.arange(1, m )))
    #should be (m-1)! --Tong
    f = f1 + f3 - f2 - np.log(2 * np.pi * v)
    return f


def precompute_binmult(N):
    #return [lg(n!) for n in range(1,N)], 第一项忽略
    # precompute all potential bin factorials
    fbin = np.zeros(int(N) + 1)
    for i in range(2, int(N) + 1):  # n=0 -> define log(n)=0, n=1, log(n)=0, so no need to compute n=0,1
        fbin[i] = fbin[i - 1] + np.log(i)
    return fbin

def compute_S(epoch_file,Tlist,w,m,fi):
    n = compute_bin(Tlist, m, w, fi)
    tao=get_T_in_mbins(epoch_file,w,m,fi)
    s_wfi=tao/(sum(tao)/m)
    ln_S_wfi=-sum(n*np.log(s_wfi))
    #print(ln_S_wfi)
    return ln_S_wfi

def compute_W_scaled(Tlist, epoch_file,m, w, fi, factor, fbin):
    # compute the scaled multiplicity 1/Wm(w,fi) (see eqn. 5.25)
    # actually, it compute only n1!n2!...nm!,cause N! is been reduced in compute_factor --Tong
    # note that for large arrival time numbers the multiplicity becomes
    # excessively large. Since W_m(fi,w) is never needed explicitly,
    # we use the scaled version.
    # input factor is the log of the actual factor
    n = compute_bin(Tlist, m, w, fi)  # find bin histogram for a given bin number m, frequency w and fi
    f = 0
    for i in range(0, m):
        # f=f+np.sum(np.log(np.arange(2,n[i]+1)))
        f = f + fbin[n[i]]
    ln_S_wfi=compute_S(epoch_file,Tlist,w,m,fi)
    y = np.exp(f + factor+ln_S_wfi)
    return y


def compute_Om1(w, epoch_file,Tlist, m, factor, fbin, ni):
    # compute  specific odds-ratios (eqn. 5.25)
    p = np.arange(0, ni) / float(ni) * 2 * np.pi / m
    #p 即为输入的fi
    # intgration range, only integrate over a single bin, as values of the integral repeat over bins
    y = np.zeros(np.size(p), 'float')  # array to hold values of W_scaled over the integration range
    for i in range(0, np.size(y)):
        y[i] = compute_W_scaled(Tlist,epoch_file, m, w, p[i], factor, fbin)
    #print(y)
    #print(np.trapz(y, p) * m)
    return np.trapz(y, p) * m  # return intregrated W_Scaled


def compute_Om1wPar(Tlist, epoch_file,m_max, w, fa, fbin, ni):  # compute odds-ratios for bins and frequencies
    # parallel version
    Om1w = np.zeros((m_max, np.size(w)), 'float')  # odds ratio matrix
    pool = mp.Pool()  # use all workers

    for m in range(0, m_max):
        Om1w[m, :] = pool.map(functools.partial(compute_Om1, epoch_file=epoch_file,Tlist=Tlist, m=(m + 1), factor=fa[m], fbin=fbin, ni=ni), w)
    #print(Om1w)
    return Om1w


def compute_Om1w(Tlist, m_max, w, fa, fbin, ni):  # compute odds-ratios for bins and frequencies
    Om1w = np.zeros((m_max, np.size(w)), 'float')  # odds ratio matrix
    for m in range(0, m_max):
        for wi in range(0, np.size(w)):
            Om1w[m, wi] = compute_Om1(w[wi], epoch_file,Tlist, m + 1, fa[m], fbin, ni)
    return Om1w


def compute_GL(Tlist,epoch_file, m_max=12, w_range=None, ni=10, parallel=False):
    # initialize output values
    O_period = None
    p_period = None
    m_opt = None
    S = None
    w = None
    w_peak = None
    w_mean = None
    w_conf = None
    N = float(np.size(Tlist))  # need float value to avoid int/int
    if N > 0:
        # compute GL algorithm
        fbin = precompute_binmult(N)
        v = m_max - 1
        T = float(np.max(Tlist))  # duration of the observation
        if w_range is None:  # use default frequencies
            w_hi = np.pi * N / T  # max default frequency
            w_lo = np.minimum(20, N / 10) * np.pi / T  # min default frequency
            dw = np.pi / T / 10  # step size
            w = np.arange(w_lo, w_hi, dw)
            if np.size(w) < 2:
                print
                "error "
                raise ValueError('bad arrival time list')
        else:  # use user supplied frequency vector
            w = w_range
            w_hi = np.max(w_range)
            w_lo = np.min(w_range)
        if w_lo == 0:
            # cannot have w_lo =0
            print("minimum frequency cannot be 0!\n")
            return

        fa = np.zeros(m_max)
        for m in range(0, m_max):  # precompute factors for each m
            fa[m] = compute_factor(N, m + 1, v)
        if parallel:
            Om1w = compute_Om1wPar(Tlist,epoch_file, m_max, w, fa, fbin, ni)
        else:
            Om1w = compute_Om1w(Tlist,epoch_file, m_max, w, fa, fbin, ni)
        pw = 1 / w / np.log(w_hi / w_lo)
        O1m = np.zeros(m_max)
        for i in range(0, m_max):  # intgreate odd ratios across frequencies
            O1m[i] = np.trapz(pw * Om1w[i], w)
        O_period = np.sum(O1m[1:])  # integrated odds ratio
        # print(Om1w[1:])
        # print(O1m[1:])
        # print(O_period)
        p_period = O_period / (1 + O_period)  # likelihood of periodic event

        m_opt = np.argmax(O1m)  # find optimum bin number, i.e largest odds-ratio
        S = Om1w[m_opt] / w  # compute Spectral probability
        m_opt = m_opt +1  # start bin index with 1
        C = np.trapz(S, w)  # compute normalization
        S = S / C  # normalized probability
        #
        S_up = [S for i in range(m_max)]
        S_up=np.array(S_up)
        for i in range(1,m_max+1):
            S_up[i-1]=Om1w[i-1]/w
            C_temp= np.trapz(S_up[i-1], w)
            S_up[i-1]=S_up[i-1]/C_temp
        pMm=O1m/sum(O1m)
        S_final=[0 for i in range(len(S_up[0]))]
        for i in range(m_max):
            S_final+=pMm[i]*S_up[i]

        cdf = np.array(S)
        for i in range(0, np.size(S)):
            cdf[i] = np.trapz(S[0:i], w[0:i])
        #print(cdf)
        wr = np.extract(np.logical_and(cdf > .025, cdf < .975), w)
        w_peak = w[np.argmax(S)]
        w_mean = np.trapz(S * w, w)
        #print(S)

        cdf_f = np.array(S_final)
        for i in range(0, np.size(S_final)):
            cdf_f[i] = np.trapz(S_final[0:i], w[0:i])
        wr = np.extract(np.logical_and(cdf_f > .025, cdf_f < .975), w)

        w_peak = w[np.argmax(S_final)]
        w_mean = np.trapz(S_final * w, w)
        #print(S_final)
        O_per_w=0

        for m in range(0, m_max):
            O_per_w+=Om1w[m][np.argmax(S_final)]

        p_per_w =O_per_w/(1.+O_per_w)

        if np.size(wr) > 0:
            w_conf = [np.min(wr), np.max(wr)]
        else:
            w_conf = [w_peak, w_peak]

        # plt.subplot(221)
        # plt.step(w, S_final)
        # plt.subplot(222)
        # plt.step(w,cdf_f)
        # plt.subplot(223)
        # plt.step(w, S)
        # plt.subplot(224)
        # plt.step(w, cdf)
        # plt.show()
        return O_period, p_period, m_opt, S, w, w_peak, w_mean, w_conf, cdf,O_per_w,p_per_w

    else:# throw an error
        print("No valid arrival time array provided!\n")
        return O_period, p_period, m_opt, S, w, w_peak, w_mean, w_conf,cdf,O_per_w,p_per_w

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
path_out='/Users/baotong/nustar/simulation/'
path='/Users/baotong/nustar/simulation/evt_1/'
# print(sum(get_T_in_mbins(epoch_file,2*np.pi/55000.,10,0.6)))
def write_result(dataname):
    # time=np.loadtxt(path+ str(dataname)+'.txt')[:,0]
    # epoch_file = path + 'epoch_src_' + str(dataname) + '.txt'
    time=np.loadtxt(path+ 'src1_'+str(dataname)+'.txt')[:,0]
    epoch_file = path + 'epoch_src_1' + '.txt'
    w_range = 2 * np.pi * np.arange(1. / 500, 1. /300., 1.e-6)
    starttime = datetime.datetime.now()
    GL_R = compute_GL(time, epoch_file,w_range=w_range, m_max=12, parallel=True)
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

def get_result_fromid(dataname):
    dataname = int(dataname)
    # cts_num = float(cts_num)
    # amp_num = float(amp_num)
    res = write_result(dataname = dataname)
    result_srcid=res[0]
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
    np.savetxt(path_out+'res_src1/result_src1_{0}.txt'.format(str(dataname)), result,
               fmt='%10d %10.2f %10.5f %10.10f %10.5f %10d %10.10f %10.10f %10d')

if __name__ == '__main__':
    for i in range(11,20):
        get_result_fromid(i)

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
