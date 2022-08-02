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
import hawkeye as hawk

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

def filter_energy(time,energy,band):
    T=time
    E=energy
    i=0
    if len(time)!=len(energy):
        print('error')
        return None
    else:
        while i <len(E):
            if E[i]<band[0] or E[i]>band[1]:
                E=np.delete(E,i)
                T=np.delete(T,i)
            else:
                i+=1
    return T

def precompute_binmult(N):
    #return [lg(n!) for n in range(1,N)], 第一项忽略
    # precompute all potential bin factorials
    fbin = np.zeros(int(N) + 1)
    for i in range(2, int(N) + 1):  # n=0 -> define log(n)=0, n=1, log(n)=0, so no need to compute n=0,1
        fbin[i] = fbin[i - 1] + np.log(i)
    return fbin

def compute_S(epoch_info,Tlist,w,m,fi):
    n = compute_bin(Tlist, m, w, fi)
    tao=get_T_in_mbins(epoch_info,w,m,fi)
    s_wfi=tao/(sum(tao)/m)
    ln_S_wfi=-sum(n*np.log(s_wfi))
    return ln_S_wfi

def compute_W_scaled(Tlist, m, w, fi, factor, fbin,epoch_info):
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
    if epoch_info.ndim==2:
        ln_S_wfi=compute_S(epoch_info,Tlist,w,m,fi)
    elif epoch_info==None:ln_S_wfi=0

    y = np.exp(f + factor+ln_S_wfi)
    #if y >1:
        #print(y)
        #print(f)
        #print(factor)

    return y

def compute_Om1(w, Tlist, m, factor, fbin, ni,epoch_info):
    # compute  specific odds-ratios (eqn. 5.25)
    p = np.arange(0, ni) / float(ni) * 2 * np.pi / m
    #p 即为输入的fi
    # intgration range, only integrate over a single bin, as values of the integral repeat over bins
    y = np.zeros(np.size(p), 'float')  # array to hold values of W_scaled over the integration range
    for i in range(0, np.size(y)):
        y[i] = compute_W_scaled(Tlist, m, w, p[i], factor, fbin,epoch_info)
    return np.trapz(y, p) * m  # return intregrated W_Scaled

def compute_Om1wPar(Tlist, m_max, w, fa, fbin, ni,epoch_info):  # compute odds-ratios for bins and frequencies
    # parallel version
    Om1w = np.zeros((m_max, np.size(w)), 'float')  # odds ratio matrix
    pool = mp.Pool()  # use all workers

    for m in range(0, m_max):
        Om1w[m, :] = pool.map(functools.partial(compute_Om1, Tlist=Tlist, m=(m + 1), factor=fa[m], fbin=fbin, ni=ni,epoch_info=epoch_info), w)
    return Om1w

def compute_Om1w(Tlist, m_max, w, fa, fbin, ni,epoch_info):  # compute odds-ratios for bins and frequencies
    Om1w = np.zeros((m_max, np.size(w)), 'float')  # odds ratio matrix
    for m in range(0, m_max):
        for wi in range(0, np.size(w)):
            Om1w[m, wi] = compute_Om1(w[wi], Tlist, m + 1, fa[m], fbin, ni,epoch_info)
    return Om1w

def compute_GL(Tlist,epoch_info, m_max=20, w_range=None, ni=10, parallel=False):
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
            w_hi = 5*np.pi * N / T  # max default frequency

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
            dw=w[2]-w[1]
        if w_lo == 0:
            # cannot have w_lo =0
            print("minimum frequency cannot be 0!\n")
            return

        fa = np.zeros(m_max)
        for m in range(0, m_max):  # precompute factors for each m
            fa[m] = compute_factor(N, m + 1, v)
        if parallel:
            Om1w = compute_Om1wPar(Tlist, m_max, w, fa, fbin, ni,epoch_info)
        else:
            Om1w = compute_Om1w(Tlist, m_max, w, fa, fbin, ni,epoch_info)

        pw = 1./ (w*np.log(w_hi / w_lo))
        O1m = np.zeros(m_max)
        for i in range(0, m_max):  # intgreate odd ratios across frequencies
            O1m[i] = np.trapz(pw * Om1w[i], w)

        m_opt = np.argmax(O1m) # find optimum bin number, i.e largest odds-ratio
        S = Om1w[m_opt] / w  # compute Spectral probability
        m_opt = m_opt + 1  # start bin index with 1
        C = np.trapz(S, w)  # compute normalization
        S = S / C  # normalized probability
        O_period = np.sum(O1m[1:])  # integrated odds ratio
        p_period = O_period / (1 + O_period)  # likelihood of periodic event
        cdf = np.array(S)
        # plt.plot(w,cdf)
        # plt.show()
        for i in range(0, np.size(S)):
            cdf[i] = np.trapz(S[0:i], w[0:i])

        wr = np.extract(np.logical_and(cdf > .0025, cdf < .9975), w)
        w_peak = w[np.argmax(S)]
        w_mean = np.trapz(S * w, w)
        if np.size(wr) > 0:
            print('exist')
            w_conf = [np.min(wr), np.max(wr)]
            if np.min(wr)==np.max(wr):
                print('asymmetric')
                w_conf=[np.min(wr)-dw, np.max(wr)+dw]
        else:
            w_conf = [w_peak-dw, w_peak+dw]
        return O_period, p_period, m_opt, S, w, w_peak, w_mean, w_conf,cdf
    else:
        # throw an error
        print("No valid arrival time array provided!\n")
        return O_period, p_period, m_opt, S, w, w_peak, w_mean, w_conf

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
    intN_bin_t_start=np.floor(N_bin_t_start)+1
    intN_bin_t_end=np.floor(N_bin_t_end)

    intN_bin_t_start=intN_bin_t_start.astype(int)
    intN_bin_t_end=intN_bin_t_end.astype(int)
    for i in range(len(N_bin_t_start)):
        if intN_bin_t_end[i]>=intN_bin_t_start[i]:
            T_in_perbin+=int((intN_bin_t_end[i]-intN_bin_t_start[i])/m)*tbin
            T_in_perbin[np.mod(intN_bin_t_start[i],m)-1]+=(intN_bin_t_start[i]-N_bin_t_start[i])*tbin
            T_in_perbin[np.mod(intN_bin_t_end[i],m)]+=(N_bin_t_end[i]-intN_bin_t_end[i])*tbin
            rest=np.mod(intN_bin_t_end[i]-intN_bin_t_start[i],m)
            for k in range(rest):
                T_in_perbin[int(np.mod((intN_bin_t_start[i] + k), m))] += tbin
        else:
            T_in_perbin[np.mod(intN_bin_t_start[i],m)-1]+=(N_bin_t_end[i]-N_bin_t_start[i])*tbin
    return T_in_perbin

#
# path_Tuc = f'/Users/baotong/eSASS/data/raw_data/47_Tuc/txt/txt_merge_psf75_0.2_5/'
path_Tuc='/Users/baotong/Desktop/period_NGC6304/txt_all_obs_p90/'
# path_Tuc='/Users/baotong/Desktop/period_NGC3201/txt_all/txt_all_obs_p90/'
# path_Tuc = f'/Users/baotong/eSASS/data/raw_data/47_Tuc/txt/txt_merge_psf75_0.2_5/'
path = path_Tuc
#path = '/Users/baotong/xmm/M28_LMXB/0701981501/txt/'
# print(sum(get_T_in_mbins(epoch_file,2*np.pi/55000.,10,0.6)))
def filter_obs(src_evt,useid):
    src_evt_use = src_evt[np.where(src_evt[:-1] == useid[0])[0]]
    i=1
    while i < len(useid):
        id=useid[i]
        src_evt_use_temp=src_evt[np.where(src_evt[:-1]==id)[0]]
        src_evt_use = np.concatenate((src_evt_use, src_evt_use_temp))
        i+=1
    return src_evt_use

def write_result(dataname):
    #path = '/Users/baotong/xmm/M28_LMXB/0701981501/txt/'
    data_file=path + str(dataname) + '.txt'
    epoch_file = path + 'epoch_src_{0}.txt'.format(dataname)
    # epoch_file = path + 'epoch_src_' + str(dataname) + '.txt'
    epoch_info=np.loadtxt(epoch_file)
    epoch_info=epoch_info
    if epoch_info.ndim == 1:
        epoch_info=np.array([epoch_info])
    src_evt=np.loadtxt(data_file)
    if src_evt.ndim==1:
        src_evt=np.array([src_evt])

    CR=hawk.plot_longT_V(src_evt=src_evt, bkg_file=None,epoch_info=epoch_info,show=False)
    plt.close()
    print(CR)
    (useid, epoch_info_use)=hawk.choose_obs(epoch_info,flux_info=CR,
                                            flux_filter=100,expT_filter=1000,
                                            if_flux_high=0, if_expT_high=True,obsID=None)
    epoch_info = epoch_info_use  ##这里随意改

    src_evt_use =hawk.filter_obs(src_evt, useid)
    src_evt=src_evt_use

    time=src_evt[:,0]
    energy=src_evt[:,1]
    time = hawk.filter_energy(time, energy, [500, 8000])
    counts=len(time)
    w_range=2*np.pi*np.arange(1./3000,1./100,1.e-5)
    starttime = datetime.datetime.now()
    GL_R=compute_GL(time,epoch_info=epoch_info,w_range=w_range,m_max=20,parallel=True)
    endtime = datetime.datetime.now()
    srcid=dataname
    runtime=(endtime - starttime).seconds
    Prob=GL_R[1]
    wpeak=GL_R[5]
    wmean=GL_R[6]
    mopt=GL_R[2]
    wconf_lo=GL_R[7][0]
    wconf_hi=GL_R[7][1]

    return [srcid,runtime,Prob,wpeak,wmean,mopt,wconf_lo,wconf_hi,counts]


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
    result = np.column_stack((result_srcid, result_runtime, result_Prob, result_wpeak, result_period, result_mopt,
                              result_wconf_lo, result_wconf_hi,result_counts))
    # result = np.column_stack((result_runtime, result_Prob, result_wpeak,result_wmean, result_period, result_mopt,
    #                           result_wconf_lo, result_wconf_hi,result_counts))
    print(result)
    print(path)
    # np.savetxt('result_1h-3h_{0}.txt'.format(id_range[0]), result, fmt='%10d %10.2f %10.2f %10.5f %10.5f %10d %10.5f %10.5f')
    np.savetxt('/Users/baotong/Desktop/period_NGC6304/result_GL/'+'result_3k_{0}.txt'.format(id_range[0]), result,
               fmt='%10d %10.2f %10.5f %10.5f %10.5f %10d %10.5f %10.5f %10d')

if __name__ == '__main__':
    # get_result_fromid(['30'])
    a=np.arange(1,184,1)
    # a=[14,16,22,1,13,20,2,12,51,29,67,17,47,50,70,73,30,64,71,68]
    for i in a:
        get_result_fromid([i])

#choose_id(1, 3)
