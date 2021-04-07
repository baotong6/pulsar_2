import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp
import functools
import datetime
from astropy.io import fits
import sys
import os
import string
#import correct as correct
import random
from scipy import optimize
from sympy import *
from scipy.optimize import fsolve
from itertools import chain
def get_time_series(T_START,N_cts,T_stop,period,amp, model,eqw=0.1):
    #N_cts为总光子数
    #T为总观测时间
    #period为周期
    #amp为model的amplitude,反应了时变幅度
    #model现只针对sin和eclipse函数
    def get_turns(t,period):
        v=1.0/period
        turns=t*v-int(t*v)
        return turns

    if model=='eclipse':
        i=0
        delta=amp
        w=2*np.pi/period
        T=T_stop-T_START
        lam_0=N_cts/T
        t=np.random.random(N_cts)*T+T_START
        t=np.sort(t)
        t_eclipse=t
        while i<len(t_eclipse):
            if (0.5-0.5*eqw)<get_turns(t_eclipse[i],period)<(0.5+0.5*eqw):
                rand=np.random.random(1)[0]
                if rand<amp:
                    t_eclipse=np.delete(t_eclipse,i)
                else:
                    i+=1
            else:i+=1
        return t_eclipse

    if model=='sin':
        delta=amp
        w=2*np.pi/period
        T=T_stop-T_START
        lam_0=N_cts/T
        t=np.random.random(N_cts)*T+T_START
        t=np.sort(t)
        t_sin=t
        i=0
        while i <len(t):
            def f(t1):
                return t1 - delta / w * np.cos(w * t1) - t[i]
            t_sin[i]=fsolve(f,t[i])
            i+=1
        return t_sin

def get_epoch_time_series(cts_rate,period,amp, model):
    t=[]
    # epoch = np.loadtxt('SgrA_I_epoch.txt')
    #epoch = np.loadtxt('LW_const_epoch.txt')
    epoch = np.array([[240853682.79, 241353682.79, 6362, 500000.0], [241353682.79, 241853682.79, 6363, 500000.0]])
    tstart=epoch[:,0]
    tstop=epoch[:,1]
    exp_time_epoch=epoch[:,-1]

    for i in range(len(tstart)):
        cts = int(cts_rate * (tstop[i] - tstart[i]))
        t.append(get_time_series(tstart[i],cts,tstop[i],period,amp,model=model))
    t=list(chain.from_iterable(t))
    t = np.array(t)
    #t = t + np.random.rand(len(t)) * 3.2 - 1.6
    return t

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
    #ln_S_wfi=compute_S(epoch_file,Tlist,w,m,fi)
    #y = np.exp(f + factor+ln_S_wfi)
    y=f+factor
    return y

def get_y(P,w_range,sim_N=100):
    m_max=12;v=m_max-1
    for k in range(sim_N):
        path='/Users/baotong/Desktop/period_LW/simulation/test/'
        epoch_file = 'LW_const_epoch.txt'
        cts_rate_standard = 0.0001
        amp_standard = 1.0
        cts_num=1.0;amp_num=0.7
        time = get_epoch_time_series(cts_rate=cts_rate_standard * cts_num, period=P, amp=amp_standard * amp_num,
                                     model='eclipse')
        N=len(time)
        fa = np.zeros(m_max)
        for m in range(0, m_max):  # precompute factors for each m
            fa[m] = compute_factor(N, m + 1, v)

        fbin = precompute_binmult(N)
        Om1w = np.zeros((m_max, np.size(w_range)), 'float')  # odds ratio matrix
        for m in range(0, m_max):
            for wi in range(0, np.size(w_range)):
                Om1w[m,wi]=compute_W_scaled(time,epoch_file,m+1,w_range[wi],0.5,fa[m],fbin)
        np.savetxt(path+'res_{0}/res_{1}.txt'.format(int(P),k),Om1w,format('%10.5f'))
P1=50000.;P2=5000.;P3=500.
# w1 = 2 * np.pi * np.arange(1. / 55000, 1. / 45000., 1e-9)
# w2 = 2 * np.pi * np.arange(1. / 5500, 1. / 4500., 1e-8)
w3 = 2 * np.pi * np.logspace(1./5000.,1./3000.,1000)
get_y(P=P3,w_range=w3)
# get_y(P=P2,w_range=w2)
# get_y(P=P1,w_range=w1)
def compute_Om1(w, epoch_file,Tlist, m, factor, fbin, ni):
    # compute  specific odds-ratios (eqn. 5.25)
    p = np.arange(0, ni) / float(ni) * 2 * np.pi / m
    #p 即为输入的fi
    # intgration range, only integrate over a single bin, as values of the integral repeat over bins
    y = np.zeros(np.size(p), 'float')  # array to hold values of W_scaled over the integration range
    for i in range(0, np.size(y)):
        y[i] = compute_W_scaled(Tlist,epoch_file, m, w, p[i], factor, fbin)
    #print(np.trapz(y, p) * m)
    #return np.trapz(y, p) * m  # return intregrated W_Scaled
    return np.mean(y)*2*np.pi

def compute_Om1w(Tlist,epoch_file, m_max, w, fa, fbin, ni):  # compute odds-ratios for bins and frequencies
    Om1w = np.zeros((m_max, np.size(w)), 'float')  # odds ratio matrix
    for m in range(0, m_max):
        for wi in range(0, np.size(w)):
            Om1w[m, wi] = compute_Om1(w[wi], epoch_file,Tlist, m + 1, fa[m], fbin, ni)
    return Om1w

def compute_Om1wPar(Tlist, epoch_file,m_max, w, fa, fbin, ni):  # compute odds-ratios for bins and frequencies
    # parallel version
    Om1w = np.zeros((m_max, np.size(w)), 'float')  # odds ratio matrix
    pool = mp.Pool(processes=16)  # use all workers

    for m in range(0, m_max):
        Om1w[m, :] = pool.map(functools.partial(cocompute_Om1wParmpute_Om1, epoch_file=epoch_file,Tlist=Tlist, m=(m + 1), factor=fa[m], fbin=fbin, ni=ni), w)
    pool.close()
    return Om1w

def compute_GL(Tlist,epoch_file, m_max=12, w_range=None, ni=10, parallel=True):
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
                print("error ")
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
        for i in range(0, m_max):
            # intgreate odd ratios across frequencies
            O1m[i] = np.trapz(pw * Om1w[i], w)
        O_period = np.sum(O1m[1:])  # integrated odds ratio
        # print(Om1w[1:])
        # print(O1m[1:])
        # print(O_period)
        p_period = O_period / (1 + O_period)  # likelihood of periodic event
        return [np.mean(Om1w),np.var(Om1w),p_period]

def write_result(P,cts_num,amp_num,w_range):
    epoch_file='LW_const_epoch.txt'
    cts_rate_standard=0.0001
    amp_standard=1.0
    time = get_epoch_time_series(cts_rate = cts_rate_standard*cts_num, period = P, amp =amp_standard*amp_num , model = 'sin')
    GL_R = compute_GL(time, epoch_file,w_range=w_range, m_max=12, parallel=True)
    return GL_R

def compare(N_run):
    path='/Users/baotong/Desktop/period_LW/simulation/test/'
    P1=50000.;P2=5000.;P3=500.
    res1=[];res2=[];res3=[]
    w1 = 2 * np.pi * np.arange(1. / 55000, 1. / 45000., 1e-9)
    w2= 2 * np.pi * np.arange(1. / 5500, 1. / 4500., 1e-8)
    w3 = 2 * np.pi * np.arange(1. / 550, 1. / 450., 1e-7)
    for i in range(N_run):
        res1.append(write_result(P1,1.0,0.7,w1))
        res2.append(write_result(P2,1.0,0.7,w2))
        res3.append(write_result(P3,1.0,0.7,w3))
    np.savetxt(path+'res1.txt',res1)
    np.savetxt(path+'res2.txt', res2)
    np.savetxt(path+'res3.txt', res3)


#compare(20)