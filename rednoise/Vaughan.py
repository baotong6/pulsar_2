#!/bin/bash
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 19 18:13:40 2022
@author: wafels
modified by Tong
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.stats import poisson_conf_interval
from scipy import interpolate
from scipy.optimize import curve_fit
import astropy.units as u
import astropy.constants as c
import stingray as sr
from stingray.events import EventList
from stingray.lightcurve import Lightcurve
from stingray import Lightcurve, Crossspectrum, sampledata,Powerspectrum,AveragedPowerspectrum
from stingray.simulator import simulator, models
import rednoise
import emcee
from astropy.io import fits
import corner
from rednoise import useful_functions as func

font1=func.font1
font2=func.font2

def func_forCV(x,p):
    """One dimensional smoothly broken power law model.

    Parameters
    ----------
    amplitude : float
        Model amplitude at the break point.
    x_break : float
        Break point.
    Cp : float
        Poisson noise
    """
    # p[1]=1e-3
    return p[0]*x**(-1)*(1+(x/3e-2)**4)**(-1/4.)+p[1]

def powerlaw(x,p):
    """One dimensional smoothly broken power law model.
    Parameters
    ----------
    amplitude : float
        Model amplitude at the break point.
    x_break : float
        Break point.
    Cp : float
        Poisson noise
    """
    # p[1]=1e-3

    return p[0]*x**(-1)+p[1]

def vaughan_2010_T_ISS(iobs, S):
    """Vaughan, 2010, MNRAS, 402, 307. Eq. 17.
    Returns a test statistic measuring the difference between
    the observed power spectrum iobs and another power spectrum S."""
    out=np.sum(iobs / S+np.log(S))
    if np.isnan(out):out=1e8
    return -out

def compute_sigma_level(trace1, trace2, nbins = 20):
    """From a set of traces, bin by number of standard deviations"""
    L, xbins, ybins = np.histogram2d(trace1, trace2, nbins)
    L[L == 0] = 1E-16
    logL = np.log(L)

    shape = L.shape
    L = L.ravel()

    # obtain the indices to sort and unsort the flattened array
    i_sort = np.argsort(L)[::-1]
    i_unsort = np.argsort(i_sort)

    L_cumsum = L[i_sort].cumsum()
    L_cumsum /= L_cumsum[-1]

    xbins = 0.5 * (xbins[1:] + xbins[:-1])
    ybins = 0.5 * (ybins[1:] + ybins[:-1])

    return xbins, ybins, L_cumsum[i_unsort].reshape(shape)

def plot_MCMC_trace(ax, xdata, ydata, trace, scatter = False, **kwargs):
    """Plot traces and contours"""
    xbins, ybins, sigma = compute_sigma_level(trace[0], trace[1])
    ax.contour(xbins, ybins, sigma.T, levels = [0.683, 0.955,0.9973], **kwargs)
    if scatter:
        ax.plot(trace[0], trace[1], ',k', alpha = 0.1)
    ax.set_xlabel('a0')
    ax.set_ylabel('a1')

def plot_MCMC_trace2(ax, xdata, ydata, trace, scatter = False, **kwargs):
    xbins2, ybins2, sigma2 = compute_sigma_level(trace[1], trace[2])
    ax.contour(xbins2, ybins2, sigma2.T, levels = [0.683, 0.955, 0.9973], **kwargs)
    if scatter:
        ax.plot(trace[1], trace[2], ',k', alpha = 0.1)
    ax.set_xlabel('a1')
    ax.set_ylabel('a2')

def plot_MCMC_model(ax, xdata, ydata, trace,CR=None,show=0):
    """Plot the linear model and 2sigma contours"""
    # ax.plot(xdata, ydata, 'ok')
    ax.step(xdata, ydata,label='PSD')

    p0, p1= trace[:2]
    xfit = xdata
    # yfit = p0[:,None]*xdata**(-1)*(1+(xdata/1e-2)**4)**(-1/4.)+p1[:,None]
    yfit = p0[:, None] * xdata ** (-p1[:,None]) +2/CR
    mu = yfit.mean(0)
    sig = 1 * yfit.std(0)

    ax.plot(xfit, mu, '-k',label=r'$\rm Model: P(\nu)=N \nu^{-1}(1+(\frac{\nu}{\nu_0})^4)^{-1/4}+C_p$')
    ax.fill_between(xfit, mu - sig, mu + sig, color = 'lightgray')
    ax.plot(xfit,np.zeros(len(xfit))+2/CR,'--',color='orange',label='Poisson noise')
    # ax.plot(np.zeros(len(yfit))+2012.0,yfit,'--',color='yellow')
    ax.set_xlabel('Frequency (Hz)', font2)
    ax.set_ylabel(r'$\rm Power~([rms/mean]^2 Hz^{-1})$', font2)
    ax.legend()
    ax.tick_params(labelsize=16)
    ax.loglog()
    figurepath = '/Users/baotong/Desktop/aas/pXS_Tuc_mod1/figure/rednoise/'
    # plt.savefig(figurepath + '{0}.pdf'.format('312_2738_psd'), bbox_inches='tight', pad_inches=0.0)
    if show:plt.show()
    else:plt.close()

def plot_MCMC_results(xdata, ydata, trace, colors = 'k',CR=None,show=0):
    """Plot both the trace and the model together"""
    # fig, ax = plt.subplots(1, 2, figsize = (10, 5))
    fig,ax=plt.subplots(1,1,figsize=(8,7))
    # plot_MCMC_trace(ax[0], xdata, ydata, trace, True, colors = colors)
    # plot_MCMC_trace2(ax[1], xdata, ydata, trace, True, colors = colors)
    # plot_MCMC_model(ax[1], xdata, ydata, trace,CR)
    plot_MCMC_model(ax, xdata, ydata, trace, CR,show=show)
    # plt.savefig('DEC_3b.eps')


# Define our posterior using Python functions
# for clarity, I've separated-out the prior and likelihood
# but this is not necessary. Note that emcee requires log-posterior

# 用Python函数定义后验，我把先验和似然估计分开写了，其实没必要，主要是显得更简洁。
# 注意emcee需要对数后验证

def log_prior(p):
    a0, a1= p
    if (0< a0 < 1e-1  and a1>0) :
        return 0.0  # log(0)
    else:
        #return -1.5 * np.log(1 + a1 ** 2) - np.log(sigma)
        return -np.inf

def log_likelihood(func,p,x, y,yerr=None):
    sigma = yerr
    y_model = func(x,p)
    print(p)
    if np.min(y_model)<0:
        print('p=',p[np.argmin(y_model)])
        print('x=',x[np.argmin(y_model)])
    # sigma = yerr ** 2 + y_model** 2 * np.exp(2 * log_f)
    S=vaughan_2010_T_ISS(iobs=y,S=y_model)
    return S

def log_posterior(p, x, y,model,yerr=None):

    return log_prior(p) + log_likelihood(model,p, x, y,yerr)

def mcmcfit(xdata,ydata,model,yerr=None,CR=None,show=0):
    # Here we'll set up the computation. emcee combines multiple "walkers",
    # each of which is its own MCMC chain. The number of trace results will
    # be nwalkers * nsteps
    # 设置计算参数。emcee组合了多个"walkers"，每个都有自己的MCMC链。
    # 跟踪结果的数量为 nwalkers * nsteps
    ndim = 2  # number of parameters in the model
    nwalkers = 100  # number of MCMC walkers
    nburn = 100  # "burn-in" period to let chains stabilize
    nsteps = 50  # number of MCMC steps to take
    # set theta near the maximum likelihood, with
    np.random.seed(0)
    # starting_guesses = np.random.random((nwalkers, ndim))
    starting_guesses=np.zeros((nwalkers,ndim))
    starting_guesses[:,0]=np.zeros(nwalkers)+1e-2+np.random.random(nwalkers)
    starting_guesses[:,1]=2/CR+np.random.random(nwalkers)
    # print(starting_guesses)
    # starting_guesses[:,1]=np.random.random(nwalkers)
    # print(starting_guesses)
    # Here's the function call where all the work happens:
    # we'll time it using IPython's %time magic
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, args=[xdata, ydata, model,yerr])
    sampler.run_mcmc(starting_guesses, nsteps)
    print("done")
    # sampler.chain is of shape (nwalkers, nsteps, ndim)
    # we'll throw-out the burn-in points and reshape:
    # sampler.chain返回数组维度为(nwalkers, nsteps, ndim)
    sampler.chain
    emcee_trace = sampler.chain[:, nburn:, :].reshape(-1, ndim).T
    plot_MCMC_results(xdata, ydata, emcee_trace,CR=CR,show=show)
    flat_samples = sampler.get_chain(discard=100, thin=15, flat=True)

    # corner.corner(flat_samples,truths=[0.01,2/CR])
    a=corner.corner(flat_samples,quantiles=[0.16, 0.5, 0.84],show_titles=True,title_kwargs={"fontsize": 12})
    mu1=corner.quantile(flat_samples[:,0],q=[0.5])[0]
    mu2=corner.quantile(flat_samples[:,1],q=[0.5])[0]
    if show:
        plt.show()
    else:
        plt.close()

    return (mu1,mu2)
def gogogo():
    # from timing_comb import get_lc_frombkgimg
    # ecf=75
    # path = f'/Users/baotong/eSASS/data/raw_data/47_Tuc/txt/txt_merge_psf{ecf}_0.2_5/'
    # (src_evt_use,epoch_info_use)=rednoise.singleobs_psd.load_data(dataname='481',ecf=90,ifobsid=[700011,700163,700013,700014],ifexpT=0,path_provide=path)

    # t_sim = rednoise.simulate_const(src_evt_use, epoch_info_use)
    # expT = np.sum(epoch_info_use[:, 3])
    # # lc= rednoise.get_hist_withbkg(t=src_evt_use[:, 0], len_bin=100, tstart=epoch_info_use[:, 0][0], tstop=epoch_info_use[:, 1][-1])
    # lc = get_lc_frombkgimg(int('481'), src_evt_use, epoch_info_use, ecf=75, bin_len=100)
    # lc.gti=[[lc.gti[0][0],lc.gti[-1][-1]]]
    # CR = np.mean(lc.counts) / lc.dt
    # psd = rednoise.plot_psd(lc, norm='frac', show=1, ifexpTfilter=expT)
    id=229
    (src_evt_use,epoch_info_use)=rednoise.singleobs_psd.load_data(id,ifobsid=[2737])
    expT=np.sum(epoch_info_use[:,3])
    lc = rednoise.get_hist(t=src_evt_use[:, 0], len_bin=100, tstart=epoch_info_use[:, 0][0], tstop=epoch_info_use[:, 1][-1])
    CR = np.mean(lc.counts) / lc.dt
    psd = rednoise.plot_psd(lc, norm='frac', show=0, ifexpTfilter=expT)
    xdata=np.array(psd.freq);ydata=np.array(psd.power)
    maskfreq=0
    if maskfreq:
        mask_index = np.argmin(np.abs(xdata - maskfreq))
        y_mask = np.concatenate((ydata[0:mask_index - 1], ydata[mask_index + 1:]))
        x_mask = np.concatenate((xdata[0:mask_index - 1], xdata[mask_index + 1:]))
        x_mask = np.concatenate((xdata[0:mask_index - 1], xdata[mask_index + 1:]))
        result_mu=mcmcfit(x_mask,y_mask,CR=CR)
    else:
        result_mu=mcmcfit(xdata,ydata,CR=CR)

    print(result_mu)
    return result_mu

def apply_Vaughan(lc,epoch_info,model,maskfreq=0,show=0):
    if epoch_info.ndim==1:
        epoch_info=np.array([epoch_info])
        expT=float(epoch_info[:,3])
        print('expT=',expT)
    else:expT = np.sum(epoch_info[:, 3])
    CR = np.mean(lc.counts) / lc.dt
    # print(CR)
    psd = rednoise.plot_psd(lc, norm='frac', show=1, ifexpTfilter=expT)
    print(psd.df)
    xdata=np.array(psd.freq);ydata=np.array(psd.power)
    if maskfreq:
        mask_index = np.argmin(np.abs(xdata - maskfreq))
        y_mask = np.concatenate((ydata[0:mask_index - 1], ydata[mask_index + 1:]))
        x_mask = np.concatenate((xdata[0:mask_index - 1], xdata[mask_index + 1:]))
        x_mask = np.concatenate((xdata[0:mask_index - 1], xdata[mask_index + 1:]))
        result_mu=mcmcfit(x_mask,y_mask,model,CR=CR,show=show)
    else:
        result_mu=mcmcfit(xdata,ydata,model,CR=CR,show=show)
        print('sadasd')
    sigma_range=poisson_conf_interval(lc.counts)
    sigma=(sigma_range[1:,]-sigma_range[0,:])/2
    mse=np.mean(sigma**2)
    frac_rms=np.sqrt((np.var(lc.counts) -mse) /np.mean(lc.counts) ** 2)
    print('frac rms=', frac_rms)
    print('result_mu:', result_mu)
    return (result_mu,psd)

if __name__=='__main__':
    id = 148
    # (src_evt_use, epoch_info_use) = rednoise.singleobs_psd.load_data(id,path_provide='/Users/baotong/Desktop/period_NGC6397/txt_all_obs_p90/', ifobsid=[7460])
    # expT = np.sum(epoch_info_use[:, 3])
    # lc = rednoise.get_hist(t=src_evt_use[:, 0], len_bin=500, tstart=epoch_info_use[:, 0][0],
    #                        tstop=epoch_info_use[:, 1][-1])
    # print('2/CR=',2/np.sum(lc.counts)*expT)
    # apply_Vaughan(lc, epoch_info=epoch_info_use, model=powerlaw, maskfreq=None)
    path = '/Users/baotong/Desktop/XMMcentral/all_lc/'
    a = fits.open(path + '0801681301_267.32_28.56_4690_pnsrc_2_10keV.lc')
    rate = a[1].data['RATE']
    time = a[1].data['TIME']
    T = time[-1] - time[0]
    lc = Lightcurve(time=time, counts=rate * (time[1] - time[0]))
    epoch = np.array([time[0], time[-1],11111, T])
    apply_Vaughan(lc, epoch_info=epoch, model=powerlaw, maskfreq=None)
