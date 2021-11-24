import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import pandas as pd
import sys
import os
from scipy import interpolate
from scipy.optimize import curve_fit
import astropy.units as u
import astropy.constants as c
from scipy import interpolate
import scipy.stats
import stingray as sr
from stingray.events import EventList
from stingray.lightcurve import Lightcurve
from stingray import Lightcurve, Crossspectrum, sampledata,Powerspectrum,AveragedPowerspectrum
from stingray.simulator import simulator, models
from stingray.modeling import PSDParEst
import matplotlib.font_manager as font_manager
from astropy.timeseries import LombScargle
import warnings
from functools import reduce
from CDFS.CDFS_startover import useful_functions as func
from CDFS.CDFS_startover import sim_psd as sim
from scipy import integrate
from astropy.modeling import models
from astropy.modeling.fitting import _fitter_to_model_params
from stingray.modeling import PSDLogLikelihood,PSDPosterior
def test_somefunc():
    lenbin=100
    CR=1e-3

    epoch_89='/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep3/epoch_src_test_single.txt'
    (TSTART, TSTOP, OBSID, exptime) = func.read_epoch(epoch_89)
    # (TSTART, TSTOP, OBSID, exptime)=([100],[100+100000],[12043],[100000])

    epoch_89=(TSTART, TSTOP, OBSID, exptime)
    w=sim.make_freq_range(dt=lenbin,epoch_file=epoch_89)
    # psd_model=build_psd(x=w,p=[2.3e-3, 3.4, 0., 4.3e-4],x_2=w,
    #                     p_2=[4.01e-4, 4.01e-4 / 16, 100, 2],type='bendp+lorentz')
    p_1=[4.3e-4, 3.4, 0., 2.3e-3];p_2=[1/18000, 16, 0.05, 2]
    psd_model=sim.build_psd(x=w,p=p_1,x_2=w,
                        p_2=p_2,type='bendp+lorentz')
    psd_model=sim.build_psd(x=w,p=p_1,type='bendp')
    frms=integrate.quad(func.bendp_lorentz,w[0],w[-1],args=(p_1,p_2))[0]

    plt.loglog()
    plt.plot(w,psd_model)
    plt.xlabel('Frequency (Hz)',func.font2)
    plt.ylabel('r$ Power [rms/mean]^2 Hz^{-1}$',func.font2)
    plt.tick_params(labelsize=16)
    plt.show()
    lc=sim.make_lc_from_psd(psd=psd_model,cts_rate=CR*lenbin,dt=lenbin,epoch_file=epoch_89,frms=1)
    print(CR*(lc.time[-1]-lc.time[0]))
    ps_org=sim.plot_psd(lc)
    ps_org = Powerspectrum(lc, norm='frac')
    print('counts={0}'.format(np.sum(lc.counts)))
    lc_evt=Lightcurve(time=lc.time,counts=lc.counts,dt=lc.dt,gti=lc.gti)
    lc_evt.counts=np.random.poisson(lc_evt.counts)
    CR=np.sum(lc.counts)/(lc.time[-1]-lc.time[0])
    # ev_all = EventList()
    # ev_all.time = func.sim_evtlist(lc)
    # lc_evt = ev_all.to_lc(dt=lenbin, tstart=ev_all.time[0] - 0.5 * lenbin,
    #                       tseg=ev_all.time[-1] - ev_all.time[0])
    ps_real = sim.plot_psd(lc_evt)
    # freq = np.arange(1 / exptime[0], 0.5 / lenbin, 1 / (5 * exptime[0]))
    # freq = freq[np.where(freq > 1 / 10000.)]
    optimize_psdmodel(ps_real,whitenoise=2/CR)

    # plt.subplot(121)
    plt.xlabel('Time',func.font2)
    plt.ylabel('Counts/bin',func.font2)
    plt.tick_params(labelsize=16)
    plt.plot(lc.time,lc.counts,color='red')
    # plt.subplot(122)
    plt.show()
    plt.plot(lc_evt.time,lc_evt.counts,color='green')
    print('counts={0}'.format(np.sum(lc_evt.counts)))
    plt.xlabel('Time',func.font2)
    plt.ylabel('Counts/bin',func.font2)
    plt.tick_params(labelsize=16)
    plt.show()
    evt=EventList()
    EventList.simulate_times(evt,lc=lc)
    print('counts={0}'.format(len(evt.time)))

def optimize_psdmodel(ps,whitenoise=100,show=0,label='test',save=0,figurepath=None,outfigname=None):
    # define power law component
    pl = models.PowerLaw1D()
    # fix x_0 of power law component
    pl.x_0.fixed = True
    # define constant
    c = models.Const1D()
    # make compound model
    plc = pl + c
    # parest = PSDParEst(ps, fitmethod="L-BFGS-B", max_post=True)
    # print(whitenoise)
    # loglike = PSDLogLikelihood(ps.freq, ps.power, plc, m=ps.m)
    # starting_pars = [1e-8, 2, 1]
    # print(loglike.npar)
    # res = parest.fit(loglike, starting_pars)
    # [amplitude, alpha, white_noise]=res.p_opt
    # print(res.p_opt)
    # _fitter_to_model_params(plc, [amplitude, alpha, white_noise])
    # print(plc)
    # psd_shape = plc(ps.freq)

    # flat prior for the power law index
    p_alpha = lambda alpha: ((0 <= alpha) & (alpha <= 4))

    # flat prior for the power law amplitude
    p_amplitude = lambda amplitude: ((1e-3 <= amplitude) & (amplitude <= 10.0))

    # normal prior for the white noise parameter
    p_whitenoise = lambda white_noise: scipy.stats.norm(whitenoise, 0.05).pdf(white_noise)

    priors = {}
    priors["alpha_0"] = p_alpha
    priors["amplitude_0"] = p_amplitude
    priors["amplitude_1"] = p_whitenoise
    lpost = PSDPosterior(ps.freq, ps.power, plc, priors=priors, m=ps.m)
    # test_pars = [0.002, 0.5,0.29]
    # starting_pars = [0.002, 0.5, 0.29]
    #
    test_pars = [0.002, 0.5,whitenoise]
    starting_pars = [0.002, 0.5,whitenoise]
    print("log-prior: " + str(lpost.logprior(test_pars)))
    print("log-likelihood: " + str(lpost.loglikelihood(test_pars)))
    print("log-posterior: " + str(lpost(test_pars)))
    parest = PSDParEst(ps,  fitmethod="L-BFGS-B", max_post=False)
    res = parest.fit(lpost, starting_pars)
    [amplitude, alpha, white_noise] = res.p_opt
    print(res.err)
    _fitter_to_model_params(plc, [amplitude, alpha, white_noise])
    print(plc)
    psd_shape = plc(ps.freq)
    plt.title(label,func.font2)
    plt.loglog(ps.freq, ps.power, ds="steps-mid", label="periodogram realization",color='green')
    plt.loglog(ps.freq, psd_shape, label="power spectrum",color='red')
    plt.loglog(ps.freq, np.zeros(len(ps.freq))+white_noise, label="Poisson noise", color='grey',linestyle='--')
    plt.xlabel('Frequency (Hz)', func.font2)
    plt.ylabel(r'$\rm Power~([rms/mean]^2 Hz^{-1})$', func.font2)
    plt.tick_params(labelsize=16)
    plt.legend(['PSD',r'$\rm Model: P(\nu)=N \nu^{-\alpha}+C$','Poisson noise'])
    if save:
        plt.savefig(figurepath + '{0}.pdf'.format(outfigname), bbox_inches='tight', pad_inches=0.0)
    if show:
        plt.show()
    return psd_shape

def sim_lc_onesrc(srcid,ep):
    lenbin=100
    path='/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep{0}/'.format(ep)
    epoch_info=path+'epoch_src_{1}.txt'.format(ep,srcid)
    srcevt=np.loadtxt(path+'{0}.txt'.format(srcid))
    bkgevt=np.loadtxt(path+'{0}_bkg.txt'.format(srcid))
    (TSTART, TSTOP, OBSID, exptime) = func.read_epoch(epoch_info)
    # for i in range(len(OBSID)):
    for i in [24]:
        t_src = srcevt[:, 0][np.where(srcevt[:, 2] == OBSID[i])]
        t_bkg = bkgevt[:, 0][np.where(bkgevt[:, 2] == OBSID[i])]
        lc=func.get_hist_withbkg(t_src,t_bkg,len_bin=lenbin,tstart=TSTART[i],tstop=TSTOP[i])
        CR=np.sum(lc.counts)/(TSTOP[i]-TSTART[i])
        print(CR)
        ps_org = Powerspectrum(lc, norm='frac')
        figurepath = '/Users/baotong/Desktop/aas/AGN_CDFS_mod1/figure/'
        ps_shape=optimize_psdmodel(ps_org, whitenoise=2/CR,show=1,label=srcid,save=1,figurepath=figurepath,outfigname='242_psd')

        # CR=(len(t_src)-len(t_bkg)/12.)/exptime[i]
        epoch_temp=([TSTART[i]], [TSTOP[i]], [OBSID[i]], [exptime[i]])
        lc = sim.make_lc_from_psd(psd=ps_shape, cts_rate=CR * lenbin, dt=lenbin, epoch_file=epoch_temp, frms=1)
        lc_evt = Lightcurve(time=lc.time, counts=lc.counts, dt=lc.dt, gti=lc.gti)
        lc_evt.counts = np.random.poisson(lc_evt.counts)

        if i == 0:
            lc_all = lc_evt
        else:
            lc_all = lc_all.join(lc_evt)

    return lc_all

def read_SAS_lc():
    # 1,2,3分别代表mos1,mos2,pn的light curve，也可以加起来用，记为_all;
    # 根据实际情况来决定lomb-scargle的输入

    # 要注意的是，如果在run_XMMproducts_spectra.py中没有自行指定tmin和tmax，这里就不能直接对三个detector的lc进行加减
    dt = 100
    path='/Users/baotong/xmm/0506440101/cal/'
    os.chdir(path)
    mode=['mos1','mos2','pn']
    filename1=mode[0]+'_lccorr_bin{0}.lc'.format(int(dt));filename2=mode[1]+'_lccorr_bin{0}.lc'.format(int(dt));filename3=mode[2]+'_lccorr_bin{0}.lc'.format(int(dt))
    # lc1=fits.open(filename1);
    # lc2=fits.open(filename2);
    lc3=fits.open(filename3)
    # time1=lc1[1].data['TIME'][0:-100];rate1=lc1[1].data['RATE'][0:-100]
    # time2=lc2[1].data['TIME'][0:-100];rate2=lc2[1].data['RATE'][0:-100]
    time3=lc3[1].data['TIME'];rate3=lc3[1].data['RATE']
    rate3=np.nan_to_num(rate3)
    rate3[np.where(rate3<0)]=0
    freq=np.arange(1./10000,0.5/dt,1./(10*100000))

    time_all=time3;rate_all=rate3
    index_gti=np.where(rate_all>0)
    time_all=time_all[index_gti];rate_all=rate_all[index_gti]
    lc=Lightcurve(time=time_all,counts=rate_all)
    print('counts=',np.sum(lc.counts))
    CR = np.sum(lc.counts)*lc.dt / (lc.time[-1] - lc.time[0])
    print('CR= ',CR)
    ps_org = Powerspectrum(lc, norm='frac')
    figurepath = '/Users/baotong/Desktop/aas/AGN_CDFS_mod1/figure/'
    ps_shape = optimize_psdmodel(ps_org, whitenoise=2 / CR, show=1, label='RE J1034+396', save=1, figurepath=figurepath,
                                 outfigname='REJ1034_psd')
if __name__=='__main__':

    # test_somefunc()
    #
    outpath='/Users/baotong/Desktop/CDFS/simulation/sim_src/'
    source_id='242';ep=4
    sim_lc_onesrc(source_id,ep)
    # read_SAS_lc()
    #
    # num_trials=900
    # with open(outpath+'src{0}_ep{1}.txt'.format(source_id,ep),'a+') as f:
    #     k_trial = 0
    #     while k_trial < num_trials:
    #         lc_long=sim_lc_onesrc(source_id,ep=ep)
    #         T_exp=lc_long.time[-1]-lc_long.time[0]
    #         freq = np.arange(1 / T_exp, 0.5 / lc_long.dt, 1 / (5 * T_exp))
    #         freq = freq[np.where(freq > 1 / 20000.)]
    #         temp = func.get_LS(lc_long.time, lc_long.counts, freq=freq)
    #         f.writelines((str(temp[0]) + '        ' + str(temp[1]) + '        ' + str(temp[2]) + '  ' + '\n'))
    #         k_trial+=1
    # sim_lc_onesrc(source_id,ep)