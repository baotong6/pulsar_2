import numpy as np
import matplotlib.pyplot as plt
import hawkeye as hawk
import rocket as rocket
import rednoise as rednoise
from scipy import integrate
from scipy.optimize import curve_fit
import scipy.stats
from lmfit import Model
from astropy.modeling import models
from stingray.events import EventList
from stingray.modeling import PSDParEst
from astropy.modeling.fitting import _fitter_to_model_params
from stingray.lightcurve import Lightcurve
from stingray.modeling import PSDLogLikelihood,PSDPosterior
from rednoise import useful_functions as func
import rednoise
from numpy import exp, linspace, random

font1=func.font1
font2=func.font2

def fit_powerspectrum(ps, model, starting_pars, max_post=False, priors=None,
                      fitmethod="L-BFGS-B"):

    if priors:
        lpost = PSDPosterior(ps.freq,ps.power, model, priors=priors,m=ps.m)
    else:
        lpost = PSDLogLikelihood(ps.freq, ps.power, model, m=ps.m)

    parest = PSDParEst(ps, fitmethod=fitmethod, max_post=max_post)
    res = parest.fit(lpost, starting_pars, neg=True)

    return parest, res

def optimize_psdmodel(ps,whitenoise=100,show=0,label='test',save=0,figurepath=None,outfigname=None):
    smoothbkplc=models.SmoothlyBrokenPowerLaw1D()+models.Const1D()
    print(smoothbkplc.param_names)
    smoothbkplc.alpha_1_0.fixed = True
    smoothbkplc.alpha_2_0.fixed = True
    smoothbkplc.delta_0.fixed = True
    smoothbkplc.x_break_0.fixed = True
    # smoothbkplc.amplitude_1.fixed = True

    smoothbkplc.alpha_1_0=1
    smoothbkplc.delta_0=1/4.
    smoothbkplc.alpha_2_0=2.
    smoothbkplc.x_break_0=1
    p_xb = lambda xb: ((1e-5 <= xb) & (xb <= 1e-1))
    p_amplitude = lambda amplitude: ((1e-3 <= amplitude) & (amplitude <= 1000.0))
    p_whitenoise = lambda white_noise: scipy.stats.norm(whitenoise, 0.05).pdf(white_noise)
    priors = {}
    priors["x_break"] = p_xb
    priors["amplitude_0"] = p_amplitude
    # priors["amplitude_1"] = p_whitenoise
    # lpost = PSDPosterior(ps, smoothbkplc, priors=priors)
    parest, res = fit_powerspectrum(ps, model=smoothbkplc, starting_pars=[10,whitenoise], max_post=False, priors=False,
                      fitmethod="BFGS")
    print(res.p_opt)
    print(res.err)
    print(smoothbkplc)
    psd_shape = smoothbkplc(ps.freq)
    plt.title(label,func.font2)
    plt.loglog(ps.freq, ps.power, ds="steps-mid", label="periodogram realization",color='green')
    plt.loglog(ps.freq, res.mfit, label="power spectrum",color='red')
    plt.loglog(ps.freq, np.zeros(len(ps.freq))+whitenoise, label="Poisson noise", color='grey',linestyle='--')
    plt.xlabel('Frequency (Hz)', func.font2)
    plt.ylabel(r'$\rm Power~([rms/mean]^2 Hz^{-1})$', func.font2)
    plt.tick_params(labelsize=16)
    plt.legend(['PSD',r'$\rm Model: P(\nu)=N \nu^{-1}(1+(\frac{\nu}{\nu_0})^4)^{-1/4}+C$','Poisson noise'])
    if save:
        plt.savefig(figurepath + '{0}.pdf'.format(outfigname), bbox_inches='tight', pad_inches=0.0)
    if show:
        plt.show()
    return smoothbkplc

def model_curvefit(ps,ifperiod=0,whitenoise=100,show=0,label='test',save=0,figurepath=None,outfigname=None):
    x=ps.freq;y=ps.power
    print(y)
    def break_po_c(x, delta,alpha,gamma,N,N_p):
        ## by default for CVs in (Revnivtsev+,2010)
        ## p=(3.36e-2Hz,4,-1/4,notsure)
        """
        Parameters
        ----------
        :param x:numpy.ndarray
            non-zero frequencies
        :param p:
        p[0] = break frequency,delta
        p[1] = alpha,即 power law index、
        p[2] = constant,gamma
        p[3] = normalization N or beta
        p[4]= constant (poisson noise)
        :return:
        """
        return N * x ** (-1) * (1 + (x / delta) ** alpha) ** gamma+N_p
    gmodel = Model(break_po_c,nan_policy='raise')
    print(f'parameter names: {gmodel.param_names}')
    print(f'independent variables: {gmodel.independent_vars}')
    print(f'parameter names: {gmodel.param_names}')
    # init_vals = [2e-3,alpha]  # for [amp, cen, wid]
    # best_vals, covar = curve_fit(break_po_c, x, y, p0=init_vals)

    params = gmodel.make_params(delta=2e-3,alpha=4,gamma=-1/4,N=0,N_p=whitenoise)
    params['delta'].vary=False;params['alpha'].vary=False;params['gamma'].vary=False
    gmodel.set_param_hint('N', min=1e-5,max=1000)
    gmodel.set_param_hint('N_p', min=1e-5,max=whitenoise*100)

    gmodel.set_param_hint('delta', min=1e-3,max=1e-2)
    gmodel.set_param_hint('alpha', min=4,max=4+1e-2)
    gmodel.set_param_hint('gamma', min=-1/4,max=-1/4+1e-2)
    result = gmodel.fit(y, x=x, delta=2e-3,alpha=4,gamma=-1/4,N=2,N_p=whitenoise)

    # print(result.fit_report())
    plt.step(x, y,label='PSD')
    # plt.plot(x, result.init_fit, '--', label='initial fit')
    plt.plot(x, result.best_fit, '-',label=r'$\rm Model: P(\nu)=N \nu^{-1}(1+(\frac{\nu}{\nu_0})^4)^{-1/4}+C$')
    # plt.plot([1/ifperiod,1/ifperiod],[0,1000],'->',label='Period')
    plt.annotate("", xy=(1/ifperiod, y.max() * 0.5),
                xytext=(1/ifperiod, y.max()), color="red",
                weight="bold",
                arrowprops=dict(arrowstyle="->", connectionstyle="arc3", color="red"))

    plt.loglog(ps.freq, np.zeros(len(ps.freq))+whitenoise, label="Poisson noise", color='grey',linestyle='--')
    plt.xlabel('Frequency (Hz)', func.font2)
    plt.ylabel(r'$\rm Power~([rms/mean]^2 Hz^{-1})$',font=func.font2)
    plt.tick_params(labelsize=16)
    plt.legend()
    # plt.legend(['PSD',r'$\rm Model: P(\nu)=N \nu^{-1}(1+(\frac{\nu}{\nu_0})^4)^{-1/4}+C$','Poisson noise'])
    plt.loglog()
    plt.legend()
    if save:
        plt.savefig(figurepath + '{0}.pdf'.format(outfigname), bbox_inches='tight', pad_inches=0.0)
    if show:
        plt.show()
    return (gmodel,result)


def test_simple():
    lenbin=100
    T_exp=63468.73
    CR=5e-3
    epoch_89='/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep3/epoch_src_test_single.txt'
    w=rednoise.make_freq_range(dt=lenbin,epoch_info=epoch_89)

    # p_1=[10,2e-3,1,2,1/4]
    p_1 = [2e-3,4,-1/4,5]
    psd_model=rednoise.build_psd(x=w,p=p_1,type='breakp')
    frms = integrate.quad(func.break_po, w[0], w[-1], args=(p_1))[0]
    frms = np.sqrt(frms)
    print(frms)

    plt.loglog()
    plt.plot(w, psd_model)
    plt.xlabel('Frequency (Hz)', font2)
    plt.ylabel('Power', font2)
    plt.tick_params(labelsize=16)
    plt.show()
    lc = rednoise.make_lc_from_psd(psd=psd_model, cts_rate=CR * lenbin, dt=lenbin, epoch_info=epoch_89)
    ps_org = rednoise.plot_psd(lc)
    print('counts={0}'.format(np.sum(lc.counts)))
    lc_evt = Lightcurve(time=lc.time, counts=lc.counts, dt=lc.dt, gti=lc.gti)
    # lc_evt.counts = np.random.poisson(lc_evt.counts)
    ev_all = EventList()
    ev_all.time = func.sim_evtlist(lc)
    lc_evt = ev_all.to_lc(dt=lenbin, tstart=ev_all.time[0] - 0.5 * lenbin,
                          tseg=ev_all.time[-1] - ev_all.time[0])
    ps_real = rednoise.plot_psd(lc_evt)

    # optimize_psdmodel(ps_real, whitenoise=2./CR, show=1, label='test', save=0, figurepath=None, outfigname=None)
    print(len(ps_real.freq))
    print(len(ps_real.power))
    model_curvefit(ps_real.freq,ps_real.power)
    print('white_noise=',2/CR)


def check_fD():
    N=500
    path='/Users/baotong/Desktop/aas/pXS_Tuc/figure/rednoise/232_sim_byconst/'
    list_id=[];list_prob=[];list_period=[]
    for i in range(N):
        filename=f'result_10h_{i+1}.txt'
        res=np.loadtxt(path+filename)
        list_id.append(res[0])
        list_prob.append(res[2])
        list_period.append(res[3])
    list_prob=np.array(list_prob)
    fD=len(np.where(list_prob>0.99)[0])/N
    print(fD)

if __name__=='__main__':
    # (src_evt_use,epoch_info_use)=rednoise.load_data(217)
    # print(epoch_info_use)
    # lc=rednoise.get_hist(t=src_evt_use[:,0],len_bin=50,tstart=epoch_info_use[:,0][0],tstop=epoch_info_use[:,1][-1])
    # psd=rednoise.plot_psd(lc,norm='frac')
    # optimize_psdmodel(psd, whitenoise=100, show=1, label='test', save=0, figurepath=None, outfigname=None)
    # test_simple()
    check_fD()