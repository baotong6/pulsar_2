#!/bin/bash
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 19 18:13:40 2022
@author: baotong
The main part a test version for getting single source PSD
the function "load_data" is very important, notice for the parameters
"""

import numpy as np
import matplotlib.pyplot as plt
import hawkeye as hawk
from scipy import signal
import rocket as rocket
import rednoise as rednoise
from scipy import integrate
import stingray as sr
from astropy.modeling import models
from stingray.modeling import PSDParEst
from astropy.modeling.fitting import _fitter_to_model_params
from stingray.modeling import PSDLogLikelihood,PSDPosterior

path='/Users/baotong/Desktop/period_Tuc/txt_startover/txt_all_obs_p90/'

def load_data(dataname,ecf=90,ifobsid=None,ifexpT=0):
    path = '/Users/baotong/Desktop/period_Tuc/txt_startover/txt_all_obs_p90/'
    dataname = '{0}.txt'.format(dataname)
    epoch_file = path + 'epoch_src_' + dataname
    # dataname = 'transient_evt.txt'
    # epoch_file = path + 'SgrA_S_epoch.txt'
    src_evt=np.loadtxt(path+dataname)
    epoch_info=np.loadtxt(epoch_file)
    if epoch_info.ndim==1:epoch_info=np.array([epoch_info])
    CR=hawk.plot_longT_V(src_evt=src_evt, bkg_file=None,epoch_info=epoch_info)
    plt.close()
    # print(epoch_info[:,2])
    CR/=ecf/100.

    (useid, epoch_info_use)=hawk.choose_obs(epoch_info,flux_info=CR,
                                            flux_filter=2,expT_filter=ifexpT,
                                            if_flux_high=0, if_expT_high=1,obsID=ifobsid)
    # [78,953,955,956, 2735, 3384, 2736, 3385, 2737, 3386, 2738, 3387,16527,15747,16529,17420,15748,16528]
    src_evt_use =hawk.filter_obs(src_evt, useid)
    # print(useid)
    return (src_evt_use,epoch_info_use)


if __name__=='__main__':
    figurepath='/Users/baotong/Desktop/aas/pXS_Tuc/figure/rednoise/'
    idlist=[217,185,414,366,423,263,331,273,317,162,252,232,283,290,198,312,229]
    period=[205.02498,8517.88756,8646.77907,10311.93607,
            10384.21599,11485.01206,15232.63368,16824.25385,
            23583.10595,31200.2746,44642.85714,44883.3034,45310.3761,
            46082.9493,46151.00609,48780.49,95731.2]
    for i in range(len(idlist)):
    # for i in [14]:
        id=idlist[i]
        (src_evt_use,epoch_info_use)=load_data(id,ifobsid=[2738])
        t_sim=rednoise.simulate_const(src_evt_use,epoch_info_use)
        expT=np.sum(epoch_info_use[:,3])
        # lc=rednoise.get_longhist_withgti(src_evt_use,epoch_info_use,len_bin=2000)
        # lc.apply_gtis()
        # print(lc.counts)
        # print(lc.counts_err)
        # (frms, frms_err) = sr.excess_variance(lc, normalization='fvar')
        # print('frms=', frms)

        lc=rednoise.get_hist(t=src_evt_use[:,0],len_bin=100,tstart=epoch_info_use[:,0][0],tstop=epoch_info_use[:,1][-1])
        CR=np.mean(lc.counts)/lc.dt
        psd=rednoise.plot_psd(lc,norm='frac',show=0,ifexpTfilter=expT)
        # print(psd.df)
        # psd.rebin(df=1/expT)
        plt.close()
        # smoothbkplc=rednoise.optimize_psdmodel(psd, whitenoise=2/CR, show=1, label=str(id)+'_2737', save=1, figurepath=figurepath, outfigname=str(id)+'_2737')
        (model,fitresult)=rednoise.model_curvefit(psd,ifperiod=period[i], whitenoise=2/CR,
                                                  label=str(id)+'_2738', show=1, save=1, figurepath=figurepath, outfigname=str(id)+'_2738')
        print('CR=',CR)
        # print(fitresult.param_names)
        np.savetxt(figurepath+str(id)+'_2738.txt',np.array(list(fitresult.best_values.values())))

        # param=open(figurepath+str(id)+'_2738.txt','w+')
        # paramnames=[f'   {x}   'for x in smoothbkplc.param_names]
        # paramnums=[f'    {x}  'for x in smoothbkplc.parameters ]
        # param.writelines(paramnames)
        # param.writelines('\n')
        # param.writelines(paramnums)