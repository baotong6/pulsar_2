import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import pandas as pd
import sys
import os
from scipy import interpolate
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from scipy.optimize import curve_fit
import astropy.units as u
import astropy.constants as c
import libwwz
from astropy.timeseries import LombScargle
import stingray as sr
from stingray.events import EventList
from stingray.lightcurve import Lightcurve
from stingray import Lightcurve, Crossspectrum, sampledata,Powerspectrum,AveragedPowerspectrum
from stingray.simulator import simulator, models
import hawkeye as hawk

bin_len=2000.
# path_Tuc='/Users/baotong/Desktop/period_omg/txt_all_obs_p{0}/'.format(ecf)
# path_Tuc='/Users/baotong/Downloads/'
# path_Tuc='/Users/baotong/Desktop/period_NGC3201/txt_all/txt_all_obs_p{0}/'.format(ecf)
# path_Tuc = f'/Users/baotong/eSASS/data/raw_data/47_Tuc/txt/txt_merge_psf{ecf}_0.2_5/'
path = path='/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep2/'
dataname=876
dataname = '{0}.txt'.format(dataname)
epoch_file = path + 'epoch_src_' + dataname
# dataname = 'transient_evt.txt'
# epoch_file = path + 'SgrA_S_epoch.txt'
src_evt=np.loadtxt(path+dataname)
epoch_info=np.loadtxt(epoch_file)
if epoch_info.ndim==1:epoch_info=np.array([epoch_info])

(useid, epoch_info_use)=hawk.choose_obs(epoch_info,flux_info=None,
                                        flux_filter=1,expT_filter=1000,
                                        if_flux_high=0, if_expT_high=1,obsID=[8592])
src_evt_use =hawk.filter_obs(src_evt, useid)
time=hawk.filter_energy(src_evt_use[:,0],src_evt_use[:,1],[500,8000])

lc_new=hawk.get_hist(time,len_bin=bin_len,tstart=epoch_info_use[:,0][0],tstop=epoch_info_use[:,1][-1])
print(useid)
T_tot = epoch_info_use[:, 1][-1] - epoch_info_use[:, 0][0]
print(T_tot)
# freq = np.arange(1 / (0.5*T_tot), 0.5 / bin_len, 1 / (1 * T_tot))
# freq = freq[np.where(freq > 1 / 30000.)]
freq = np.arange(1 / (0.5*T_tot), 0.5 / bin_len, 1 / (1 * T_tot))

# hawk.phase_fold(time=time,epoch_info=epoch_info_use,net_percent=0.99,p_test=7065,outpath=None,bin=15,shift=0.,
#                     label=dataname,text='Seq.{}'.format(dataname),save=0,show=1)
hawk.plot_singleobs_lc(lc_new,period=7065,ifsin=0,figurepath='/Users/baotong/Desktop/aas/pXS_Tuc/figure/',
                       dataname=dataname,save=0,show=1)
# [Tau, Freq, WWZ, AMP, COEF, NEFF]=libwwz.wwt(lc_new.time,lc_new.counts,time_divisions=1000,
#                                              freq_params=[freq[0],freq[-1],freq[-1]-freq[-2],False],decay_constant=1e-3)
# fig, ax = plt.subplots()
# ax.set_title('test title')
# libwwz.plot_methods.linear_plotter(ax,Tau,Freq,WWZ)
# plt.show()