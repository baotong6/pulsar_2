import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import string
from datetime import datetime
from scipy.interpolate import lagrange
from scipy import interpolate
from scipy.optimize import curve_fit
import pandas as pd
import linecache
from astropy.timeseries import LombScargle
from astropy.stats import poisson_conf_interval
import stingray as sr
from stingray.events import EventList
from stingray.lightcurve import Lightcurve
from stingray import Lightcurve, Crossspectrum, sampledata,Powerspectrum,AveragedPowerspectrum
from stingray.simulator import simulator, models
import useful_functions as func

font1 = {'family': 'Normal',
         'weight': 'normal',
         'size': 18, }

source_id=np.array([495,175,716,479,949,855,291,730,711,938,567,208,220,242,
                    200,19,856,168,979,557,996,89,13,785,485,941,306,643,100,
                    943,28,876,102,988,796,947,507,98,81])

def LSP_band(srcevt,bkgevt,epoch,band,epochnum=[0,31]):
    bin_len=100
    if len(epoch)==0 or len(srcevt)<4:return (1,0)
    if len(np.shape(epoch)) == 1:
        epoch=np.array([epoch])
    if len(bkgevt)<2:bkgevt=[]
    else:
        useid = epoch[:, 2][epochnum[0]:epochnum[1]]
        time=[];bkg_time=[];energy=[];bkg_energy=[]
        srcevt[:, -1] = srcevt[:, -1].astype('int');
        bkgevt[:, -1] = bkgevt[:, -1].astype('int')
        for k in range(len(useid)):
            time =np.concatenate((time, srcevt[:, 0][np.where(srcevt[:, -1] == useid[k])]))
            energy=np.concatenate((energy, srcevt[:, 1][np.where(srcevt[:, -1] == useid[k])]))
            bkg_time = np.concatenate((bkg_time,bkgevt[:, 0][np.where(bkgevt[:, -1] == useid[k])]))
            bkg_energy=np.concatenate((bkg_energy,bkgevt[:, 1][np.where(bkgevt[:, -1] == useid[k])]))
        ## 看一下不同波段的lc如何 ##
        (time,energy)=func.filter_energy(time,energy,band)
        (bkg_time,bkg_energy)=func.filter_energy(bkg_time,bkg_energy,band)
        src_cts=len(time);bkg_cts=len(bkg_time)
        if len(time) < 4:return (1, 0)
        if len(bkg_time) < 2: bkg_time = []
        T_tot=time[-1]-time[0]
        freq=np.arange(1/T_tot,0.5/bin_len-0.0002,1/(5*T_tot))
        freq=freq[np.where(freq > 1 / 20000.)]
        lc=func.get_hist_withbkg(time,bkg_time,bin_len)
        counts=np.sum(lc.counts)
        # cts_rate=(counts/exptime)
        x=lc.time;flux=lc.counts
        (FP, out_period) = func.get_LS(x, flux, freq)
    return (FP,out_period)
def LS_source_among_obs(source_id,k_num):
    path='/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep{0}/'.format(k_num)
    src_evt=np.loadtxt(path+'{0}.txt'.format(source_id))
    bkg_evt=np.loadtxt(path+'{0}_bkg.txt'.format(source_id))
    epoch=np.loadtxt(path+'epoch_src_{0}.txt'.format(source_id))
    epoch_id_num=len(epoch)
    FP=np.zeros((epoch_id_num,epoch_id_num))+1
    out_period=np.zeros((epoch_id_num,epoch_id_num))
    if epoch_id_num==0:print('No observation')
    else:
        for i in range(epoch_id_num):
            for j in range(i,epoch_id_num):
                (FP[i][j],out_period[i][j])=LSP_band(src_evt,bkg_evt,epoch,[500,8000],epochnum=[i,j])
        np.savetxt('/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep{0}/bright_color/FP_obs_color_{1}.txt'.format(k_num,source_id),FP)
        np.savetxt('/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep{0}/bright_color/period_obs_color_{1}.txt'.format(k_num,source_id),out_period)

def plot_color_obs_src(source_id,k_num):
    FP_file='/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep{0}/bright_color/FP_obs_color_{1}.txt'.format(k_num,source_id)
    Period_file='/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep{0}/bright_color/period_obs_color_{1}.txt'.format(k_num,source_id)
    if os.path.exists(FP_file):
        FP = np.loadtxt(FP_file)
        epoch_id_num = len(FP)
        period = np.loadtxt(Period_file)
        im=plt.contourf(np.linspace(1,epoch_id_num,epoch_id_num),np.linspace(1,epoch_id_num,epoch_id_num),FP.T,levels=np.linspace(0,0.01,11),cmap="OrRd_r")
        cand_id=np.where(FP<0.01)
        print(FP[cand_id])
        print(period[cand_id])
        plt.title('XID {0}'.format(source_id))
        plt.colorbar(im)
        plt.show()
    else:
        return None


def plot_candidate_info(source_id,k_num):
    FP_file='/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep{0}/bright_color/FP_obs_color_{1}.txt'.format(k_num,source_id)
    Period_file='/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep{0}/bright_color/period_obs_color_{1}.txt'.format(k_num,source_id)
    if os.path.exists(FP_file):
        FP = np.loadtxt(FP_file)
        epoch_id_num = len(FP)
        period = np.loadtxt(Period_file)
        cand_id=np.where(FP<0.01)
        start=cand_id[0];stop=cand_id[1]
        for i in range(len(start)):
            plt.scatter(start[i],FP[start[i]][stop[i]],color='r')
            plt.scatter(stop[i], FP[start[i]][stop[i]],color='r')
            plt.plot([start[i],stop[i]],[FP[start[i]][stop[i]],FP[start[i]][stop[i]]],color='black')
            plt.text((start[i]+stop[i])/2,FP[start[i]][stop[i]],'P={0}'.format(np.round(period[start[i]][stop[i]],2)))
        plt.semilogy()
        plt.xlim(0,epoch_id_num)
        plt.tick_params(labelsize=16)
        plt.ylabel('FAP',font1)
        plt.xlabel('Obs',font1)
        plt.show()


def select_candidata():
    FP_file='/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep{0}/bright_color/FP_obs_color_{1}.txt'.format(k_num,source_id)
    Period_file='/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep{0}/bright_color/period_obs_color_{1}.txt'.format(k_num,source_id)
    if os.path.exists(FP_file):
        FP=np.loadtxt(FP_file)
        epoch_id_num=len(FP)
        period=np.loadtxt(Period_file)

if __name__=='__main__':
    # LS_source_among_obs(str(source_id[1]), '4')
    # for id in source_id:
    #     for k in [1,2,3]:
    #         LS_source_among_obs(str(id), str(k))
    # LS_source_among_obs('643', '4')
    # for id in source_id:
    #     plot_color_obs_src(id, 4)
    # plot_color_obs_src('89', '3')
    plot_candidate_info('643','4')