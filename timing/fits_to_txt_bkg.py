#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import string
# import correct as correct
from datetime import datetime
from scipy.interpolate import lagrange
from scipy import interpolate
from scipy.fftpack import fft, ifft
import scipy.signal as ss
import random
import pandas as pd
from tkinter import _flatten

obs_ID_I = np.loadtxt('SgrA_I_epoch.txt')[:, 2]
obs_ID_I=obs_ID_I.astype('int')
obs_ID_G = np.loadtxt('SgrA_G_epoch.txt')[:, 2]
obs_ID_G=obs_ID_G.astype('int')

def get_txt(regname,obsid):
    #path = '/Volumes/pulsar/SgrA/merge_data/spectra/'
    path='/Volumes/pulsar/SgrAGRT/merge_data/spectra/'
    evt_list = 'all_bcc_{0}_reproj_evt.fits'.format(obsid)

    reg_file = []
    hdul_evt = fits.open(path + evt_list)
    x=hdul_evt[1].data.field(10)
    y=hdul_evt[1].data.field(11)
    #energy=hdul_evt[1].data.field(14)
    energy = hdul_evt[1].data.field(15)
    time=hdul_evt[1].data.field(0)
    obs_ID = np.array([obsid for i in range(len(time))])
    #
    # x = hdul_evt[1].data.field(0)
    # y = hdul_evt[1].data.field(1)
    # energy = hdul_evt[1].data.field(2)
    # time = hdul_evt[1].data.field(3)
    # obs_ID = hdul_evt[1].data.field(11)

    def read_region_bkg_minus(regname):
        #with open(path + 'spectra_p/'+regname, 'r') as file_to_read:

        with open(path + regname, 'r') as file_to_read:
            while True:
                lines = file_to_read.readline()  # 整行读取数据
                reg_file.append(lines)
                if not lines:
                    break
                    pass
        annulus=reg_file[3][8:-2]
        circle=[]
        i=4
        while reg_file[i]!='':
            if reg_file[i][0]=='-':
                circle.append(reg_file[i][8:-2])
            i+=1
        an_x,an_y,an_r1,an_r2=[float(i) for i in annulus.split(',')]
        an_out=[an_x,an_y,an_r1,an_r2]
        if len(circle)==1:
            ###-----大于1的以后再说吧------###
            cir_x,cir_y,cir_r=[float(i) for i in circle[0].split(',')]
            cir_out=[cir_x,cir_y,cir_r]
            return [an_out,[cir_out]]
        elif len(circle)==0:
            return [an_out,[[]]]

    reg = read_region_bkg_minus(regname)

    def where_region_annulus(x, y, reg):
        print(reg[1])
        if len(reg[1][0])==0:
            an=reg[0]
            r = np.array((x - an[0], y - an[1]))
            len_r = np.sqrt(r[0] ** 2 + r[1] ** 2)
            temp_in= len_r - an[2]
            temp_out=len_r-an[3]

            return np.where((temp_in >= 0)&(temp_out <= 0))
        elif len(reg[1])==1:
            an = reg[0]
            r = np.array((x - an[0], y - an[1]))
            len_r = np.sqrt(r[0] ** 2 + r[1] ** 2)
            temp_in = len_r - an[2]
            temp_out = len_r - an[3]
            temp1=np.where((temp_in >= 0) & (temp_out <= 0))
            print(temp1)

            cir=reg[1][0]
            r = np.array((x - cir[0], y - cir[1]))
            len_r = np.sqrt(r[0] ** 2 + r[1] ** 2)
            temp_cir = len_r - cir[2]
            temp2=np.where(temp_cir <= 0)
            print(temp2)


            return np.setdiff1d(temp1,temp2)


    src_index = where_region_annulus(x, y, reg)
    src_x = x[src_index]
    src_y = y[src_index]
    src_t = time[src_index]
    src_E = energy[src_index]
    src_ID = obs_ID[src_index]

    def delete_photon_ID(time, energy, ID):
        i = 0
        while i < len(energy):
            if energy[i] > 8000 or energy[i] < 2000:
                energy = np.delete(energy, i)
                time = np.delete(time, i)
                ID = np.delete(ID, i)
                i = i - 1
            i = i + 1
        return [time, energy, ID]

    [src_t, src_E, src_ID] = delete_photon_ID(src_t, src_E, src_ID)

    src_t = src_t.astype('float')
    src_E = src_E.astype('float')
    src_ID = src_ID.astype('int')

    src_txt = np.column_stack((src_t, src_E, src_ID))
    src_txt = src_txt[src_txt[:, 0].argsort()]

    path_out= '/Volumes/pulsar/WR/1671/'
    np.savetxt(path_out + 'txt/bkg_txt/' + regname[0:-4] + '_{0}.txt'.format(obsid), src_txt, fmt = "%.7f  %5.3f  %d")


item=1671

# for obsid in obs_ID_I:
#     get_txt(str(item) + '_bkg.reg',obsid)

# for obsid in obs_ID_G:
#     get_txt(str(item) + '_bkg.reg', obsid)

def merge_txt(src_id):
    path='/Volumes/pulsar/WR/1671/txt/'
    epoch_file='epoch_src_{0}_G.txt'.format(src_id)

    epoch_all = np.loadtxt(path + epoch_file)
    obs_tstart=epoch_all[:,0]
    obs_tstop = epoch_all[:,1]
    obs_ID_all=epoch_all[:,2]
    obs_expt=epoch_all[:,-1]

    obs_ID_all=obs_ID_all.astype(int)

    res_t=[]
    res_E=[]
    res_ID=[]
    epoch_ID=[]
    epoch_start=[]
    epoch_stop=[]
    epoch_expt=[]
    for i in range(len(obs_ID_all)):
        res_temp=np.loadtxt(path+'bkg_txt/'+str(src_id)+'_bkg_{0}.txt'.format(obs_ID_all[i]))
        if len(res_temp)==0:
            continue
        elif type(res_temp[0])==type(np.array([1.2])[0]):
            ##判断是否只有1个光子，数据类型的bug,1.2只是随便取的一个值，任意float均可
            epoch_ID.append(obs_ID_all[i])
            epoch_start.append(obs_tstart[i])
            epoch_stop.append(obs_tstop[i])
            epoch_expt.append(obs_expt[i])

            res_t.append(list([res_temp[0]]))
            res_E.append(list([res_temp[1]]))
            res_ID.append(list(([res_temp[2]])))
        else:
            epoch_ID.append(obs_ID_all[i])
            epoch_start.append(obs_tstart[i])
            epoch_stop.append(obs_tstop[i])
            epoch_expt.append(obs_expt[i])

            res_t.append(list(res_temp[:,0]))
            res_E.append(list(res_temp[:,1]))
            res_ID.append(list((res_temp[:,2])))

    res_t=list(_flatten(res_t))
    res_E = list(_flatten(res_E))
    res_ID = list(_flatten(res_ID))
    result=np.column_stack((res_t,res_E,res_ID))
    epoch_info=np.column_stack((epoch_start,epoch_stop,epoch_ID,epoch_expt))
    np.savetxt(path+'epoch_bkg_'+str(src_id)+'_G.txt',epoch_info,fmt='%15.2f %15.2f %10d %20.2f')
    np.savetxt(path+str(src_id)+'_bkg_G.txt',result,fmt="%.7f  %5.3f  %d")
merge_txt('1671')