#!/bin/bash
# -*- coding: utf-8 -*-
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
import read_csv as data
type=['NSC','ND','LW']
mode='NSC'
def read_region(regname):
    reg_file=[]
    with open(regname, 'r') as file_to_read:
        while True:
            lines = file_to_read.readline()  # 整行读取数据
            reg_file.append(lines)
            if not lines:
                break
                pass
    region = reg_file[-2][7:-2]
    reg_x, reg_y, reg_r = [float(i) for i in region.split(',')]
    return [reg_x, reg_y, reg_r]
if mode=='ND':
    obs_ID=np.loadtxt('ACIS-I_epoch.txt')[:,2]
    path = '/Volumes/pulsar/GC/spectra/'
    os.chdir(path)
    os.system('mkdir aprates')
    src_ID=data.ID_ND
    for ID in src_ID:
        ID=int(ID)
        for obs in obs_ID:
            pi_file_name=str(ID)+'_'+str(int(obs))+'.pi'
            b_pi_file_name = 'b_' + pi_file_name
            expmap_name='reproj_evt2_sou_'+str(int(obs))+'_i4.fits'
            reg = read_region(str(int(ID)) + '.reg')
            if os.path.exists(pi_file_name):
                pi=fits.open(pi_file_name)
                b_pi=fits.open(b_pi_file_name)
                backscale=pi[1].header['BACKSCAL']
                b_backscale=b_pi[1].header['BACKSCAL']
                backscale*=8192**2
                b_backscale*=8192**2

                pi_data=pi[1].data
                b_pi_data=b_pi[1].data
                src_cts=sum(pi_data['COUNTS'][137:548])
                bkg_cts=sum(b_pi_data['COUNTS'][137:548])
                print(src_cts,bkg_cts)
                expmap=fits.open(expmap_name)
                exptime_all = expmap[0].data
                exptime_all = exptime_all.T
                img_x=reg[0]-2896
                img_y=reg[1]-2896
                r=reg[2]**0.707
                ##为了简便，取个圆的内接正方形做平均吧
                exp_src=exptime_all[np.arange(int(img_x-r),int(img_x+r))]
                exp_src=exp_src[:,np.arange(int(img_y-r),int(img_y+r))]
                exp_s=np.mean(exp_src)
                aprates_text='aprates n={0} m={1} A_s={2} A_b={3} alpha=0.9 beta=0.02 T_s=1 ' \
                             'E_s={4} eng_s=1 flux_s=1 T_b=1 E_b={5} eng_b=1 flux_b=1 clobber=yes ' \
                             'outfile={6} conf=0.9973'.format(src_cts,bkg_cts,backscale,b_backscale,exp_s,exp_s,str(int(ID))+'_'+str(int(obs))+'_out.par')
                with open('aprates/'+'run_'+str(int(ID))+'_'+str(int(obs))+'.e', 'w+') as f:
                    f.writelines(aprates_text)
            else:
                print('attention')
                continue

if mode=='LW':
    obs_ID=np.loadtxt('LW_epoch.txt')[:,2]
    path = '/Volumes/pulsar/LimWin_damage/merge_data/spectra'
    os.chdir(path)
    os.system('mkdir aprates')
    src_LW=data.ID_LW
    for ID in src_LW:
        if str(int(ID))[-3:]=='001' or str(int(ID))[-3:]=='002':
            ID=str(ID)[:-3]
        for obs in obs_ID:
            pi_file_name=str(ID)+'_'+str(int(obs))+'.pi'
            b_pi_file_name = 'b_' + pi_file_name
            expmap_name='reproj_evt2_sou_'+str(int(obs))+'_e1.fits'
            reg = read_region(str(int(ID)) + '.reg')
            if os.path.exists(pi_file_name) and os.path.exists(b_pi_file_name):
                pi=fits.open(pi_file_name)
                b_pi=fits.open(b_pi_file_name)
                backscale=pi[1].header['BACKSCAL']
                b_backscale=b_pi[1].header['BACKSCAL']
                backscale*=8192**2
                b_backscale*=8192**2

                pi_data=pi[1].data
                b_pi_data=b_pi[1].data
                src_cts=sum(pi_data['COUNTS'][69:548])
                bkg_cts=sum(b_pi_data['COUNTS'][69:548])
                print(src_cts,bkg_cts)
                expmap=fits.open(expmap_name)
                exptime_all = expmap[0].data
                exptime_all = exptime_all.T
                img_x=reg[0]-2896
                img_y=reg[1]-2896
                r=reg[2]**0.707
                ##为了简便，取个圆的内接正方形做平均吧
                exp_src=exptime_all[np.arange(int(img_x-r),int(img_x+r))]
                exp_src=exp_src[:,np.arange(int(img_y-r),int(img_y+r))]
                exp_s=np.mean(exp_src)
                aprates_text='aprates n={0} m={1} A_s={2} A_b={3} alpha=0.9 beta=0.02 T_s=1 ' \
                             'E_s={4} eng_s=1 flux_s=1 T_b=1 E_b={5} eng_b=1 flux_b=1 clobber=yes ' \
                             'outfile={6} conf=0.9973'.format(src_cts,bkg_cts,backscale,b_backscale,exp_s,exp_s,str(int(ID))+'_'+str(int(obs))+'_out.par')
                with open('aprates/'+'run_'+str(int(ID))+'_'+str(int(obs))+'.e', 'w+') as f:
                    f.writelines(aprates_text)
            else:
                print('attention')
                continue

if mode=='NSC':
    obs_ID=np.loadtxt('SgrA_I_epoch.txt')[:,2]
    path = '/Volumes/pulsar/SgrA/merge_data/spectra/spectra_p'
    os.chdir(path)
    os.system('mkdir aprates')
    src_NSC=data.ID_NSC
    for ID in src_NSC:
        if str(int(ID))[-3:]=='001' or str(int(ID))[-3:]=='002':
            ID=str(ID)[:-3]
        for obs in obs_ID:
            pi_file_name=str(ID)+'_'+str(int(obs))+'.pi'
            b_pi_file_name = 'b_' + pi_file_name
            expmap_name='reproj_evt2_sou_'+str(int(obs))+'_i4.fits'
            reg = read_region(str(int(ID)) + '.reg')
            if os.path.exists(pi_file_name) and os.path.exists(b_pi_file_name):
                pi=fits.open(pi_file_name)
                b_pi=fits.open(b_pi_file_name)
                backscale=pi[1].header['BACKSCAL']
                b_backscale=b_pi[1].header['BACKSCAL']
                backscale*=8192**2
                b_backscale*=8192**2

                pi_data=pi[1].data
                b_pi_data=b_pi[1].data
                src_cts=sum(pi_data['COUNTS'][137:548])
                bkg_cts=sum(b_pi_data['COUNTS'][137:548])
                print(src_cts,bkg_cts)
                expmap=fits.open(expmap_name)
                exptime_all = expmap[0].data
                exptime_all = exptime_all.T
                img_x=reg[0]
                img_y=reg[1]
                r=reg[2]**0.707
                ##为了简便，取个圆的内接正方形对exposure map平均吧
                exp_src=exptime_all[np.arange(int(img_x-r),int(img_x+r))]
                exp_src=exp_src[:,np.arange(int(img_y-r),int(img_y+r))]
                exp_s=np.mean(exp_src)
                aprates_text='aprates n={0} m={1} A_s={2} A_b={3} alpha=0.9 beta=0.02 T_s=1 ' \
                             'E_s={4} eng_s=1 flux_s=1 T_b=1 E_b={5} eng_b=1 flux_b=1 clobber=yes ' \
                             'outfile={6} conf=0.9973'.format(src_cts,bkg_cts,backscale,b_backscale,exp_s,exp_s,str(int(ID))+'_'+str(int(obs))+'_out.par')
                with open('aprates/'+'run_'+str(int(ID))+'_'+str(int(obs))+'.e', 'w+') as f:
                    f.writelines(aprates_text)
            else:
                print('attention')
                continue




