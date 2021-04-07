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
type=['NSC','ND','LW','NSC_G','Tuc']
mode='CDFS'
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
    src_LW =['90','202']
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
    #src_NSC=data.ID_NSC
    #src_NSC = ['1906','2176','2178','2276','3236']
    #src_NSC=['1620','1612','1626','1906','2176','2178','2276','3236']
    src_NSC=['3370','147','1769','2574','1529','1941','1525',
             '6','3483','790','1674','442','3242','350','1748',
             '2268','2478','1854','1873','1080','1083']

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
                ##特别注意，这里的expmap，physical坐标与image坐标一致(忘了是为啥了)
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
if mode=='NSC_G':
    path = '/Volumes/pulsar/SgrAGRT/merge_data/spectra/'
    os.chdir(path)
    obs_ID = np.loadtxt('SgrA_G_epoch.txt')[:, 2]
    os.system('mkdir aprates')
    #src_NSC=data.ID_NSC
    # src_NSC=['2560','1502','2532','1206','3067','1624','2508',
    #          '3357','2841','1219','2672','2422','1853','3120','1133','2730',
    #          '1084','2525','2157','2187','2344','2199','1677','1634','973']
    #src_NSC=['1620','1612','1626','1906','2176','2178','2276','3236']
    src_NSC=['3370','147','1769','2574','1529','1941','1525',
             '6','3483','790','1674','442','3242','350','1748',
             '2268','2478','1854','1873','1080','1083']
    #src_NSC=['1180','1182','1628','1538','2961','1487','1514']

    for ID in src_NSC:
        if str(int(ID))[-3:]=='001' or str(int(ID))[-3:]=='002':
            ID=str(ID)[:-3]
        for obs in obs_ID:
            pi_file_name=str(ID)+'_'+str(int(obs))+'.pi'
            b_pi_file_name = 'b_' + pi_file_name
            expmap_name='reproj_evt2_sou_'+str(int(obs))+'_i4.fits'

            if os.path.exists(pi_file_name) and os.path.exists(b_pi_file_name) and os.path.exists(str(int(ID)) + '.reg'):
                reg = read_region(str(int(ID)) + '.reg')
                pi=fits.open(pi_file_name)
                b_pi=fits.open(b_pi_file_name)
                backscale=pi[1].header['BACKSCAL']
                b_backscale=b_pi[1].header['BACKSCAL']
                backscale*=8192**2
                b_backscale*=8192**2

                pi_data=pi[1].data
                b_pi_data=b_pi[1].data
                # src_cts=sum(pi_data['COUNTS'][137:273])
                # bkg_cts=sum(b_pi_data['COUNTS'][137:273])
                src_cts=sum(pi_data['COUNTS'][137:548])
                bkg_cts=sum(b_pi_data['COUNTS'][137:548])
                print(src_cts,bkg_cts)
                expmap=fits.open(expmap_name)
                exptime_all = expmap[0].data
                exptime_all = exptime_all.T
                img_x=reg[0]-2896.
                img_y=reg[1]-2896.
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

if mode=='Tuc':
    path = '/Volumes/pulsar/omega_cen/merge_data/spectra/'
    os.chdir(path)
    obs_ID = np.loadtxt('omg_cen_epoch.txt')[:, 2]
    os.system('mkdir aprates')
    # src_NSC=data.ID_NSC
    # src_NSC = ['1906','2176','2178','2276','3236']
    # src_NSC=['1620','1612','1626','1906','2176','2178','2276','3236']
    #src_ID = ['38', '121', '5', '53', '7', '97']  #M28
    #src_ID=['12','87']  #omg_cen
    # src_ID=['98','283','282','8','377','299','33','91','74','60','117','84','37','85',
    #         '128','89','78','372','374','66','43','83','375','14','292','352','55','294'] # terzan5
    # src_ID=['245','453','407','182','304','224','206','223','211',
    #         '462','402','294','485','549','258','5','508','520',
    #         '292','350','283','314','284','62','367','345','364']  # 47Tuc
    # src_ID = ['22','92','131','321','320','38','366','134','91','13','85','82','251','326']  #M28
    src_ID = ['69', '52', '36', '87', '25']  #omg_cen
    for ID in src_ID:
        if str(int(ID))[-3:] == '001' or str(int(ID))[-3:] == '002':
            ID = str(ID)[:-3]
        for obs in obs_ID:
            pi_file_name = str(ID) + '_' + str(int(obs)) + '.pi'
            b_pi_file_name = 'b_' + pi_file_name
            expmap_name = 'reproj_evt2_sou_' + str(int(obs)) + '_i5.fits'
            reg = read_region(str(int(ID)) + '.reg')
            if os.path.exists(pi_file_name) and os.path.exists(b_pi_file_name):
                pi = fits.open(pi_file_name)
                b_pi = fits.open(b_pi_file_name)
                backscale = pi[1].header['BACKSCAL']
                b_backscale = b_pi[1].header['BACKSCAL']
                backscale *= 8192 ** 2
                b_backscale *= 8192 ** 2

                pi_data = pi[1].data
                b_pi_data = b_pi[1].data
                src_cts = sum(pi_data['COUNTS'][35:548])
                bkg_cts = sum(b_pi_data['COUNTS'][35:548])
                print(src_cts, bkg_cts)
                expmap = fits.open(expmap_name)
                exptime_all = expmap[0].data
                exptime_all = exptime_all.T
                img_x = reg[0]-2896.
                img_y = reg[1]-2896.
                ##特别注意，这里的expmap，physical坐标与image坐标一致(忘了是为啥了)
                r = reg[2] ** 0.707
                ##为了简便，取个圆的内接正方形对exposure map平均吧
                exp_src = exptime_all[np.arange(int(img_x - r), int(img_x + r))]
                exp_src = exp_src[:, np.arange(int(img_y - r), int(img_y + r))]
                exp_s = np.mean(exp_src)
                aprates_text = 'aprates n={0} m={1} A_s={2} A_b={3} alpha=0.9 beta=0.02 T_s=1 ' \
                               'E_s={4} eng_s=1 flux_s=1 T_b=1 E_b={5} eng_b=1 flux_b=1 clobber=yes ' \
                               'outfile={6} conf=0.9973'.format(src_cts, bkg_cts, backscale, b_backscale, exp_s, exp_s,
                                                                str(int(ID)) + '_' + str(int(obs)) + '_out.par')
                with open('aprates/' + 'run_' + str(int(ID)) + '_' + str(int(obs)) + '.e', 'w+') as f:
                    f.writelines(aprates_text)
            else:
                print('attention')
                continue

if mode=='CDFS':
    path = '/Volumes/pulsar/CDFS/merge_data/spectra/'
    os.chdir(path)
    obs_ID = np.loadtxt('CDFS_epoch.txt')[:, 2]
    os.system('mkdir aprates')
    src_ID = ['711']
    for ID in src_ID:
        for obs in obs_ID:
            pi_file_name = str(ID) + '_' + str(int(obs)) + '.pi'
            b_pi_file_name = 'b_' + pi_file_name
            expmap_name = 'reproj_evt2_sou_' + str(int(obs)) + '_i3.fits'
            reg = read_region(str(int(ID)) + '.reg')
            if os.path.exists(pi_file_name) and os.path.exists(b_pi_file_name):
                pi = fits.open(pi_file_name)
                b_pi = fits.open(b_pi_file_name)
                backscale = pi[1].header['BACKSCAL']
                b_backscale = b_pi[1].header['BACKSCAL']
                backscale *= 8192 ** 2
                b_backscale *= 8192 ** 2

                pi_data = pi[1].data
                b_pi_data = b_pi[1].data
                src_cts = sum(pi_data['COUNTS'][35:548])
                bkg_cts = sum(b_pi_data['COUNTS'][35:548])
                print(src_cts, bkg_cts)
                expmap = fits.open(expmap_name)
                exptime_all = expmap[0].data
                exptime_all = exptime_all.T
                img_x = reg[0]-2896.
                img_y = reg[1]-2896.
                ##特别注意，这里的expmap，physical坐标与image坐标一致(忘了是为啥了)
                r = reg[2] ** 0.707
                ##为了简便，取个圆的内接正方形对exposure map平均吧
                exp_src = exptime_all[np.arange(int(img_x - r), int(img_x + r))]
                exp_src = exp_src[:, np.arange(int(img_y - r), int(img_y + r))]
                exp_s = np.mean(exp_src)
                aprates_text = 'aprates n={0} m={1} A_s={2} A_b={3} alpha=0.9 beta=0.02 T_s=1 ' \
                               'E_s={4} eng_s=1 flux_s=1 T_b=1 E_b={5} eng_b=1 flux_b=1 clobber=yes ' \
                               'outfile={6} conf=0.9973'.format(src_cts, bkg_cts, backscale, b_backscale, exp_s, exp_s,
                                                                str(int(ID)) + '_' + str(int(obs)) + '_out.par')
                with open('aprates/' + 'run_' + str(int(ID)) + '_' + str(int(obs)) + '.e', 'w+') as f:
                    f.writelines(aprates_text)
            else:
                print('attention')
                continue

