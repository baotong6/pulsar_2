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
if mode=='ND':
    path='/Volumes/pulsar/GC/spectra/aprates'
    obs_ID=np.loadtxt('ACIS-I_epoch.txt')[:,2]
    os.chdir(path)
    src_ID=data.ID_ND
    for ID in src_ID:
        for obs in obs_ID:
            os.system('./run_{0}_{1}.e'.format(str(int(ID)),str(int(obs))))
if mode=='LW':
    path = '/Volumes/pulsar/LimWin_damage/merge_data/spectra/aprates'
    obs_ID=np.loadtxt('LW_epoch.txt')[:,2]
    os.chdir(path)
    src_ID=data.ID_LW
    src_ID=['90','202']
    for ID in src_ID:
        if str(int(ID))[-3:]=='001' or str(int(ID))[-3:]=='002':
            ID=str(ID)[:-3]
        for obs in obs_ID:
            os.system('./run_{0}_{1}.e'.format(str(int(ID)),str(int(obs))))
if mode=='NSC':
    path = '/Volumes/pulsar/SgrA/merge_data/spectra/spectra_p/aprates'
    obs_ID=np.loadtxt('SgrA_I_epoch.txt')[:,2]
    os.chdir(path)
    #src_ID=data.ID_NSC
    #src_ID = ['1180', '1182', '1628', '1538', '2961', '1487', '1514']
    src_ID = ['3370', '147', '1769', '2574', '1529', '1941', '1525',
              '6', '3483', '790', '1674', '442', '3242', '350', '1748',
              '2268', '2478', '1854', '1873', '1080', '1083']
    #src_ID=['1671']
    for ID in src_ID:
        if str(int(ID))[-3:]=='001' or str(int(ID))[-3:]=='002':
            ID=str(ID)[:-3]
        for obs in obs_ID:
            os.system('./run_{0}_{1}.e'.format(str(int(ID)),str(int(obs))))
if mode=='NSC_G':
    path = '/Volumes/pulsar/SgrAGRT/merge_data/spectra/aprates'
    os.chdir(path)
    obs_ID = np.loadtxt('SgrA_G_epoch.txt')[:, 2]
    #src_ID=data.ID_NSC
    #src_ID = ['1180', '1182', '1628', '1538', '2961', '1487', '1514']

    # src_ID = ['2560','1502','2532','1206','3067','1624','2508',
    #          '3357','2841','1219','2672','2422','1853','3120','1133','2730',
    #          '1084','2525','2157','2187','2344','2199','1677','1634','973']
    src_ID = ['3370', '147', '1769', '2574', '1529', '1941', '1525',
               '6', '3483', '790', '1674', '442', '3242', '350', '1748',
               '2268', '2478', '1854', '1873', '1080', '1083']
    for ID in src_ID:
        if str(int(ID))[-3:]=='001' or str(int(ID))[-3:]=='002':
            ID=str(ID)[:-3]
        for obs in obs_ID:
            os.system('./run_{0}_{1}.e'.format(str(int(ID)),str(int(obs))))
if mode=='Tuc':
    path = '/Volumes/pulsar/omega_cen/merge_data/spectra/aprates'
    os.chdir(path)
    #os.system("cp ../47Tuc_epoch.txt ./")
    obs_ID = np.loadtxt('omg_cen_epoch.txt')[:, 2]
    # src_ID=data.ID_NSC
    # src_ID = ['1180', '1182', '1628', '1538', '2961', '1487', '1514']

    # src_ID = ['2560','1502','2532','1206','3067','1624','2508',
    #          '3357','2841','1219','2672','2422','1853','3120','1133','2730',
    #          '1084','2525','2157','2187','2344','2199','1677','1634','973']
    # src_ID = ['3370', '147', '1769', '2574', '1529', '1941', '1525',
    #           '6', '3483', '790', '1674', '442', '3242', '350', '1748',
    #           '2268', '2478', '1854', '1873', '1080', '1083']
    #src_ID=['12','87'] #omg_cen
    # src_ID = ['98', '283', '282', '8', '377', '299', '33', '91', '74', '60', '117', '84', '37', '85',
    #           '128', '89', '78', '372', '374', '66', '43', '83', '375', '14', '292', '352', '55', '294']  # terzan5
    # src_ID=['245','453','407','182','304','224','206','223','211',
    #         '462','402','294','485','549','258','5','508','520',
    #         '292','350','283','314','284','62','367','345','364']  # 47Tuc

    # src_ID = ['22', '92', '131', '321', '320', '38', '366', '134', '91', '13', '85', '82', '251', '326']  # M28
    src_ID = ['69', '52', '36', '87', '25']  #omg_cen
    for ID in src_ID:
        if str(int(ID))[-3:] == '001' or str(int(ID))[-3:] == '002':
            ID = str(ID)[:-3]
        for obs in obs_ID:
            os.system('./run_{0}_{1}.e'.format(str(int(ID)), str(int(obs))))
if mode=='CDFS':
    path = '/Volumes/pulsar/CDFS/merge_data/spectra/'
    os.chdir(path+'aprates/')
    #os.system("cp ../47Tuc_epoch.txt ./")
    obs_ID = np.loadtxt(path+'CDFS_epoch.txt')[:, 2]
    src_ID = ['711']
    for ID in src_ID:
        for obs in obs_ID:
            os.system('./run_{0}_{1}.e'.format(str(int(ID)), str(int(obs))))