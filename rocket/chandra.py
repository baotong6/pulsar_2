#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import string
#import correct as correct
from datetime import datetime
from scipy.interpolate import lagrange
from scipy import interpolate
import random
import pandas as pd
from tkinter import _flatten
from astropy.wcs import WCS
import warnings

source_info = fits.open('ngc6752_catalog.fits')

ra = source_info[1].data['RAdeg']
dec = source_info[1].data['DEdeg']


def delete_photon_ID(time, energy, ID):

    i = 0
    while i < len(energy):
        if energy[i] > 8000 or energy[i] < 500:
            # if energy[i] > 8000 or energy[i] < 1000:
            energy = np.delete(energy, i)
            time = np.delete(time, i)
            ID = np.delete(ID, i)
            i = i - 1
        i = i + 1
    return [time, energy, ID]

def read_region(regname):
    reg_file=[]
    with open(regname, 'r') as file_to_read:
        while True:
            lines = file_to_read.readline() # 整行读取数据
            reg_file.append(lines)
            if not lines:
                break
                pass
    region=reg_file[-2][7:-2]
    reg_x,reg_y,reg_r=[float(i) for i in region.split(',')]
    return [reg_x,reg_y,reg_r]

def make_region_each_obs(path_in,path_out,ra,dec,wcsimage,obs_ID_all,multiple_src=1,single_name=0):

    os.chdir(path_in)
    fitsname = wcsimage
    w = WCS(path_in + fitsname)
    src_x, src_y = w.all_world2pix(ra, dec, 1)
    binx,biny=np.shape(fits.open(path_in+fitsname)[0].data)
    index=np.union1d(np.where(src_x>binx-1),np.where(src_y>biny-1))
    src_x[index]=binx-1;src_y[index]=biny-1
    phy_x = src_x + 4096-binx/2
    phy_y = src_y + 4096-biny/2
    src_x = np.rint(src_x)
    src_y = np.rint(src_y)
    src_x = src_x.astype(np.int)
    src_y = src_y.astype(np.int)

    for i in range(len(obs_ID_all)):
        os.chdir(path_out)
        os.system('mkdir region_{0}'.format(obs_ID_all[i]))
        p90_list='reproj_psf90_{0}_b5.fits'.format(obs_ID_all[i])
        hdul_p90 = fits.open(path_in + p90_list)

        p90_data = hdul_p90[0].data
        p90_data = p90_data.T
        src_radius = p90_data[src_x, src_y]
        src_radius *= 2.032521
        os.chdir(path_out+'region_{0}'.format(obs_ID_all[i]))
        os.system('mkdir region_90')

        if multiple_src:
            for i in range(len(phy_x)):
                with open('./region_90/{0}.reg'.format(i + 1), 'w+') as f1:
                    reg = 'circle(' + str(phy_x[i]) + ',' + str(phy_y[i]) + ',' + str(src_radius[i]) + ')'
                    f1.writelines(reg)
                with open('./region_90/all.reg', 'a+') as f2:
                    f2.writelines(reg + '\n')
        elif single_name:
            for i in range(len(phy_x)):
                with open('./region_90/{0}.reg'.format(single_name), 'w+') as f1:
                    reg = 'circle(' + str(phy_x[i]) + ',' + str(phy_y[i]) + ',' + str(src_radius[i]) + ')'
                    f1.writelines(reg)
        else:
            print("Please offer the outfile name!")

        return None

def get_txt(path_in,path_out,reg_name,obs_id):

    evt_list='all_bcc_{0}_reproj_evt.fits'.format(obs_id)
    hdul_evt= fits.open(path_in+evt_list)
    x=hdul_evt[1].data['x']
    y=hdul_evt[1].data['y']
    # for non-grating obs, energy in column 14
    #energy = hdul_evt[1].data.field(14)
    # #be careful for grating obs, energey in column 15
    energy=hdul_evt[1].data['energy']
    time=hdul_evt[1].data['time']
    obs_ID=np.array([obs_id for i in range(len(time))])

    def where_region(x,y,reg):
        r=np.array((x-reg[0],y-reg[1]))
        len_r=np.sqrt(r[0]**2+r[1]**2)
        temp=len_r-reg[2]
        return np.where(temp<=0)

    os.chdir(path_out)
    os.system('mkdir txt_{0}'.format(obs_id))

    for item in reg_name:
        reg = read_region(path_in+'region_{0}/'.format(obs_id)+'region_90/'+str(item)+'.reg')
        src_index=where_region(x,y,reg)
        src_x=x[src_index]
        src_y=y[src_index]
        src_t=time[src_index]
        src_E=energy[src_index]
        src_ID=obs_ID[src_index]

        [src_t,src_E,src_ID]=delete_photon_ID(src_t,src_E,src_ID)
        src_t=src_t.astype('float')
        src_E =src_E.astype('float')
        src_ID=src_ID.astype('int')
        src_txt=np.column_stack((src_t,src_E,src_ID))
        src_txt = src_txt[src_txt[:,0].argsort()]

        np.savetxt(path_out+'txt_{0}/'.format(obs_id)+str(item)+'.txt',src_txt,fmt="%.7f  %5.3f  %d")

    return None
