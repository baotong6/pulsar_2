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
warnings.filterwarnings("ignore")
path='/Volumes/pulsar/CDFS/merge_data/timing/'
obs_ID=[581,1431,441,582,2406,2405,2312,1672,2409,2313,2239,8591,9593,9718,8593,8597,8595,8592,8596,
9575,9578,8594,9596,12043,12123,12044,12128,12045,12129,12135,12046,12047,12137,12138,12055,12213,12048,
12049,12050,12222,12219,12051,12218,12223,12052,12220,12053,12054,12230,12231,12227,12233,12232,12234,16183,
16180,16456,16641,16457,16644,16463,17417,17416,16454,16176,16175,16178,16177,16620,16462,17535,17542,16184,
16182,16181,17546,16186,16187,16188,16450,16190,16189,17556,16179,17573,17633,17634,16453,16451,16461,16191,
16460,16459,17552,16455,16458,17677,18709,18719,16452,18730,16185]
epoch_file='CDFS_epoch.txt'
obs_ID_all=np.loadtxt(path+epoch_file)[:,2]
obs_ID_all=obs_ID_all.astype(int)
obs_ID_all=obs_ID_all.astype(str)
def make_region_each_obs():
    os.chdir(path)
    fitsname = 'CDFS_3.fits'
    w = WCS(path + fitsname)
    source_info = fits.open('7Ms_catalog.fit')
    ra = source_info[1].data['RAJ2000']
    dec = source_info[1].data['DEJ2000']
    src_x, src_y = w.all_world2pix(ra, dec, 1)
    phy_x = src_x + 2896
    phy_y = src_y + 2896
    #print(phy_x,phy_y)

    src_x = np.rint(src_x)
    src_y = np.rint(src_y)
    src_x = src_x.astype(np.int)
    src_y = src_y.astype(np.int)

    src_x[np.where(src_x > 2399)] = 2399
    src_x[np.where(src_x < 0)] = 0
    src_y[np.where(src_y > 2399)] = 2399
    src_y[np.where(src_y < 0)] = 0

    for i in range(len(obs_ID_all)):
        os.chdir(path)
        os.system('mkdir region_{0}'.format(obs_ID_all[i]))
        p90_list = 'reproj_psf90_{0}_b3.fits'.format(obs_ID_all[i])
        hdul_p90 = fits.open(path + p90_list)

        p90_data = hdul_p90[0].data
        p90_data = p90_data.T
        src_radius = p90_data[src_x, src_y]
        # src_radius *= 2.032521*0.6
        # src_radius=[1.2*2.032521]
        os.chdir(path + 'region_{0}'.format(obs_ID_all[i]))
        os.system('mkdir region_90')

        for i in range(len(phy_x)):
            with open('./region_90/{0}_bkg.reg'.format(i + 1), 'w+') as f1:
                reg = 'annulus(' + str(phy_x[i]) + ',' + str(phy_y[i]) + ',' + str(2*src_radius[i]) + ','+str(4*src_radius[i])+')'
                f1.writelines(reg)
            with open('./region_90/all_bkg.reg', 'a+') as f2:
                f2.writelines(reg + '\n')
# make_region_each_obs()
def get_txt(obs_id):
    # for CDFS
    source_id = np.linspace(1, 1055, 1055)
    source_id = source_id.astype(int)

    os.chdir(path)
    os.system('mkdir txt_{0}_0.5_8'.format(obs_id))

    evt_list='all_bcc_{0}_reproj_evt.fits'.format(obs_id)
    hdul_evt= fits.open(path+evt_list)
    x=hdul_evt[1].data['x']
    y=hdul_evt[1].data['y']
    energy=hdul_evt[1].data['energy']
    time=hdul_evt[1].data['time']
    obs_ID=np.array([obs_id for i in range(len(time))])

    def read_region_bkg(regname):
        reg_file = []
        with open(path+'region_{0}/'.format(obs_id)+'region_90/'+regname, 'r') as file_to_read:
            while True:
                lines = file_to_read.readline() # 整行读取数据
                reg_file.append(lines)
                if not lines:
                    break
                    pass
        annulus = reg_file[0][8:-2]
        an_x, an_y, an_r1, an_r2 = [float(i) for i in annulus.split(',')]
        an_out = [an_x, an_y, an_r1, an_r2]
        return an_out

    def where_region_annulus(x, y, reg):
        an=reg
        r = np.array([x - an[0], y - an[1]])
        len_r = np.sqrt(r[0] ** 2 + r[1] ** 2)
        temp_in= len_r - an[2]
        temp_out=len_r-an[3]
        return np.where((temp_in >= 0)&(temp_out <= 0))

    def delete_photon_energy(x,y,time,energy,ID):
        i=0
        while i < len(energy):
            if energy[i]>8000 or energy[i]<500:
            #if energy[i] > 8000 or energy[i] < 2000:
                energy=np.delete(energy,i)
                time=np.delete(time,i)
                ID=np.delete(ID,i)
                x=np.delete(x,i);y=np.delete(y,i)
                i=i-1
            i=i+1
        return [x,y,time,energy,ID]

    def delete_photon_overlap(x,y,time,energy,ID):
        path = '/Volumes/pulsar/CDFS/merge_data/timing/'
        os.chdir(path)
        fitsname = 'CDFS_3.fits'
        w = WCS(path + fitsname)
        source_info = fits.open('7Ms_catalog.fit')
        ra = source_info[1].data['RAJ2000']
        dec = source_info[1].data['DEJ2000']
        src_x, src_y = w.all_world2pix(ra, dec, 1)

        phy_x = src_x + 2896
        phy_y = src_y + 2896

        src_x = np.rint(src_x)
        src_y = np.rint(src_y)
        src_x = src_x.astype(np.int)
        src_y = src_y.astype(np.int)
        src_x[np.where(src_x > 2399)] = 2399
        src_x[np.where(src_x < 0)] = 0
        src_y[np.where(src_y > 2399)] = 2399
        src_y[np.where(src_y < 0)] = 0

        p90_list = 'reproj_psf90_{0}_b3.fits'.format(obs_id)
        hdul_p90 = fits.open(path + p90_list)
        p90_data = hdul_p90[0].data
        p90_data = p90_data.T
        src_radius = p90_data[src_x, src_y]
        del_index=[]
        for i in range(len(src_radius)):
            r=((x-phy_x[i])**2+(y-phy_y[i])**2)**0.5
            del_index=np.union1d(del_index,np.where(r<2*src_radius[i]))
        del_index=del_index.astype('int64')
        x=np.delete(x,del_index);y=np.delete(y,del_index)
        time=np.delete(time,del_index);energy=np.delete(energy,del_index)
        ID=np.delete(ID,del_index)
        return [x, y, time, energy, ID]


    for item in source_id:
        reg = read_region_bkg(str(item)+'_bkg.reg')
        src_index=where_region_annulus(x,y,reg)
        src_x=x[src_index]
        src_y=y[src_index]
        src_t=time[src_index]
        src_E=energy[src_index]
        src_ID=obs_ID[src_index]


        [src_x,src_y,src_t,src_E,src_ID]=delete_photon_energy(src_x,src_y,src_t,src_E,src_ID)
        [src_x, src_y, src_t, src_E, src_ID] = delete_photon_overlap(src_x, src_y, src_t, src_E, src_ID)

        src_t=src_t.astype('float')
        src_E =src_E.astype('float')
        src_ID=src_ID.astype('int')
        src_txt=np.column_stack((src_t,src_E,src_ID))
        src_txt = src_txt[src_txt[:,0].argsort()]

        np.savetxt(path+'txt_{0}_0.5_8/'.format(obs_id)+str(item)+'_bkg.txt',src_txt,fmt="%.7f  %5.3f  %d")

# for i in range(len(obs_ID_all)):
#     get_txt(obs_ID_all[i])

def merge_txt(src_id):
    def read_region(obs_id,regname):
        reg_file = []
        with open(path+'region_{0}/'.format(obs_id)+'region_90/'+regname, 'r') as file_to_read:
            while True:
                lines = file_to_read.readline() # 整行读取数据
                reg_file.append(lines)
                if not lines:
                    break
                    pass
        region=reg_file[-2][7:-2]
        reg_x,reg_y,reg_r=[float(i) for i in region.split(',')]
        return [reg_x,reg_y,reg_r]

    #epoch_all=np.loadtxt(path+'txt_all_obs/'+'ACIS-I_epoch.txt')
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
        res_temp=np.loadtxt(path+'txt_{0}_0.5_8/'.format(obs_ID_all[i])+str(src_id)+'_bkg.txt')
        if len(res_temp)==0:
            expmap_name = 'reproj_evt2_sou_{0}_i3.fits'.format(obs_ID_all[i])
            expmap = fits.open(path + expmap_name)
            exptime_all = expmap[0].data
            exptime_all = exptime_all.T
            reg = read_region(obs_ID_all[i], str(src_id) + '.reg')
            src_x=int(np.rint(reg[0])-2896);src_y=int(np.rint(reg[1])-2896)
            #print(src_x,src_y)
            #src_x = int(np.rint(reg[0]));src_y = int(np.rint(reg[1]))  #for NSC_I
            if src_x>2399:src_x=2399
            if src_y>2399:src_y=2399
            if src_x<0:src_x=0
            if src_y<0:src_y = 0


            if exptime_all[src_x][src_y]>0:
                #print('keep')
                epoch_ID.append(obs_ID_all[i])
                epoch_start.append(obs_tstart[i])
                epoch_stop.append(obs_tstop[i])
                epoch_expt.append(obs_expt[i])
            else:
                #print('none')
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
    np.savetxt(path+'txt_all_obs_0.5_8/'+'epoch_src_'+str(src_id)+'_bkg.txt',epoch_info,fmt='%15.2f %15.2f %10d %20.2f')
    np.savetxt(path+'txt_all_obs_0.5_8/'+str(src_id)+'_bkg.txt',result,fmt="%.7f  %5.3f  %d")

source_id=np.linspace(1,1055,1055)
source_id=source_id.astype(int)
for item in source_id:
    merge_txt(item)