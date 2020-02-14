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
#path='/Volumes/pulsar/GC/merge_data/timing/'
#path='/Volumes/pulsar/SgrA/merge_data/timing/'
path='/Volumes/pulsar/SgrAGRT/merge_data/timing/'
#path='/Volumes/pulsar/LimWin_damage/merge_data/timing/'
# obs_ID=['17239','945','14897', '17236', '17237',
#         '18852','17240','17238','20118',
#         '17241','20807', '20808']
#obs_ID=['945','14897','17236','17239','17237','18852','17240','17238','20118','17241','20807','20808']
# obs_ID=['242' ,'2943' ,'2951' ,'2952' ,'2953' ,'2954' ,'3392' ,'3393' ,'3549'
#       ,'3663' ,'3665' ,'4683' ,'4684' ,'5360' ,'5950' ,'5951' ,'5952' ,'5953' ,'5954' \
#       ,'6113' ,'6363' ,'6639','6640' ,'6641' ,'6642' ,'6643' ,'6644' ,'6645' ,'6646' \
#       ,'7554' ,'7555' ,'7556' ,'7557' ,'7558' ,'7559' ,'9169' ,'9170','9171' ,'9172' \
#       ,'9173' ,'9174' ,'10556' ,'11843' ,'15611' ,'15612' ,'13016' ,'13017' ,'14941','14942']
#obs_ID=['5934','6362','6365','9500','9501','9502','9503','9504','9505','9854','9855','9892','9893']

# obs_ID=['13838','13839','13840','13841','13842','13843','13844','13845','13846','13847' \
#          ,'13848','13849','13850','13851','13852','13853','13854','13855','13856','13857' \
#          ,'14392','14393','14394','14413','14414','14427','14432','14438','14439','14460' \
#          ,'14461','14462','14463','14465','14466','14468','15568','15570']
# epoch_file='ACIS-I_epoch.txt'
#epoch_file='SgrA_I_epoch.txt'
epoch_file='SgrA_G_epoch.txt'
#epoch_file='LW_epoch.txt'
obs_ID_all=np.loadtxt(path+epoch_file)[:,2]
obs_ID_all=obs_ID_all.astype(int)
obs_ID_all=obs_ID_all.astype(str)

def make_region_each_obs():
    os.chdir(path)
    ###for nuclear disk
    # source_info=np.loadtxt('combineobs_info_box.txt')
    # phy_x=source_info[:,2]
    # phy_y=source_info[:,3]
    #
    # phy_x=np.rint(phy_x)
    # phy_y=np.rint(phy_y)
    # phy_x_int=phy_x.astype(np.int)
    # phy_y_int=phy_y.astype(np.int)
    #
    # src_x=phy_x_int-2896
    # src_y=phy_y_int-2896

    ##for NSC
    fitsname = 'SgrA_2000:8000.fits'
    fitsname='SgrAall_2000:8000.fits'
    w = WCS(path + fitsname)

    source_info=fits.open('zhu18_3.fits')
    ra = source_info[1].data.field(0)
    dec = source_info[1].data.field(1)

    #for single source
    ##cnball
    ra=[266.4390417]
    dec=[-28.9747778]

    src_x, src_y = w.all_world2pix(ra, dec, 1)

    src_x = np.rint(src_x)
    src_y = np.rint(src_y)

    src_x = src_x.astype(np.int)
    src_y = src_y.astype(np.int)
    phy_x = src_x + 2896
    phy_y = src_y + 2896

    ##for limiting window
    # fitsname='LW_e1.fits'
    # w = WCS(path + fitsname)
    # source_info = np.loadtxt(path+'catalog_LW.txt')
    # ra = source_info[:,1]
    # dec = source_info[:,2]
    # src_x, src_y = w.all_world2pix(ra, dec, 1)
    #
    # src_x = np.rint(src_x)
    # src_y = np.rint(src_y)
    #
    # src_x = src_x.astype(np.int)
    # src_y = src_y.astype(np.int)
    # phy_x = src_x + 2896
    # phy_y = src_y + 2896



    #将physical坐标转换为img坐标
    for i in range(len(obs_ID_all)):
        os.chdir(path)
        os.system('mkdir region_{0}'.format(obs_ID_all[i]))
        p90_list='reproj_psf90_{0}_e2000:8000.fits'.format(obs_ID_all[i])
        #p90_list = 'reproj_psf90_{0}_b2000:8000.fits'.format(obs_ID_all[i])
        #p90_list = 'reproj_psf90_{0}_b1000:8000.fits'.format(obs_ID_all[i])
        hdul_p90 = fits.open(path + p90_list)

        p90_data = hdul_p90[0].data
        p90_data = p90_data.T
        src_radius = p90_data[src_x, src_y]
        src_radius *= 2.032521
        os.chdir(path+'region_{0}'.format(obs_ID_all[i]))
        os.system('mkdir region_90')

        for i in range(len(phy_x)):
            #with open('./region_90/{0}.reg'.format(i + 1), 'w+') as f1:
            with open('./region_90/{0}.reg'.format('cnball'), 'w+') as f1:
                reg = 'circle(' + str(phy_x[i]) + ',' + str(phy_y[i]) + ',' + str(src_radius[i]) + ')'
                f1.writelines(reg)
            # with open('./region_90/all.reg', 'a+') as f2:
            #     f2.writelines(reg + '\n')
#make_region_each_obs()

def get_txt(obs_id):
    ##for NSC
    # source_id = np.linspace(1, 3619, 3619)
    # source_id = source_id.astype(int)

    ##for LW
    source_id = np.linspace(1, 847, 847)
    source_id = source_id.astype(int)

    ##for single source
    source_id=['cnball']


    os.chdir(path)
    os.system('mkdir txt_{0}'.format(obs_id))
    reg_file=[]
    evt_list='all_bcc_{0}_reproj_evt.fits'.format(obs_id)
    hdul_evt= fits.open(path+evt_list)
    x=hdul_evt[1].data.field(10)
    y=hdul_evt[1].data.field(11)
    energy=hdul_evt[1].data.field(15)
    #be careful for grating obs, energey in column 15
    time=hdul_evt[1].data.field(0)
    obs_ID=np.array([obs_id for i in range(len(time))])

    # x = hdul_evt[1].data.field(0)
    # y = hdul_evt[1].data.field(1)
    # energy = hdul_evt[1].data.field(2)
    # time = hdul_evt[1].data.field(3)
    # obs_ID = hdul_evt[1].data.field(11)

    def read_region(regname):
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

    def where_region(x,y,reg):
        r=np.array((x-reg[0],y-reg[1]))
        len_r=np.sqrt(r[0]**2+r[1]**2)
        temp=len_r-reg[2]
        return np.where(temp<=0)

    def delete_photon_ID(time,energy,ID):
        i=0
        while i < len(energy):
            if energy[i]>8000 or energy[i]<2000:
            #if energy[i] > 8000 or energy[i] < 1000:
                energy=np.delete(energy,i)
                time=np.delete(time,i)
                ID=np.delete(ID,i)
                i=i-1
            i=i+1
        return [time,energy,ID]

    for item in source_id:
        reg = read_region(str(item)+'.reg')
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

        np.savetxt(path+'txt_{0}/'.format(obs_id)+str(item)+'.txt',src_txt,fmt="%.7f  %5.3f  %d")
def merge_txt(src_id):

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
        res_temp=np.loadtxt(path+'txt_{0}/'.format(obs_ID_all[i])+str(src_id)+'.txt')
        if len(res_temp)==0:
            continue
        elif type(res_temp[0])==type(np.array([1.2])[0]):
            ##判断是否只有1个光子，数据类型的bug
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
    np.savetxt(path+'txt_all_obs/'+'epoch_src_'+str(src_id)+'.txt',epoch_info,fmt='%15.2f %15.2f %10d %20.2f')
    np.savetxt(path+'txt_all_obs/'+str(src_id)+'.txt',result,fmt="%.7f  %5.3f  %d")


def merge_I_with_GRT():
    path1= '/Users/baotong/Desktop/period/txt_all_obs_I/'
    path2= '/Users/baotong/Desktop/period/txt_all_obs_G/'
    path_out='/Users/baotong/Desktop/period/txt_all_obs_IG/'
    src_listfile='G_src.txt'
    src_list=np.loadtxt(path2+src_listfile)
    src_list=src_list.astype(int)
    for i in range(len(src_list)):
        ACIS_I_list=np.loadtxt(path1+str(src_list[i])+'.txt')
        grt_list=np.loadtxt(path2+str(src_list[i])+'.txt')
        epoch_I_list=np.loadtxt(path1+'epoch_src_'+str(src_list[i])+'.txt')
        epoch_G_list = np.loadtxt(path2 + 'epoch_src_' + str(src_list[i]) + '.txt')
        # print(ACIS_I_list)
        # print(grt_list)
        combine_list=np.concatenate((ACIS_I_list,grt_list))
        combine_list=combine_list[np.lexsort(combine_list[:, ::-1].T)]

        combine_list_epoch =np.concatenate((epoch_I_list,epoch_G_list))
        combine_list_epoch =combine_list_epoch[np.lexsort(combine_list_epoch[:, ::-1].T)]

        np.savetxt(path_out+str(src_list[i])+'.txt',combine_list,fmt="%20.7f  %10.3f  %10d")
        np.savetxt(path_out+'epoch_src_' + str(src_list[i]) + '.txt',combine_list_epoch,fmt="%20.7f  %20.7f %10d %15.2f ")
#merge_I_with_GRT()

#merge_result_cluster()

#make_region_each_obs()

# ##for ND
# source_id=np.linspace(1,518,518)
# source_id=source_id.astype(int)
#
# ##for NSC
# source_id=np.linspace(1,3619,3619)
# source_id=source_id.astype(int)

##for LW
source_id=np.linspace(1,847,847)
source_id=source_id.astype(int)
# #

for i in range(len(obs_ID_all)):
    get_txt(obs_ID_all[i])

merge_txt('cnball')
# for item in source_id:
#     merge_txt(item)

