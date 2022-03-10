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
from astropy.stats import poisson_conf_interval
import scipy
import rocket
import chandra_evt as evt_input

path_in='/Volumes/pulsar/47Tuc/merge_data/timing/'
path_out_reg='/Volumes/pulsar/47Tuc/merge_data/timing/region_startover/'
path_out_txt='/Volumes/pulsar/47Tuc/merge_data/timing/txt_startover/'
obs_ID_all=[78,953,954,955,956,2735,3384, 2736,3385,2737,3386,2738,3387,16527,15747,16529,17420,15748,16528]

ra_center = [6.0223292];
dec_center = [-72.0814444]
inter_radius = 21.61363636363636
def input_srcinfo():
    cat = fits.open(path_in + 'xray_properties-592.fits')
    ra = cat[1].data['RAdeg']
    dec = cat[1].data['DEdeg']
    srcID_list = np.arange(1, len(ra) + 1, 1)
    return (srcID_list, ra, dec)

def select_src_reg(ra_center,dec_center,inter_radius):
    (srcID_list, ra, dec) = input_srcinfo()
    inter_srcID=rocket.select_src_bypos(srcID_list, ra, dec, ra_c=ra_center, dec_c=dec_center,
                            inter_radius=inter_radius,outpath=path_out_txt+'txt_all_obs_p50/',
                            outname='inter_src_60arcsec.txt')
    return inter_srcID

def delete_ID(a,b):
    for i in range(len(b)):
        a=np.delete(a,np.where(a==b[i]))
    return a

def read_regionfile(regname,ifpara=0):
    reg_file=[]
    with open(regname, 'r') as file_to_read:
        while True:
            lines = file_to_read.readline() # 整行读取数据
            reg_file.append(lines)
            if not lines:
                break
                pass
    if ifpara:
        region = reg_file[-2][7:-2]
        reg_x, reg_y, reg_r = [float(i) for i in region.split(',')]
        return [reg_x,reg_y,reg_r]
    else: return reg_file[0]

def make_stack_bkg_reg(path,srcid):
    os.system(f'rm {path}all_select_bkg_smooth.reg')
    core_reg_file = read_regionfile(path + 'CORE.reg')
    with open(path + 'all_select_bkg_smooth.reg', 'a+') as f2:
        f2.writelines(core_reg_file + '\n')
        for i in range(len(srcid)):
            if srcid[i]==224 or 182 or 245:
                print('birght source caution!')
                [reg_x,reg_y,reg_r] = read_regionfile(path + '{0}.reg'.format(srcid[i]),ifpara=1)
                reg_file=f'circle({reg_x},{reg_y},{1.5*reg_r})'
            else:
                reg_file = read_regionfile(path + '{0}.reg'.format(srcid[i]))
            f2.writelines('-'+reg_file + '\n')

def make_stack_reg(path,srcid):
    os.system(f'rm {path}all_select.reg')
    os.system(f'rm {path}all_select_bkg.reg')
    for i in range(len(srcid)):
        reg_file=read_regionfile(path + '{0}.reg'.format(srcid[i]))
        reg_back_file=read_regionfile(path + '{0}_bkg.reg'.format(srcid[i]))
        with open(path+'all_select.reg', 'a+') as f2:
            f2.writelines(reg_file+ '\n')
        with open(path+'all_select_bkg.reg', 'a+') as f2:
            f2.writelines(reg_back_file + '\n')
    return None
if __name__=='__main__':
    inter_srcID=select_src_reg(ra_center,dec_center,inter_radius)
    del_src=[224,245,182,261,211,304,
             294,235,337,254,238,352,236,357,285,163,410,439,194,384,458,540,376,126,366,313,495,16,99]
    inter_srcID_clean=delete_ID(inter_srcID,del_src)
    print(len(inter_srcID_clean),len(inter_srcID))
    for i in range(len(obs_ID_all)):
        path=path_out_reg+f'region_{obs_ID_all[i]}/region_90/'
        make_stack_bkg_reg(path, srcid=inter_srcID)
        make_stack_reg(path,srcid=inter_srcID_clean)
