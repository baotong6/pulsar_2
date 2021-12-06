#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import rocket

path_in='/Volumes/pulsar/47Tuc/merge_data/timing/'
path_out_reg='/Volumes/pulsar/47Tuc/merge_data/timing/region_startover/'
path_out_txt='/Volumes/pulsar/47Tuc/merge_data/timing/txt_startover/'
obs_ID_all=[78,953,954,955,956,2735,3384, 2736,3385,2737,3386,2738,3387,16527,15747,16529,17420,15748,16528]
ra_center = 6.022318;dec_center = -72.081443;inter_radius = 60  ##角秒

def input_srcinfo():
    cat=fits.open(path_in+'xray_properties-592.fits')
    ra = cat[1].data['RAdeg']
    dec = cat[1].data['DEdeg']
    srcID_list=np.arange(1,len(ra)+1,1)
    return (srcID_list,ra,dec)

def main_process():
    # rocket.make_region_each_obs(path_in,path_out_reg,ra,dec,wcsimage=wcsimage,obs_ID_all=obs_ID_all,ecf=90,srcid=None,multiple_src=1,single_name=0)
    # rocket.make_epoch_file(obsid=obs_ID_all,inpath=path_in,outpath=path_out_txt,outname='47Tuc_epoch')
    epoch_info=np.loadtxt(path_in+'47Tuc_epoch.txt')
    # for obsid in obs_ID_all:
    #     rocket.get_txt(path_in, path_out=path_out_txt, reg_name=srcID_list, obs_id=obsid,suffix='_p90')
    for srcid in srcID_list:
        rocket.merge_txt(srcid,epoch_info,inpath=path_out_txt,outpath=path_out_txt,suffix='_p75',outname='txt_all_obs_p75')

if __name__=='__main__':
    (srcID_list, ra, dec) = input_srcinfo()
    wcsimage = '47Tuc_4.fits'
    main_process()
    # rocket.select_src_bypos(srcID_list,ra,dec,ra_c=ra_center,dec_c=dec_center,
    #                         inter_radius=inter_radius,outpath=path_out_txt+'txt_all_obs_p50/',
    #                         outname='inter_src_60arcsec.txt')