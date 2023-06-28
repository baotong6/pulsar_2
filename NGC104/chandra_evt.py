#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import rocket

path_in='/Volumes/pulsar/terzan5/merge_data/timing/'
path_out_reg='/Volumes/pulsar/terzan5/merge_data/timing/region_startover/'
path_out_txt='/Volumes/pulsar/terzan5/merge_data/timing/txt_startover/'
obs_ID_all=[3798,10059,13225,13252,13705,14339,13706,
            14475,14476,14477,14625,15615,14478,
            14479,16638,15750,17779,18881]
# obs_ID_all=[ 2677,15761]
# obs_ID_all=[13850,14392,14394,14393,13856,13857,13854,14413,13855,
#             14414,13847,14427,13848,13849,13846,14438,13845,14460,
#             13844,14461,13853,13841,14465,14466,13842,13839,13840,
#             14432,13838,13852,14439,14462,14463,13851,15568,13843,
#             15570,14468]
ra_center = 6.0236250 ;dec_center =-72.0812833;inter_radius =3.17*60  ##角秒
# obs_ID_all=[11031]
# ra_center = 154.40343;dec_center = -46.41248;inter_radius = 60  ##角秒

def input_srcinfo():
    cat=fits.open(path_in+'cheng2019_terzan.fit')
    ra = cat[1].data['RAJ2000']
    dec = cat[1].data['DEJ2000']

    srcID_list=np.arange(1,len(ra)+1,1)
    return (srcID_list,ra,dec)

def main_process():
    # rocket.make_region_each_obs(path_in,path_out_reg,ra=ra,dec=dec,wcsimage=wcsimage,obs_ID_all=obs_ID_all,
    #                             ecf=50,srcid=srcID_list,multiple_src=1,bkg=1)
    # rocket.make_region_each_obs(path_in,path_out_reg,ra=[267.0195083],dec=[-24.7809444],wcsimage=wcsimage,obs_ID_all=obs_ID_all,
    #                             ecf=90,srcid=['TerO'],multiple_src=0,bkg=1,single_name='TerO',single_srcradius=0.7)
    # rocket.make_epoch_file(obsid=obs_ID_all,inpath=path_in,outpath=path_out_txt,outname='NGC6304_epoch')
    epoch_info=np.loadtxt(path_out_txt+'terzan5_epoch.txt')
    ecf=90
    for id in obs_ID_all:
    #     rocket.get_txt(path_in=path_in,path_in_reg=path_out_reg,path_out=path_out_txt,srcid_list=srcID_list,obs_id=id,ecf=ecf,suffix=f'_p{ecf}')
        rocket.extract_evtlist_bkg(path_in=path_in,path_in_reg=path_out_reg, path_out=path_out_txt, obs_id=id,
                                   srcid_list=srcID_list, ecf=ecf,suffix=f'_p{ecf}')
    for srcid in srcID_list:
        rocket.merge_txt(srcid,epoch_info,inpath=path_out_txt,outpath=path_out_txt,bkg=1,suffix=f'_p{ecf}',outname=f'txt_all_obs_p{ecf}')

if __name__=='__main__':
    # (srcID_list, ra, dec) = input_srcinfo()
    srcID_list=['TerO']
    # srcID_list=[55293];ra=[5.5125976];dec=[-72.0692052]
    wcsimage = 'terzan5_5.fits'
    main_process()
    # rocket.select_src_bypos(srcID_list,ra,dec,ra_c=ra_center,dec_c=dec_center,
    #                         inter_radius=1*60,outpath=path_out_txt+'txt_all_obs_p90/',
    #                         outname='inter_src_1arcm.txt',save=1)