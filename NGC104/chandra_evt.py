#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import rocket

path_in='/Volumes/pulsar/NGC6266/merge_data/timing/'
path_out_reg='/Volumes/pulsar/NGC6266/merge_data/timing/region_startover/'
path_out_txt='/Volumes/pulsar/NGC6266/merge_data/timing/txt_startover/'
# obs_ID_all=[78,953,954,955,956,2735,3384, 2736,3385,2737,3386,2738,3387,16527,15747,16529,17420,15748,16528]
obs_ID_all=[ 2677,15761]
# obs_ID_all=[13850,14392,14394,14393,13856,13857,13854,14413,13855,
#             14414,13847,14427,13848,13849,13846,14438,13845,14460,
#             13844,14461,13853,13841,14465,14466,13842,13839,13840,
#             14432,13838,13852,14439,14462,14463,13851,15568,13843,
#             15570,14468]
# ra_center =265.17539;dec_center =-53.67433;inter_radius = 21.61363636363636  ##角秒
# obs_ID_all=[11031]
# ra_center = 154.40343;dec_center = -46.41248;inter_radius = 60  ##角秒

def input_srcinfo():
    cat=fits.open(path_in+'NGC6266_p50_i5_src_1_2_4_8.fits')
    ra = cat[1].data['RA']
    dec = cat[1].data['DEC']

    srcID_list=np.arange(1,len(ra)+1,1)
    return (srcID_list,ra,dec)

def main_process():
    rocket.make_region_each_obs(path_in,path_out_reg,ra=ra,dec=dec,wcsimage=wcsimage,obs_ID_all=obs_ID_all,
                                ecf=90,srcid=srcID_list,multiple_src=1,bkg=1)
    rocket.make_epoch_file(obsid=obs_ID_all,inpath=path_in,outpath=path_out_txt,outname='NGC6266_epoch')
    epoch_info=np.loadtxt(path_out_txt+'NGC6266_epoch.txt')
    for id in obs_ID_all:
        rocket.get_txt(path_in=path_in,path_in_reg=path_out_reg,path_out=path_out_txt,srcid_list=srcID_list,obs_id=id,ecf=90,suffix='_p90')
        rocket.extract_evtlist_bkg(path_in=path_in,path_in_reg=path_out_reg, path_out=path_out_txt, obs_id=id, srcid_list=srcID_list, ecf=90,suffix='_p90')
    for srcid in srcID_list:
        rocket.merge_txt(srcid,epoch_info,inpath=path_out_txt,outpath=path_out_txt,bkg=1,suffix='_p90',outname='txt_all_obs_p90')

if __name__=='__main__':
    (srcID_list, ra, dec) = input_srcinfo()
    # srcID_list=[55293];ra=[5.5125976];dec=[-72.0692052]
    wcsimage = 'NGC6266_5.fits'
    main_process()
    # rocket.select_src_bypos(srcID_list,ra,dec,ra_c=ra_center,dec_c=dec_center,
    #                         inter_radius=inter_radius,outpath=path_out_txt+'txt_all_obs_p50/',
    #                         outname='inter_src_60arcsec.txt')