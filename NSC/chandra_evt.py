'''
Author: baotong && baotong@smail.nju.edu.cn
Date: 2021-12-05 22:29:53
LastEditors: baotong && baotong@smail.nju.edu.cn
LastEditTime: 2024-04-02 09:47:39
FilePath: /pulsar/NSC/chandra_evt.py
Description: 

Copyright (c) 2024 by baotong, All Rights Reserved. 
'''
#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import rocket
import warnings
# 禁止所有UserWarning
warnings.filterwarnings("ignore", category=UserWarning)


path_in='/Volumes/pulsar/SgrA/merge_data/timing/'
path_out_reg='/Volumes/pulsar/SgrA/merge_data/timing/region_startover/'
path_out_txt='/Volumes/pulsar/SgrA/merge_data/timing/txt_startover/'
obs_ID_all=[242,15611,15612, 2951, 2952, 2953, 2954, 2943, 3663, 3392, 
            3393, 3665, 3549, 4683, 4684, 5360, 6113, 5950, 5951, 5952, 
            5953, 5954, 6639, 6640, 6641, 6642, 6363, 6643, 6644, 6645, 
            6646, 7554, 7555, 7556, 7557, 7558, 7559, 9169, 9170, 9171, 
            9172, 9174, 9173,10556,11843,13016,13017,14941,14942]
# ra_center = 6.0236250 ;dec_center =-72.0812833;inter_radius =3.17*60  ##角秒

def input_srcinfo():
    cat=fits.open(path_in+'zhu18_3.fits')
    ra = cat[1].data['_RAJ2000']
    dec = cat[1].data['_DEJ2000']
    srcID_list=np.arange(1,len(ra)+1,1)
    return (srcID_list,ra,dec)
def main_process():
    # rocket.make_region_each_obs(path_in,path_out_reg,ra=ra,dec=dec,wcsimage=wcsimage,obs_ID_all=obs_ID_all,
    #                             ecf=90,srcid=srcID_list,multiple_src=1,bkg=1)
    # rocket.make_region_each_obs(path_in,path_out_reg,ra=[267.0195083],dec=[-24.7809444],wcsimage=wcsimage,obs_ID_all=obs_ID_all,
    #                             ecf=90,srcid=['TerO'],multiple_src=0,bkg=1,single_name='TerO',single_srcradius=0.7)
    # rocket.make_epoch_file(obsid=obs_ID_all,inpath=path_in,outpath=path_out_txt,outname='NGC6304_epoch')
    epoch_info=np.loadtxt(path_out_txt+'SgrA_I_epoch.txt')
    ecf=90
    # for id in obs_ID_all[33:]:
    #     rocket.get_txt(path_in=path_in,path_in_reg=path_out_reg,path_out=path_out_txt,srcid_list=srcID_list,obs_id=id,ecf=ecf,suffix=f'_p{ecf}')
    #     rocket.extract_evtlist_bkg(path_in=path_in,path_in_reg=path_out_reg, path_out=path_out_txt, obs_id=id,
    #                                srcid_list=srcID_list, ecf=ecf,suffix=f'_p{ecf}')
    for srcid in srcID_list[384:]:
        rocket.merge_txt(srcid,epoch_info,inpath=path_out_txt,outpath=path_out_txt,bkg=1,suffix=f'_p{ecf}',outname=f'txt_all_obs_p{ecf}')


def merge_I_with_GRT():
    path1= '/Users/baotong/Desktop/period/txt_startover_I/txt_all_obs_p90/'
    path2= '/Users/baotong/Desktop/period/txt_startover_G/txt_all_obs_p90/'
    path_out='/Users/baotong/Desktop/period/txt_startover_IG/txt_all_obs_p90/'
    # src_listfile='G_src.txt'
    # src_list=np.loadtxt(path2+src_listfile)
    src_list=np.linspace(1,3619,3619)
    src_list=src_list.astype(int)
    for i in range(len(src_list)):
        ACIS_I_list=np.loadtxt(path1+str(src_list[i])+'.txt')
        grt_list=np.loadtxt(path2+str(src_list[i])+'.txt')
        epoch_I_list=np.loadtxt(path1+'epoch_src_'+str(src_list[i])+'.txt')
        epoch_G_list = np.loadtxt(path2 + 'epoch_src_' + str(src_list[i]) + '.txt')
        # print(ACIS_I_list)
        #print(len(grt_list))
        if len(epoch_G_list)>0:
            # src_list_both.append(src_list[i])
            if len(grt_list)==3 and type(grt_list[0])==type(np.array([1.2])[0]):
                grt_list=[grt_list]
            if len(epoch_G_list)==4 and type(epoch_G_list[0])==type(np.array([1.2])[0]) :
                epoch_G_list=[epoch_G_list]

            combine_list=np.concatenate((ACIS_I_list,grt_list))
            combine_list=combine_list[np.lexsort(combine_list[:, ::-1].T)]

            combine_list_epoch =np.concatenate((epoch_I_list,epoch_G_list))
            combine_list_epoch =combine_list_epoch[np.lexsort(combine_list_epoch[:, ::-1].T)]

            np.savetxt(path_out+str(src_list[i])+'.txt',combine_list,fmt="%20.7f  %10.3f  %10d")
            np.savetxt(path_out+'epoch_src_' + str(src_list[i]) + '.txt',combine_list_epoch,fmt="%20.7f  %20.7f %10d %15.2f ")
        elif len(epoch_G_list)==0:
            combine_list=ACIS_I_list
            combine_list_epoch=epoch_I_list
            np.savetxt(path_out + str(src_list[i]) + '.txt', combine_list, fmt="%20.7f  %10.3f  %10d")
            np.savetxt(path_out + 'epoch_src_' + str(src_list[i]) + '.txt', combine_list_epoch,
                       fmt="%20.7f  %20.7f %10d %15.2f ")



if __name__=='__main__':
    (srcID_list, ra, dec) = input_srcinfo()
    # srcID_list=['TerO']
    # srcID_list=[55293];ra=[5.5125976];dec=[-72.0692052]
    wcsimage = 'SgrA_5.fits'
    # main_process()
    merge_I_with_GRT()
    # rocket.select_src_bypos(srcID_list,ra,dec,ra_c=ra_center,dec_c=dec_center,
    #                         inter_radius=1*60,outpath=path_out_txt+'txt_all_obs_p90/',
    #                         outname='inter_src_1arcm.txt',save=1)