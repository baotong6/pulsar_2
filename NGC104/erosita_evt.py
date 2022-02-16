#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import rocket

def main_process():
    # rocket.make_region_each_obs(path_in,path_out_reg,ra,dec,wcsimage=wcsimage,obs_ID_all=obs_ID_all,ecf=90,srcid=None,multiple_src=1,single_name=0)
    # rocket.make_epoch_file(obsid=obs_ID_all,inpath=path_in,outpath=path_out_txt,outname='47Tuc_epoch')
    epoch_info=np.loadtxt(path_in+'47Tuc_epoch.txt')
    for obsid in obs_ID_all:
        rocket.get_txt(path_in, path_in_reg=path_out_reg,path_out=path_out_txt, reg_name=srcID_list, obs_id=obsid,suffix='_p90')
    for srcid in srcID_list:
        rocket.merge_txt(srcid,epoch_info,inpath=path_out_txt,outpath=path_out_txt,suffix='_p90',outname='txt_all_obs_p90')