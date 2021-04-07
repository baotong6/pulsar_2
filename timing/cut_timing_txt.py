#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
#path='/Users/baotong/Desktop/period_LW/txt_all_obs'

### aprates文件夹中的文件要在plot_aprates_xxx.py中来做 ###
path='/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8'
def exclude(src_ID,obs_ID,mode):
    obs_ID=np.array(obs_ID)
    os.chdir(path)
    event=np.loadtxt('{0}.txt'.format(src_ID))
    epoch=np.loadtxt('epoch_src_{0}.txt'.format(src_ID))
    time=event[:,0]
    index_temp=np.array([])
    index_id=[]
    for i in range(len(obs_ID)):
        obs=obs_ID[i]
        tstart=epoch[:,0][np.where(epoch[:,2]==obs)]
        tstop =epoch[:,1][np.where(epoch[:,2]==obs)]
        index_id.append(np.where(epoch[:,2]==obs)[0][0])
        if mode=='exclude':
            index=np.where((time<tstart)|(time>tstop))[0]
        if mode=='keep':
            index=np.where((time > tstart) & (time < tstop))[0]
        if i==0:
            index_temp =index
        elif mode=='exclude':
            index_temp=np.intersect1d(index_temp,index)
        elif mode=='keep':
            index_temp=np.union1d(index_temp,index)

    event = event[index_temp]

    if mode=='exclude':
        index_id_exclude=np.setdiff1d(np.linspace(0,len(epoch)-1,len(epoch)),np.array(index_id))
        index_id_exclude=index_id_exclude.astype('int')
        #print(index_id_exclude)
        epoch=epoch[index_id_exclude]
    if mode == 'keep':
        epoch=epoch[index_id]

    np.savetxt('{0}_{1}.txt'.format(src_ID,mode),event,fmt="%20.7f  %10.3f  %10d")
    np.savetxt('epoch_src_{0}_{1}.txt'.format(src_ID, mode), epoch, fmt = "%20.7f  %10.3f  %10d %10.2f")

#exclude('1677',[6641,6363,6643,6646,7554,7555,7556,7557,7558,7559],'exclude')
exclude('872',[16453],'keep'),
# obsid_1671=[  242. ,2952.,  2953.  ,6640. , 6641.,  7556., 13850., 14392. ,14394., 14393.,
#  13853., 13841., 14465., 14466., 13842. ,13839., 13840., 14432. ,13838. ,13852.,
#  14439.]
#exclude('202',obsid_1671,'keep')