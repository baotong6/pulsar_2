get_ipython().run_line_magic("matplotlib", " widget      ")


import warnings
warnings.filterwarnings('ignore')


import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import pandas as pd
import sys
import os
from scipy import interpolate
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from tkinter import _flatten
from astropy.stats import poisson_conf_interval


srcid=np.linspace(1,1055,1055)
srcid=srcid.astype(int)


path='/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8/'


def plot_longT_V(data_file,epoch_file):
    epoch_info = np.loadtxt(epoch_file)
    t_start = epoch_info[:, 0]
    t_end = epoch_info[:, 1]
    obsID = epoch_info[:, 2]
    expT = epoch_info[:, 3]
    time = np.loadtxt(data_file)
    cts=[]
    for i in range(len(obsID)):
        cts.append(len(np.where(time[:,2]==obsID[i])[0]))
    cts=np.array(cts)
    CR=cts/expT
#     print(obsID[np.where(CR>0.0006)])
    if np.min(CR)get_ipython().getoutput("=0:")
        VI=np.max(CR)/np.min(CR)
    else:VI=0
    plt.title(dataname[0:-4]+', VI={0}'.format(VI))
    plt.scatter(t_start,CR,marker='+')
    plt.show()
    


def get_bright_id(cts_limit):
    cts_num=[]
    for i in range(len(srcid)):
        dataname='{0}.txt'.format(srcid[i])
        data_file=path+dataname
        time=np.loadtxt(data_file)[:,0]
        cts_num.append(len(time))
    cts_num=np.array(cts_num)
    return srcid[np.where(cts_num>cts_limit)]


# get_bright_id(1000)


bright_id=get_bright_id(1000)
for id in bright_id:
    dataname='{0}.txt'.format(id)
    data_file=path+dataname
    time=np.loadtxt(data_file)[:,0]
    plot_longT_V(path + dataname, path + 'epoch_src_' + dataname)


obs_ID=[581,1431,441,582,2406,2405,2312,1672,2409,2313,2239,8591,9593,9718,8593,8597,8595,8592,8596,
9575,9578,8594,9596,12043,12123,12044,12128,12045,12129,12135,12046,12047,12137,12138,12055,12213,12048,
12049,12050,12222,12219,12051,12218,12223,12052,12220,12053,12054,12230,12231,12227,12233,12232,12234,16183,
16180,16456,16641,16457,16644,16463,17417,17416,16454,16176,16175,16178,16177,16620,16462,17535,17542,16184,
16182,16181,17546,16186,16187,16188,16450,16190,16189,17556,16179,17573,17633,17634,16453,16451,16461,16191,
16460,16459,17552,16455,16458,17677,18709,18719,16452,18730,16185]

source_id=np.linspace(1,1055,1055)
source_id=source_id.astype(int)


path='/Users/baotong/Desktop/CDFS/'


# def filter_evt_byobstime(epochtime,src_evt,bkg_evt,epoch):
#     useobs=epoch[np.where((epoch[:,0]>epoch1[0])&(epoch[:,0]<epoch1[1]))]
#     print(src_evt[:,-1])
#     print(useobs)
# epoch1=[55e6,100e6]
# epoch2=[300e6,320e6]
# epoch3=[380e6,400e6]
# epoch4=[500e6,580e6]
# k=0;i=122
# src_evt=np.loadtxt(path+'txt_obs_each/txt_{0}_0.5_8/{1}.txt'.format(obs_ID[k],source_id[i]))
# bkg_evt=np.loadtxt(path+'txt_obs_each/txt_{0}_0.5_8/{1}_bkg.txt'.format(obs_ID[k],source_id[i]))
# epoch=np.loadtxt(path+'txt_all_obs_0.5_8/epoch_src_{0}.txt'.format(source_id[i]))
# filter_evt_byobstime(epoch1,src_evt,bkg_evt,epoch)


def get_evt_from_epoch(epochnum,source_id):
    src_evt_all=[];bkg_evt_all=[];epoch_all=[]
    for i in range(len(source_id)):
        epoch=np.loadtxt(path+'txt_all_obs_0.5_8/epoch_src_{0}.txt'.format(source_id[i]))
        useobs=epoch[np.where((epoch[:,0]>epochnum[0])&(epoch[:,0]<epochnum[1]))]
#         print(useobs)
        use_obsid=useobs[:,2].astype(int)
        src_evt=[];bkg_evt=[];
        for k in range(len(use_obsid)):
            src_evt_k=np.loadtxt(path+'txt_obs_each/txt_{0}_0.5_8/{1}.txt'.format(use_obsid[k],source_id[i]))
            bkg_evt_k=np.loadtxt(path+'txt_obs_each/txt_{0}_0.5_8/{1}_bkg.txt'.format(use_obsid[k],source_id[i]))
            if len(src_evt_k)==0:
                continue
            elif type(src_evt_k[0])==type(np.array([1.2])[0]):src_evt_k=[src_evt_k]
                
            if len(src_evt)==0:
                src_evt=src_evt_k
            if len(bkg_evt)==0:
                bkg_evt=bkg_evt_k
               
            if len(bkg_evt_k)==0:
                bkg_evt=bkg_evt
            elif type(bkg_evt_k[0])==type(np.array([1.2])[0]):
                bkg_evt_k=[bkg_evt_k]
                bkg_evt=np.row_stack((bkg_evt,bkg_evt_k))    
            else:
                bkg_evt=np.row_stack((bkg_evt,bkg_evt_k))   
                
            src_evt=np.row_stack((src_evt,src_evt_k))
        src_evt_all.append(src_evt);bkg_evt_all.append(bkg_evt);epoch_all.append(useobs)
    return [src_evt_all,bkg_evt_all,epoch_all]


epoch1=[55e6,100e6]
epoch2=[300e6,320e6]
epoch3=[380e6,400e6]
epoch4=[500e6,580e6]
EPALL=np.array([epoch1,epoch2,epoch3,epoch4])


source_id=srcid
for i in range(len(EPALL)):
    EVTALL=get_evt_from_epoch(EPALL[i],source_id)
    for j in range(len(source_id)):
        if len(EVTALL[0][j])==0:
            np.savetxt(path+'txt_all_obs_0.5_8_ep{0}/{1}.txt'.format(i+1,source_id[j]),
                   EVTALL[0][j])
        else:
            np.savetxt(path+'txt_all_obs_0.5_8_ep{0}/{1}.txt'.format(i+1,source_id[j]),
                   EVTALL[0][j],fmt="get_ipython().run_line_magic(".7f", "  %8.3f  %10d\")")
        if len(EVTALL[1][j])==0:    
            np.savetxt(path+'txt_all_obs_0.5_8_ep{0}/{1}_bkg.txt'.format(i+1,source_id[j]),
                   EVTALL[1][j])
        else:
            np.savetxt(path+'txt_all_obs_0.5_8_ep{0}/{1}_bkg.txt'.format(i+1,source_id[j]),
                   EVTALL[1][j],fmt="get_ipython().run_line_magic(".7f", "  %8.3f  %10d\")")
        np.savetxt(path+'txt_all_obs_0.5_8_ep{0}/epoch_src_{1}.txt'.format(i+1,source_id[j]),
                   EVTALL[2][j],fmt='get_ipython().run_line_magic("15.2f", " %15.2f %10d %20.2f')")


def make_aprates_input():
    os.chdir(path+'txt_all_obs_0.5_8_ep{0}'.format(i+1))
    os.system('mkdir aprates')
    src_cts=len(np.loadtxt('{0}.txt'.format(source_id[j])))
    bkg_cts=len(np.loadtxt('{0}_bkg.txt'.format(source_id[j])))
    epoch=np.loadtxt('epoch_src_{0}.txt'.format(source_id[j]))
    exp_s=1
    aprates_text='aprates n={0} m={1} A_s={2} A_b={3} alpha=0.9 beta=0.02 T_s=1 ' \
                             'E_s={4} eng_s=1 flux_s=1 T_b=1 E_b={5} eng_b=1 flux_b=1 clobber=yes ' \
                             'outfile={6} conf=0.9973'.format(src_cts,bkg_cts,backscale,b_backscale,exp_s,exp_s,str(int(ID))+'_'+str(int(obs))+'_out.par')
    with open('aprates/'+'run_'+str(int(ID))+'_'+str(int(obs))+'.e', 'w+') as f:
        f.writelines(aprates_text)
    


source_id=srcid
def write_4ep_info():
    for j in range(len(source_id)):
        src_cts=[];bkg_cts=[];exptime=[];
        for i in range(len(EPALL)):
            os.chdir(path+'txt_all_obs_0.5_8_ep{0}'.format(i+1))
            src_cts.append(len(np.loadtxt('{0}.txt'.format(source_id[j]))))
            bkg_cts.append(len(np.loadtxt('{0}_bkg.txt'.format(source_id[j]))))
            epoch=np.loadtxt('epoch_src_{0}.txt'.format(source_id[j]))
            if len(epoch)==0:
                exptime.append(0)
            elif type(epoch[0])==type(np.array([1.2])[0]):
                exptime.append(epoch[-1])
            else:
                exptime.append(np.sum(epoch[:,-1]))
        info=np.column_stack((src_cts,bkg_cts,exptime,obstime))
        np.savetxt(path+'ep1-4_info/{0}_4ep_info.txt'.format(source_id[j]),info,fmt='get_ipython().run_line_magic("10d", "  %10d  %10.5f')")



write_4ep_info()


  def plot_var_4ep(src_ID):
    ep4info=np.loadtxt(path+'ep1-4_info/{0}_4ep_info.txt'.format(src_ID))
    expT=ep4info[:,-1]
    src_cts=ep4info[:,0];bkg_cts=ep4info[:,1]
    err=poisson_conf_interval(src_cts, background=bkg_cts/12., confidence_level=0.90,
    interval='kraft-burrows-nousek').T  
    x=[np.mean(EPALL[i]) for i in range(len(EPALL))]
    y=(src_cts-bkg_cts/12.)
    yerr=[0,0]
    yerr[0]=y-err[:,0]
    yerr[1]=err[:,1]-y
    
    y/=expT
    yerr/=expT
    plt.errorbar(x, y, yerr = yerr, fmt = '.', capsize = 1, elinewidth = 1, ecolor = 'red')
    plt.show()


get_ipython().run_line_magic("matplotlib", " widget  ")
plot_var_4ep('716')
