#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
from astropy.io import fits
import sys
import os
from tkinter import _flatten
from astropy.wcs import WCS
import rocket.make_region as reg_func

def delete_photon_ID(time, energy, ID,e_low,e_high):
    i = 0
    while i < len(energy):
        if energy[i] > e_high or energy[i] < e_low:
            energy = np.delete(energy, i)
            time = np.delete(time, i)
            ID = np.delete(ID, i)
            i = i - 1
        i = i + 1
    return [time, energy, ID]

def make_region_each_obs(path_in,path_out,ra,dec,wcsimage,obs_ID_all,ecf=90,srcid=None,multiple_src=1,single_name=0):
    if srcid==None:
        srcid=np.arange(1,len(ra)+1,1)
    os.chdir(path_in)
    fitsname = wcsimage
    w = WCS(path_in + fitsname)
    src_x, src_y = w.all_world2pix(ra, dec, 1)
    binx,biny=np.shape(fits.open(path_in+fitsname)[0].data)

    out_index=np.union1d(np.where(src_x>binx-1),np.where(src_y>biny-1))
    src_x[out_index]=binx-1;src_y[out_index]=biny-1

    phy_x = src_x + 4096-binx/2
    phy_y = src_y + 4096-biny/2

    src_x = np.rint(src_x)
    src_y = np.rint(src_y)
    src_x = src_x.astype(np.int)
    src_y = src_y.astype(np.int)

    for i in range(len(obs_ID_all)):
        os.chdir(path_out)
        if not os.path.exists('region_{0}'.format(obs_ID_all[i])):
            os.system('mkdir region_{0}'.format(obs_ID_all[i]))
        p90_list='reproj_psf{0}_{1}_b4.fits'.format(ecf,obs_ID_all[i])
        hdul_p90 = fits.open(path_in + p90_list)

        p90_data = hdul_p90[0].data
        p90_data = p90_data.T
        src_radius = p90_data[src_x, src_y]
        src_radius *= 2.032521
        os.chdir(path_out+'region_{0}'.format(obs_ID_all[i]))
        os.system('mkdir region_{0}'.format(ecf))
        if multiple_src:singlename=None
        elif single_name:singlename=single_name

        reg_func.make_phy_reg(srcid=srcid,x=phy_x,y=phy_y,psfradii=src_radius,outpath='./region_{0}/'.format(ecf),singlename=singlename)

    return None
def make_region_merge(path_in,path_out,ra,dec,psfimage,ecf=90,srcid=None,multiple_src=1,single_name=0):
    psfmap=fits.open(path_in+psfimage)

def get_txt(path_in,path_in_reg,path_out,reg_name,obs_id,ecf=90,suffix=''):
    ##对单个obs_id##
    evt_list='all_bcc_{0}_reproj_evt.fits'.format(obs_id)
    hdul_evt= fits.open(path_in+evt_list)
    x=hdul_evt[1].data['x']
    y=hdul_evt[1].data['y']
    energy=hdul_evt[1].data['energy']
    time=hdul_evt[1].data['time']
    obs_ID=np.array([obs_id for i in range(len(time))])

    def where_region(x,y,reg):
        r=np.array((x-reg[0],y-reg[1]))
        len_r=np.sqrt(r[0]**2+r[1]**2)
        temp=len_r-reg[2]
        return np.where(temp<=0)

    os.chdir(path_out)
    if not os.path.exists('txt_{0}{1}'.format(obs_id,suffix)):
        os.system('mkdir txt_{0}{1}'.format(obs_id,suffix))

    for item in reg_name:
        reg = reg_func.read_region(path_in_reg+'region_{0}/'.format(obs_id)+'region_{0}/'.format(suffix[2:])+str(item)+'.reg')
        src_index=where_region(x,y,reg)
        src_x=x[src_index]
        src_y=y[src_index]
        src_t=time[src_index]
        src_E=energy[src_index]
        src_ID=obs_ID[src_index]

        [src_t,src_E,src_ID]=delete_photon_ID(src_t,src_E,src_ID,e_low=500,e_high=8000)
        src_t=src_t.astype('float')
        src_E =src_E.astype('float')
        src_ID=src_ID.astype('int')
        src_txt=np.column_stack((src_t,src_E,src_ID))
        src_txt = src_txt[src_txt[:,0].argsort()]

        np.savetxt(path_out+'txt_{0}{1}/'.format(obs_id,suffix)+str(item)+'.txt',src_txt,fmt="%.7f  %5.3f  %d")

    return None

def make_epoch_file(obsid,inpath,outpath,outname):
    TSTART=[]
    TSTOP=[]
    for id in obsid:
        evtfits=inpath+'all_bcc_{0}_reproj_evt.fits'.format(str(id))
        hdul=fits.open(evtfits)
        TSTART.append(hdul[1].header['TSTART'])
        TSTOP.append(hdul[1].header['TSTOP'])
    TSTART=np.array(TSTART)
    TSTOP=np.array(TSTOP)
    obs_id=np.array(obsid)
    T_exp=TSTOP-TSTART

    epoch=np.column_stack((TSTART,TSTOP,obs_id,T_exp))
    epoch = epoch[epoch[:,0].argsort()]
    np.savetxt(outpath+outname+'.txt',epoch,fmt='%15.2f %15.2f %10d %20.2f')

def merge_txt(src_id,epoch_info,inpath,outpath,suffix,outname='txt_all_obs_0.5_8'):
    ## 对单个点源 ##
    ## Watch out the name of directory  ##
    ## Must follow the style given by get_txt##

    obs_tstart=epoch_info[:,0]
    obs_tstop = epoch_info[:,1]
    obs_ID_all=epoch_info[:,2]
    obs_expt=epoch_info[:,-1]

    obs_ID_all=obs_ID_all.astype(int)

    res_t=[];res_E=[];res_ID=[]
    epoch_ID=[];epoch_start=[];epoch_stop=[];epoch_expt=[]

    for i in range(len(obs_ID_all)):
        res_temp=np.loadtxt(inpath+'txt_{0}{1}/'.format(obs_ID_all[i],suffix)+str(src_id)+'.txt')
        if len(res_temp)==0:continue

        if res_temp.ndim==1:res_temp=np.array(res_temp)
        else:
            epoch_ID.append(obs_ID_all[i])
            epoch_start.append(obs_tstart[i])
            epoch_stop.append(obs_tstop[i])
            epoch_expt.append(obs_expt[i])

            res_t.append(list(res_temp[:, 0]))
            res_E.append(list(res_temp[:, 1]))
            res_ID.append(list((res_temp[:, 2])))

    res_t = list(_flatten(res_t))
    res_E = list(_flatten(res_E))
    res_ID = list(_flatten(res_ID))
    result = np.column_stack((res_t, res_E, res_ID))
    epoch_info = np.column_stack((epoch_start, epoch_stop, epoch_ID, epoch_expt))
    if not os.path.exists(outpath+'{0}/'.format(outname)):
        os.mkdir(outpath+'{0}/'.format(outname))
    np.savetxt(outpath + '{0}/'.format(outname) + 'epoch_src_' + str(src_id) + '.txt', epoch_info,
               fmt='%15.2f %15.2f %10d %20.2f')
    np.savetxt(outpath + '{0}/'.format(outname) + str(src_id) + '.txt', result, fmt="%.7f  %5.3f  %d")