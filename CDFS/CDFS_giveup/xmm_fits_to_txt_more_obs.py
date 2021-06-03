#!/bin/bash
# -*- coding: utf-8 -*-
# written by Tong
# extract the arrival time of photons in xmm data
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import string
#import correct as correct
from scipy.interpolate import lagrange
from scipy import interpolate
from scipy.fftpack import fft,ifft
import scipy.signal as ss
import random
import pandas as pd
#path='/Volumes/pulsar/xmm_obs_16/'
# path='/Volumes/pulsar/period/'
path="/Volumes/pulsar/xmmCDFS/"

def get_txt(srcName,obsID,mode):
    hdul_evt= fits.open(path+obsID+'/cal/'+mode+'_'+srcName+'_src_filt_time_20sec.fits')
    hdul_evt_bkg=fits.open(path+obsID+'/cal/'+mode+'_'+srcName+'_bkg_filt_time_20sec.fits')
    # # Be careful! Different column if used for other satellite-----
    # x=hdul_evt[1].data['X']
    # y=hdul_evt[1].data['Y']
    energy=hdul_evt[1].data['PI']
    time=hdul_evt[1].data['TIME']

    energy_bkg=hdul_evt_bkg[1].data['PI']
    time_bkg=hdul_evt_bkg[1].data['TIME']

    tstart=hdul_evt[1].header['TSTART']
    tstop=hdul_evt[1].header['TSTOP']

    src_t=time;src_E=energy

    def delete_photon_ID(time,energy):
        i=0
        while i < len(energy):
            if energy[i]>10000 or energy[i]<200:
                energy=np.delete(energy,i)
                time=np.delete(time,i)
                i=i-1
            i=i+1
        return [time,energy]

    # [src_t,src_E]=delete_photon_ID(src_t,src_E)

    src_txt = np.column_stack((src_t,src_E))
    src_txt = src_txt[src_txt[:,0].argsort()]
    bkg_txt = np.column_stack((time_bkg, energy_bkg))
    bkg_txt = bkg_txt[bkg_txt[:, 0].argsort()]
    epoch=np.array([tstart,tstop,int(obsID),tstop-tstart])
    #print src_txt
    os.chdir(path+obsID)
    os.system('mkdir txt')
    #os.system('rm ./txt/*.txt')
    np.savetxt(path +obsID+ '/txt/' +srcName+'_'+mode+ '_20sec.txt', src_txt, fmt="%.7f  %.3f ")
    np.savetxt(path + obsID + '/txt/' + srcName + '_' + mode + '_bkg_20sec.txt', bkg_txt, fmt="%.7f  %.3f ")
    np.savetxt(path + obsID + '/txt/' +'epoch_'+ srcName + '_' + mode + '_20sec.txt', [epoch],fmt="%.2f  %.2f %15d %15.3f")

def merge_txt(srcName):
    for mode in ['mos1','mos2','pn']:
        src_evt = np.loadtxt(path + obsList[0] + '/txt/' + srcName + '_' + mode + '_20sec.txt')
        bkg_evt = np.loadtxt(path + obsList[0] + '/txt/' + srcName + '_' + mode + '_bkg_20sec.txt')
        epoch_evt=np.array([np.loadtxt(path + obsList[0] + '/txt/epoch_' + srcName + '_' + mode + '_20sec.txt')])

        src_obsid = np.zeros(len(src_evt)) + int(obsList[0])
        bkg_obsid = np.zeros(len(bkg_evt)) + int(obsList[0])

        i=1
        while i <len(obsList):
            id=obsList[i]
            src_evt_temp=np.loadtxt(path +id+ '/txt/' +srcName+'_'+mode+ '_20sec.txt')
            bkg_evt_temp=np.loadtxt(path +id+ '/txt/' +srcName+'_'+mode+ '_bkg_20sec.txt')
            epoch_evt_temp=np.array([np.loadtxt(path +id + '/txt/epoch_' + srcName + '_' + mode + '_20sec.txt')])
            src_obsid_temp=np.zeros(len(src_evt_temp))+int(id)
            bkg_obsid_temp=np.zeros(len(bkg_evt_temp)) + int(id)
            if len(src_evt_temp)>3 and len(bkg_evt_temp)>3:
                print(obsList[i])
                src_evt = np.concatenate((src_evt,src_evt_temp))
                bkg_evt = np.concatenate((bkg_evt, bkg_evt_temp))
                src_obsid = np.concatenate((src_obsid,src_obsid_temp))
                bkg_obsid = np.concatenate((bkg_obsid, bkg_obsid_temp))
                epoch_evt=np.concatenate((epoch_evt,epoch_evt_temp))
                i += 1
            else:
                i+=1
        print(epoch_evt)
        src_merge_txt=np.column_stack((src_evt,src_obsid))
        bkg_merge_txt=np.column_stack((bkg_evt,bkg_obsid))
        np.savetxt(path+'txt/{0}_{1}_all_obs_20sec.txt'.format(srcName,mode),src_merge_txt,fmt="%15.7f  %15.3f  %15d")
        np.savetxt(path + 'txt/bkg_{0}_{1}_all_obs_20sec.txt'.format(srcName,mode), bkg_merge_txt, fmt="%15.7f  %15.3f  %15d")
        np.savetxt(path + 'txt/epoch_{0}_{1}_all_obs_20sec.txt'.format(srcName, mode), epoch_evt,fmt="%15.2f  %15.2f %15d %15.3f")

obsList=["0108060401","0108060501","0108060601","0108060701","0108061801","0108061901","0108062101",
         "0108062301","0555780101","0555780201","0555780301","0555780401","0555780501","0555780601",
         "0555780701","0555780801","0555780901","0555781001","0555782301","0604960101","0604960201",
         "0604960301","0604960401","0604961101","0604961201","0604960701","0604960501","0604961301",
         "0604960601","0604960801","0604960901","0604961001","0604961801"]
#
mode=['mos1','mos2','pn']
for det in mode:
    for i in range(len(obsList)):
        get_txt('XID19',obsList[i],det)

merge_txt('XID19')
