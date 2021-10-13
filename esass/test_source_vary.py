import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
import pandas as pd
import sys
import os
from tkinter import _flatten
import funcs_fits2txt as funcs
from funcs_fits2txt import Circle
from scipy import interpolate
from astropy.coordinates import SkyCoord
from astropy import units as u

def source_photon():
    ##background 大约为100arcsec 8000个光子
    # bkgscale=4000/(np.pi*100**2)
    ecf=75
    path='/Users/baotong/eSASS/data/raw_data/47_Tuc/'
    srcid=np.arange(1,889,1)
    # srcid=np.loadtxt(path+'txt/inter_srcID.txt')
    # srcid=srcid.astype('int')
    obsIDlist = [700011,700163,700013,700014,700173,700174,700175]
    expT=[25808.86,25268.79,25208.83,25208.82,8809.16,8389.71,8330.68]
    for k in range(len(obsIDlist)):
        obsid=obsIDlist[k]
        radi=[];cts_all=[];bkg_cts_est=[];SNR=[]
        bkgarea=np.loadtxt(path+'txt/txt_psf{0}_{1}/bkg_area.txt'.format(ecf,obsid))
        srcarea=np.loadtxt(path+'txt/txt_psf{0}_{1}/src_area.txt'.format(ecf,obsid))
        backscale=bkgarea[:,1];srcscale=srcarea[:,1]
        for i in range(len(srcid)):
            srctxtfile=path+'txt/txt_psf{0}_{1}/{2}_{1}.txt'.format(ecf,obsid,srcid[i])
            bkgtxtfile=path+'txt/txt_psf{0}_{1}/{2}_bkg_{1}.txt'.format(ecf,obsid,srcid[i])
            if not (os.path.exists(srctxtfile)):
                continue
            time=np.loadtxt(srctxtfile)
            time_bkg=np.loadtxt(bkgtxtfile)
            counts=len(time);bkg_cts=len(time_bkg)
            regfile=path+'reg_{0}/region_{1}/{2}.reg'.format(obsid,ecf,srcid[i])
            [reg_x, reg_y, reg_r]=funcs.read_region(regfile)
            bkg_cts=bkg_cts/backscale[i]*srcscale[i]
            cts_all.append(counts);radi.append(reg_r*3600)
            bkg_cts_est.append(bkg_cts)
            SNR.append((counts-bkg_cts)/np.sqrt(counts))
        info=np.column_stack((srcid,cts_all,bkg_cts_est,radi,SNR))
        np.savetxt(path+'txt/txt_psf{0}_{1}/src_info.txt'.format(ecf,obsid),info,fmt='%10d %10d %10.2f %10.5f %10.5f')

def obs_src_vary():
    bkgscale=4000/(np.pi*100**2)
    ecf=75
    path='/Users/baotong/eSASS/data/raw_data/47_Tuc/'
    srcid=np.loadtxt(path+'txt/inter_srcID.txt')
    srcid=srcid.astype('int')
    obsIDlist = [700011,700163,700013,700014,700173,700174,700175]
    expT=[25808.86,25268.79,25208.83,25208.82,8809.16,8389.71,8330.68]
    SCR=[];NCR=[]

    for k in range(len(obsIDlist)):
        obsid=obsIDlist[k]
        src_info=np.loadtxt(path+'txt/txt_psf{0}_{1}/src_info.txt'.format(ecf,obsid))
        SCR.append(src_info[:,1]/expT[k])
        NCR.append((src_info[:,1]-src_info[:,2])/expT[k])
    SCR=np.array(SCR);NCR=np.array(NCR)
    vary_info=[]
    for i in range(len(srcid)):
        vary_info_temp=np.concatenate(([srcid[i]],SCR[:,i],NCR[:,i],[np.max(SCR[:,i])/np.min(SCR[:,i])]))
        vary_info.append(vary_info_temp)
    vary_info=np.array(vary_info)
    print(vary_info)
    np.savetxt(path+'txt/txt_merge_psf{0}_0.2_5/src_vary_info.txt'.format(ecf),vary_info,fmt='%10d %10.5f %10.5f %10.5f %10.5f '
                                                                                        '%10.5f %10.5f %10.5f %10.5f %10.5f '
                                                                                        '%10.5f %10.5f %10.5f %10.5f %10.5f' 
                                                                                             '%10.5f')

def plot_vary_src():
    bkgscale=4000/(np.pi*100**2)
    ecf=75
    path='/Users/baotong/eSASS/data/raw_data/47_Tuc/'
    srcid=np.loadtxt(path+'txt/inter_srcID.txt')
    srcid=srcid.astype('int')
    obsIDlist = [700011,700163,700013,700014,700173,700174,700175]
    expT=[25808.86,25268.79,25208.83,25208.82,8809.16,8389.71,8330.68]
    np.loadtxt(path+'txt/txt_merge_psf{0}_0.2_5/src_vary_info.txt'.format(ecf))

if __name__=='__main__':
    source_photon()