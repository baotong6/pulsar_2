'''
Author: baotong && baotong@smail.nju.edu.cn
Date: 2023-12-20 12:06:46
LastEditors: baotong && baotong@smail.nju.edu.cn
LastEditTime: 2023-12-28 11:17:51
FilePath: /pulsar/swift/extract_GTI_event.py
Description: 

Copyright (c) 2023 by baotong, All Rights Reserved. 
'''
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import astropy.units as u
import astropy.constants as c
import os
import stingray as sr
from stingray.events import EventList
from stingray.lightcurve import Lightcurve
from stingray import Lightcurve, Crossspectrum, sampledata,Powerspectrum,AveragedPowerspectrum
import hawkeye as hawk
path='/Users/baotong/swift/ztf_src1/output/'
id=['00035071001','00035071002','00035071003',
'00035071004','00035071005','00035071006',
'00035071007','00035071008','00035071009','00035071010','00035071011']
def make_gti_evt_and_epoch():
    epoch_all=[];all_evt=np.array([])
    for i in range(len(id)):
        epoch_oneobs=[]
        os.chdir(path+id[i])
        matching_files = [file for file in os.listdir() if 'cl_bary.evt' in file][0]
        file=fits.open(matching_files)
        GTI=file[2].data;GTI=np.array(GTI)

        file2=fits.open(path+id[i]+'/'+'src_30sec.fits')
        time=file2[1].data['TIME']
        energy=file2[1].data['PHA']

        total_array=np.array([])
        for j in range(len(GTI)):
            tstart=GTI[j][0]
            tstop=GTI[j][1]
            tlast=tstop-tstart
            index=np.where((time>tstart)&(time<tstop))[0]
            if tlast>400:
                epoch_oneobs.append([tstart,tstop,id[i],tlast])
                epoch_all.append([tstart,tstop,id[i],tlast])
                evt=np.column_stack((time[index],energy[index]))
                evt2=np.column_stack((time[index],energy[index],[id[i] for k in range(len(time[index]))]))
                total_array=np.vstack([total_array,evt]) if total_array.size else evt
                all_evt=np.vstack([all_evt,evt2]) if all_evt.size else evt2
                
        np.savetxt(path+id[i]+'/'+'epoch_GTI400.txt',epoch_oneobs,fmt='%20s %20s %20s %20s')
        np.savetxt(path+id[i]+'/'+'srcevt_GTI400.txt',total_array,fmt='%20.5f %10d')
    np.savetxt(path+'/'+'epoch_allobs_GTI400.txt',epoch_all,fmt='%20s %20s %20s %20s')
    np.savetxt(path+'/'+'evt_allobs_GTI400.txt',all_evt,fmt='%20s %20s %20s')

def plot_LS():
    path='/Users/baotong/swift/ztf_src1/output/00035071003/'
    filename='srcevt.txt'
    time=np.loadtxt(path+filename)[:,0]
    evt = EventList()
    evt.time=time
    lc_out = evt.to_lc(dt=5)
    freq=np.arange(1/1000.,1/100.,1e-6)
    hawk.get_LS(lc_out.time, lc_out.counts,freq,outpath=None,outname=None,save=False,show=True)
if __name__ == '__main__':
    make_gti_evt_and_epoch()
    # plot_LS()
    



    