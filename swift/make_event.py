'''
Author: baotong && baotong@smail.nju.edu.cn
Date: 2023-12-15 15:28:55
LastEditors: baotong && baotong@smail.nju.edu.cn
LastEditTime: 2023-12-20 12:01:56
FilePath: /pulsar/swift/make_event.py
Description: 

Copyright (c) 2023 by baotong, All Rights Reserved. 
'''
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import astropy.units as u
import astropy.constants as c
import os

path='/Users/baotong/swift/ztf_src1/output/'
id=['00035071001','00035071002','00035071003',
'00035071004','00035071005','00035071006',
'00035071007','00035071008','00035071009','00035071010','00035071011']
def extract_evt_singleobs():
    for i in range(len(id)):
        file=fits.open(path+id[i]+'/'+'src_30sec.fits')
        time=file[1].data['TIME']
        energy=file[1].data['PHA']
        evt=np.column_stack((time,energy))
        np.savetxt(path+id[i]+'/'+'srcevt.txt',evt,fmt='%20.5f %10d')
    return None

def extract_evt_allobs():
    total_array = np.array([])
    for i in range(len(id)):
        file=fits.open(path+id[i]+'/'+'src_30sec.fits')
        time=file[1].data['TIME']
        energy=file[1].data['PHA']
        obsid=np.array([id[i] for j in range(len(time))])
        evt=np.column_stack((time,energy,obsid))
        total_array=np.vstack([total_array,evt]) if total_array.size else evt
    print(total_array[0][0])
    np.savetxt(path+'srcevt_allobs.txt',total_array,fmt='%20s %20s %20s')
    # print(total_array)

def make_epoch_allobs():
    epoch=[]
    for i in range(len(id)):
        os.chdir(path+id[i])
        matching_files = [file for file in os.listdir() if 'cl_bary.evt' in file][0]
        file=fits.open(matching_files)
        GTI=file[2].data
        GTI=np.array(GTI)
        for j in range(len(GTI)):
            tstart=GTI[j][0]
            tstop=GTI[j][1]
            tlast=tstop-tstart
            epoch.append([tstart,tstop,id[i],tlast])
    epoch=np.array(epoch)
    print(epoch)
    np.savetxt(path+'epoch_allobs.txt',epoch,fmt='%20s %20s %20s %20s')
    # for i in range(len(id)):
    #     file=fits.open(path+id[i]+'/'+'src_30sec.fits')
    #     time=file[1].data['TIME']
    #     energy=file[1].data['PHA']
def plot_GTI():
    lcfile=fits.open(path+'/00035071011/sw00035071011xpcw2po_cl_bary.evt')
    GTI=lcfile[2].data
    GTI=np.array(GTI)
    GTI_array = np.array([list(x) for x in GTI.tolist()])
    t0=GTI_array[:,0];t1=GTI_array[:,1]
    evt=np.loadtxt(path+'/00035071011/'+'srcevt.txt')
    print(len(t0))
    for i in range(len(t0)):
        plt.plot([t0[i],t1[i]],[1,1],'-')
        plt.scatter(evt[:,0],np.zeros(len(evt))+1,s=5,color='black')
    plt.show()
if __name__ == '__main__':
    # make_epoch_allobs()
    plot_GTI()
