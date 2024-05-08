#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 15 19:33:00 2023

@author: helin
@author: baotong
"""
import numpy as np
import os
import glob
import subprocess
import csv
import pandas as pd
import math
from astropy.io import fits

list3=['name','R','Rl','Ru','HR','HRl','HRu']
list4=[]

BEHR=os.path.realpath("/Users/baotong/BEHR/BEHR")

# def input_cts(path,specname,ebands,ebandh):
#     specfile=fits.open(path+specname)
#     channel=specfile[1].data['CHANNEL']
#     counts=specfile[1].data['COUNTS']
#     chansl=int(ebands[0]/14.6+1)
#     chansh=int(ebands[1]/14.6+1)
#     chanhl=int(ebandh[0]/14.6+1)
#     chanhh=int(ebandh[1]/14.6+1)
#     softcts=np.sum(counts[chansl:chansh+1])
#     hardcts=np.sum(counts[chanhl:chanhh+1])
#
#     return (softcts,hardcts)

def input_parsfromtxt(path,srcid,ebands,ebandh):
    srcevt=np.loadtxt(path+f'{srcid}.txt')
    bkgevt=np.loadtxt(path+f'{srcid}_bkg.txt')
    if srcevt.ndim== 1:srcevt = np.array([srcevt])

    softsrc=len(np.where((srcevt[:,1]<ebands[1])&(srcevt[:,1]>ebands[0]))[0])
    hardsrc=len(np.where((srcevt[:,1]<ebandh[1])&(srcevt[:,1]>ebandh[0]))[0])
    if len(bkgevt)==0 or bkgevt.ndim==1:
        softbkg=0;hardbkg=0
    # elif bkgevt.ndim==1:
    #     bkgevt = np.array([bkgevt])
    else:
        softbkg=len(np.where((bkgevt[:,1]<ebands[1])&(bkgevt[:,1]>ebands[0]))[0])
        hardbkg=len(np.where((bkgevt[:,1]<ebandh[1])&(bkgevt[:,1]>ebandh[0]))[0])

    return (softsrc,softbkg,hardsrc,hardbkg)
def input_pars(path,specpath,srcid):
    specname=f'{srcid}_stack_src.pi'
    bkgspecname=f'{srcid}_stack_bkg.pi'
    # ##==from spec==##
    # (softsrc, hardsrc) = input_cts(specpath, specname, ebands=[2000, 6000], ebandh=[6000, 7000])
    # (softbkg,hardbkg)=input_cts(specpath,bkgspecname,ebands=[2000,6000],ebandh=[6000,7000])
    ##--from txt--##
    (softsrc, softbkg, hardsrc, hardbkg)=input_parsfromtxt(path+'txt_all_obs_p90/',srcid,ebands=[500, 2000], ebandh=[2000, 8000])

    all_area=np.loadtxt(path+'backscale/all_area.txt')
    srcidlist=all_area[:,0];srcarea_all=all_area[:,1];bkgarea_all=all_area[:,2]
    srcarea=srcarea_all[np.where(srcidlist==srcid)[0]][0]
    bkgarea=bkgarea_all[np.where(srcidlist==srcid)[0]][0]
    if srcarea==0 or bkgarea/srcarea<0:
        print('Warning, area should be positive')
        return None
    else:
        softarea=bkgarea/srcarea;hardarea=softarea
        print(softsrc)
        p1 = "softsrc=" + str(int(softsrc))
        p2 = "hardsrc=" + str(int(hardsrc))
        p3 = "softbkg=" + str(int(softbkg))
        p4 = "hardbkg=" + str(int(hardbkg))
        p5 = "softarea=" + str(softarea)
        p6 = "hardarea=" + str(hardarea)
        # p9="softeff="+softeff
        # p10="hardeff="+hardeff
        os.chdir(path+'HR/S26_H67_evt')

        filename = f'behr_{srcid}'
        p7 = "output=" + filename
        p8 = "outputHR=True"
        p11 = "level=68.0"
        subprocess.call([BEHR, p1, p2, p3, p4, p5, p6, p7, p8, p11])
        str_tmp = "{0:10d} {1:10f} {2:10f} {3:10f} {4:10f} {5:15.5f} {6:15.1f}".format(int(filename[5:]),softsrc,hardsrc,softbkg,hardbkg,softarea,hardarea)
        if not os.path.exists(filename + '.txt'):
            return None
        else:
            with open('all_inputpars.txt', 'a+') as f2:
                f2.writelines(str_tmp + '\n')
            with open( filename + '.txt', "r+") as f:
                data = f.readlines()
                R=data[1].split('\t')[2]
                R_lo=data[1].split('\t')[4]
                R_up=data[1].split('\t')[5]
                HR = data[2].split('\t')[2]
                HR_lo = data[2].split('\t')[4]
                HR_up = data[2].split('\t')[5]
                list4.append([filename, R,R_lo,R_up,HR, HR_lo, HR_up])
                print(filename + ':' + HR + ' ' + HR_lo + ' ' + HR_up)
            test = pd.DataFrame(columns=list3, data=list4)
            test.to_csv(f'{srcid}_HR_sup.csv', index=False)
if __name__=='__main__':
    path='/Users/baotong/Desktop/period_Tuc/'
    for i in range(1,538):
        input_pars(path,specpath=path+'spectra_startover/',srcid=i)
    cmd='cp '+f'{i}_HR_sup.csv '+'all_HR_sup.csv'
    os.system(cmd)

