#!/bin/bash
# -*- coding: utf-8 -*-
"""
Created on March 29 23:42:40 2022
@author: baotong
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib.colors import LogNorm
from astropy.wcs import WCS
from astropy.visualization import astropy_mpl_style
from tkinter import _flatten
from itertools import chain
import PyPDF2
import sys,os

# plt.style.use(astropy_mpl_style)
font1 = {'family': 'Normal',
         'weight': 'normal',
         'size': 16, }

path='/Volumes/pulsar/py_image/'

obsID_list = [581, 1431, 441, 582, 2406, 2405, 2312, 1672, 2409, 2313, 2239, 8591, 9593, 9718, 8593, 8597, 8595, 8592,
              8596,
              9575, 9578, 8594, 9596, 12043, 12123, 12044, 12128, 12045, 12129, 12135, 12046, 12047, 12137, 12138,
              12055,
              12213, 12048, 12049, 12050, 12222, 12219, 12051, 12218, 12223, 12052, 12220, 12053, 12054, 12230, 12231,
              12227,
              12233, 12232, 12234, 16183, 16180, 16456, 16641, 16457, 16644, 16463, 17417, 17416, 16454, 16176, 16175,
              16178,
              16177, 16620, 16462, 17535, 17542, 16184, 16182, 16181, 17546, 16186, 16187, 16188, 16450, 16190, 16189,
              17556,
              16179, 17573, 17633, 17634, 16453, 16451, 16461, 16191, 16460, 16459, 17552, 16455, 16458, 17677, 18709,
              18719,
              16452, 18730,16185]

def labelmaker(ax,wcs,ra,dec,binx=2400,biny=2400):
    ## By default, label region is circle
    ## By default binx and biny is 2400, caution if you want to modify those
    srcid=np.arange(1,len(ra)+1,1)
    src_x, src_y = wcs.all_world2pix(ra, dec, 1)
    for i in range(len(src_x)):
        if 0<src_x[i]<binx and 0<src_y[i]<biny:
            ax.scatter(src_x[i],src_y[i],marker='o',edgecolor='white',c=" ",s=200)
            ax.text(src_x[i]-50,src_y[i]+100,f"No.{srcid[i]}",color='white',fontsize='large')

def plot_onepage_img(inpath,imgroot_list,outname,outpath,grid=None,save=1,show=0,label=0,ra=None,dec=None):
    ## By default, out pdf is 3x2 subplots
    ## if label is True, ra and dec must be given
    figlabel = [321,322,323,324,325,326]
    fig=plt.figure(1,(10,14.143))
    i=0
    while i < len(imgroot_list):
        imgroot=imgroot_list[i]
        with fits.open(inpath + f'{imgroot}.fits') as hdu:
            data = hdu[0].data
            wcs = WCS(hdu[0].header)

        ax_temp=plt.subplot(figlabel[i],projection=wcs)
        if label:labelmaker(ax_temp,wcs,ra,dec)
        ax_temp.set_ylabel('DEC')
        ax_temp.set_xlabel('RA')
        ax_temp.set_title(f"{imgroot}",font1)
        map=ax_temp.imshow(data,cmap='gist_heat',vmin=0.01,vmax=0.15)
        i+=1
    plt.subplots_adjust(bottom=0.05, left=0.05,right=0.95, top=0.95)
    if grid:plt.grid(color='white', ls=':', alpha=0.7)
    if save:plt.savefig(outpath+f"{outname}.pdf")
    if show:plt.show()
    plt.close()

def make_all_img(obsID_list,suffix,path_parent,outname,label=0,ra=None,dec=None):
    pagelist=[]
    if not os.path.exists(path_parent+'pdf_file/'): os.mkdir(path_parent+'pdf_file/')
    imgname_list=[f"{obsid}{suffix}" for obsid in obsID_list]
    i=0
    while i < len(imgname_list):
        plot_onepage_img(inpath=path_parent+'imgfile/',outpath=path_parent+'pdf_file/',imgroot_list=imgname_list[i:i+6],
                         outname=f"{outname}_{int((i+6)/6)}",save=0,show=1,label=label,ra=ra,dec=dec)
        pagelist.append(f"{outname}_{int((i+6)/6)}.pdf")
        i+=6
    if i-len(imgname_list)>0:
        plot_onepage_img(inpath=path_parent+'imgfile/', outpath=path_parent+'pdf_file/',imgroot_list=imgname_list[i-6:],
                         outname=f"{outname}_{int(i/6)}",save=0,show=1,label=label,ra=ra,dec=dec)

    return pagelist

def merge_pdf(inpath,filenames,outname):
    os.chdir(inpath)
    merger = PyPDF2.PdfFileMerger()
    for filename in filenames:
        merger.append(PyPDF2.PdfFileReader(filename))
    os.system(f"rm {outname}.pdf")
    merger.write(f"{outname}.pdf")

if __name__=='__main__':
    band=[500,8000];block=1;binx=biny=2400
    suffix=f"_img_{band[0]}_{band[1]}_block{block}_{binx}x{biny}"
    # coord=np.loadtxt(path+'radec.txt')
    # ra=coord[:,0];dec=coord[:,1]

    filenames=make_all_img(obsID_list[0:6],suffix,path_parent=path,outname=f"page_{band[0]}_{band[1]}",label=1,ra=[53,53],dec=[-27,-27.5])
    # merge_pdf(inpath=path+'pdf_file/', filenames=filenames,outname=f"allobs_{band[0]}_{band[1]}")