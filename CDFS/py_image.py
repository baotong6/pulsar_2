#!/bin/bash
# -*- coding: utf-8 -*-
"""
Created on March 29 23:42:40 2022
@author: baotong
compute multiple events.fits to image in a PDF
Caution: all the directory must be end by '/'
Caution: this code uses ciao tools thus must not be run by IDE, see user_manual.txt for details
"""

import numpy as np
import sys
import os

def make_img_from_ciaotools(obsid,band,inpath,outpath=None,binx=3000,biny=3000,block=1,aimpoint=[4096.5,4096.5],ccd_id=None):

    ## By default, outpath=inpath/imgfile/, unless outpath is defined
    ## e.g. band=[500,8000],in units of eV
    ## binx and biny is for image size, aimpoint locates at the box center
    ## To "block" the image, use e.g. block=2 giving an image each of whose pixels corresponds to 4 pixels in the original image.
    ## by default the aimpoint is for ACIS, if HRC,[16384.5,16384.5] for HRC-I observations, and [32768.5,32768.5] for HRC-S observations

    xgrid=f"{aimpoint[0]-binx/2}:{aimpoint[0]+binx/2}:{block}"
    ygrid=f"{aimpoint[1]-biny/2}:{aimpoint[1]+biny/2}:{block}"
    event_root=f"all_bcc_{obsid}_reproj_evt.fits"
    imgroot=f"{obsid}_img_{band[0]}_{band[1]}_block{block}_{binx}x{biny}"

    os.chdir(inpath)
    if not outpath: outpath=inpath+'imgfile/'
    if not os.path.exists(outpath): os.mkdir(outpath)

    if ccd_id:
        cmd = f"dmcopy \"{event_root}[bin x={xgrid},y={ygrid}][energy={band[0]}:{band[1]}][ccd_id={ccd_id}]\" {outpath}{imgroot}.fits clobber=yes"
    else:
        cmd = f"dmcopy \"{event_root}[bin x={xgrid},y={ygrid}][energy={band[0]}:{band[1]}]\" {outpath}{imgroot}.fits clobber=yes"

    os.system(cmd)

if __name__=='__main__':
    ## e.g. when test in mypath;
    inpath = '/Volumes/pulsar/py_image/'
    band=[500,8000]
    obsID_list = [581, 1431, 441, 582, 2406, 2405, 2312, 1672, 2409, 2313, 2239, 8591, 9593, 9718, 8593, 8597, 8595, 8592, 8596,
                  9575, 9578, 8594, 9596, 12043, 12123, 12044, 12128, 12045, 12129, 12135, 12046, 12047, 12137, 12138, 12055,
                  12213, 12048,12049, 12050, 12222, 12219, 12051, 12218, 12223, 12052, 12220, 12053, 12054, 12230, 12231, 12227,
                  12233, 12232, 12234, 16183,16180, 16456, 16641, 16457, 16644, 16463, 17417, 17416, 16454, 16176, 16175, 16178,
                  16177, 16620,16462, 17535, 17542, 16184,16182, 16181, 17546, 16186, 16187, 16188, 16450, 16190, 16189, 17556,
                  16179, 17573, 17633, 17634, 16453, 16451, 16461, 16191,16460, 16459, 17552, 16455, 16458, 17677, 18709, 18719,
                  16452, 18730, 16185]
    for obsid in obsID_list:
        make_img_from_ciaotools(obsid=obsid,band=band,inpath=inpath,binx=2400,biny=2400)