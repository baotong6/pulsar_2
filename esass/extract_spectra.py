# import numpy as np
# import matplotlib.pyplot as plt
# from astropy.io import fits
# import pandas as pd
import sys
import os

def read_region(regname):
    reg_file=[]
    with open(regname, 'r') as file_to_read:
        while True:
            lines = file_to_read.readline()  # 整行读取数据
            reg_file.append(lines)
            if not lines:
                break
                pass
    region = reg_file[-2][7:-2]
    reg_x, reg_y, reg_r = [float(i) for i in region.split(',')]
    return [reg_x, reg_y, reg_r]

spectra=1

obsIDlist=[700011,700013,700014,700163,700173,700174,700175]
srclist=['715']

if spectra:
     modelist = ['50', '75', 'SNRpsf']  ##保持与前面做region时同样的路径
     for obsid in obsIDlist:
          for srcname in srclist:
               mode=modelist[0]
               # path='/Users/baotong/eSASS/data/raw_data/47_Tuc/'
               # os.chdir(path)
               path_reg='./reg_{0}/'.format(obsid)+'region_{0}/'.format(mode)
               [ra, dec, src_radius]=read_region(path_reg+srcname+'.reg')
               src_radius*=3600
               bkg_radius1=2*src_radius;bkg_radius2=4*src_radius
               evtname='pm00_{0}_020_EventList_c001.fits'.format(obsid)

               cmd= 'srctool eventfiles="{0}" '.format(evtname)+\
                    'srccoord="fk5;{0},{1}" '.format(ra,dec)+\
                    'prefix="spectra/{0}_psf{1}_" '.format(srcname,mode)+\
                    'suffix="_{0}" '.format(obsid)+\
                    'todo="SPEC ARF RMF EVENTS NOSRCGTI" '+\
                    'insts="1 2 3 4 5 6 7" '+\
                    "srcreg='icrs;circle * * {0}\"' ".format(src_radius)+\
                    "backreg='icrs;annulus *,*,{0}\" {1}\"' ".format(bkg_radius1,bkg_radius2)+\
                    'exttype="POINT" '+\
                    'tstep=0.05 '+\
                    'xgrid="0.5 1.0" '+\
                    'gtitype="GTI" '+\
                    'psftype="2D_PSF" '+\
                    'flagsel=255 '+\
                    'clobber="yes"'

               print("cmd="+str(cmd))
               os.system(cmd)
