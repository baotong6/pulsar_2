#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import string
#import correct as correct
from datetime import datetime
from scipy.interpolate import lagrange
from scipy import interpolate
from scipy.optimize import curve_fit
import pandas as pd
import read_csv as data

path_out='/Users/baotong/Desktop/aas/pCV_GC/'
path_figure='/Users/baotong/Desktop/aas/pCV_GC/figure/'
###----------------read table---------------###
label='omg_cen'
path_table = '/Users/baotong/Desktop/period_Tuc/'
result_Tuc = pd.read_excel(path_table + 'result_0.5_8_all.xlsx', label)

ID_Tuc=result_Tuc['seq']
P_Tuc=result_Tuc['P_out']
net_percent_Tuc=result_Tuc['net_percent']
###----------------read table---------------###

def get_tex(mode,ID_use):
    os.chdir(path_figure+mode)
    with open(path_out+'plot_figure_{0}.tex'.format(mode),'w+') as f:
        f.writelines(r'\documentclass{aastex63}'+'\n')
        f.writelines(r'\usepackage{graphics}'+'\n')
        f.writelines(r'\addtolength{\oddsidemargin}{-0.75in}'+'\n')
        f.writelines(r'\addtolength{\evensidemargin}{-0.75in}' + '\n')
        f.writelines(r'\addtolength{\textwidth}{1.2in}'+'\n')
        f.writelines(r'\addtolength{\topmargin}{0.1in}' + '\n')
        f.writelines(r'\pagestyle{empty}'+'\n')
        # f.writelines(r'\addtolength{\floatsep}{1.0in}'+'\n')
        f.writelines(r'\begin{document}'+'\n')
        i = 0;j=0
        while i <len(ID_use):
            if j%12== 0:
                f.writelines(r'\begin{figure*}[!ht]'+'\n')
                f.writelines(r'\centering'+'\n')

            f.writelines(r'\includegraphics[angle =0, width = 0.33\textwidth]{%s}'
                         %(path_figure+'{0}/pfold_lc_'.format(mode)+str(ID_use[i])+'.eps'))
            f.writelines('\n')
            f.writelines(r'\hfill'+'\n')
            specname=path_figure+mode+'/'+str(ID_use[i])
            #os.system('cp '+specname+'.ps '+ specname+'.eps')
            # os.system('convert '+specname+'.ps '+ specname+'.pdf')
            # os.system('sips -r 90 {0}.pdf'.format(specname))
            f.writelines(r'\includegraphics[angle =0, width = 0.33\textwidth,bb=-54 126 666 666]{%s}'
                         % (path_figure + mode + '/' + str(ID_use[i]) + '_lc.eps'))
            f.writelines('\n')
            f.writelines(r'\hfill'+'\n')
            f.writelines(r'\rotatebox[origin=c,x=12pt,y=125pt]{270}{\includegraphics[width = 0.23\textwidth]{%s}}'
                         % (path_figure + mode + '/' + str(ID_use[i]) + '.ps'))
            f.writelines('\n')
            j += 3
            if j%12== 0 or i==len(ID_use)-1:
                f.writelines(r'\end{figure*}'+'\n')
                f.writelines(r'\clearpage'+'\n')

            i+=1

        f.writelines(r'\end{document}' + '\n')


    f.close()
#get_tex('NSC',data.ID_NSC)
get_tex(label,ID_Tuc)
# get_tex('ND', data.ID_ND)