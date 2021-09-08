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

path_out='/Users/baotong/Desktop/aas/V63/'
path_figure='/Users/baotong/Desktop/CDFS/fig_LC/'


ID_use=[64,153,319,398,457,518,780,872,910,925,950,991]

def get_tex(ID_use):
    os.chdir(path_figure)
    with open(path_out+'plot_figure_star_lc.tex','w+') as f:
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

            f.writelines(r'\includegraphics[angle =0, width = 0.29\textwidth]{%s}'
                         %(path_figure+str(ID_use[i])+'.pdf'))
            f.writelines('\n')
            f.writelines(r'\hfill'+'\n')
            #os.system('cp '+specname+'.ps '+ specname+'.eps')
            # os.system('convert '+specname+'.ps '+ specname+'.pdf')
            # os.system('sips -r 90 {0}.pdf'.format(specname))
            f.writelines(r'\includegraphics[angle =0, width = 0.29\textwidth]{%s}'
                         % (path_figure  + str(ID_use[i]) + '_EP.pdf'))
            f.writelines('\n')
            f.writelines(r'\hfill'+'\n')
            f.writelines(r'\includegraphics[angle =0, width = 0.29\textwidth]{%s}}'
                         % (path_figure + str(ID_use[i]) + '_epoch.pdf'))
            f.writelines('\n')
            j += 3
            if j%12== 0 or i==len(ID_use)-1:
                f.writelines(r'\end{figure*}'+'\n')
                f.writelines(r'\clearpage'+'\n')

            i+=1

        f.writelines(r'\end{document}' + '\n')


    f.close()
#get_tex('NSC',data.ID_NSC)
get_tex(ID_use)
# get_tex('ND', data.ID_ND)