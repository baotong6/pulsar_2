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

DN_name=['HT_CAS','OY_CAR','QZ_VIR','RU_PEG','SS_AUR',
      'V893_SCO','VW_HYI','YZ_CNC','V405_PEG','WX_HYI']
IP_name=['AO_PSC','DW_CNC','FO_AQR','HT_CAM','J1509_6649',
         'J1649_3307','J1719_4100','J1817_2508','J1830_1232','XY_ARI']

period_DN=[6363.1008,5453.65,5218.56,32365.44,15793.91,
        6563.0304,6417.0144,7499.52,15348.7008,6704.64]
period_IP=[14325.12,5166.1152,17457.984,5159.1168,21202.56,
           13020.48,14420.16,5514.048,19344.96,21833.0208]
spin_IP=[805.2,2315.026,1254.284,514.6,809.42,
         571.9,1062,1660.8,1820,206.298]
path_figure='/Users/baotong/Desktop/aas/V63/figure/'
path_out='/Users/baotong/Desktop/aas/V63/'
mode='CV'
with open(path_out + 'plot_figure_{0}.tex'.format(mode), 'w+') as f:
    f.writelines(r'\documentclass{aastex63}' + '\n')
    f.writelines(r'\usepackage{graphics}' + '\n')
    f.writelines(r'\addtolength{\oddsidemargin}{-0.75in}' + '\n')
    f.writelines(r'\addtolength{\evensidemargin}{-0.75in}' + '\n')
    f.writelines(r'\addtolength{\textwidth}{1.5in}' + '\n')
    f.writelines(r'\addtolength{\topmargin}{0.6in}' + '\n')
    # f.writelines(r'\addtolength{\floatsep}{1.0in}'+'\n')
    f.writelines(r'\begin{document}' + '\n')
    i = 0;
    j = 0
    while i < len(DN_name):
        if j % 12 == 0:
            f.writelines(r'\begin{figure*}[!ht]' + '\n')
            f.writelines(r'\centering' + '\n')

        f.writelines(r'\includegraphics[angle =0, width = 0.32\textwidth]{%s}'
                     % (path_figure + '{0}/pfold_lc_'.format(mode) + str(DN_name[i]) + '.eps'))
        f.writelines('\n')
        f.writelines(r'\hfill' + '\n')
        # os.system('cp '+specname+'.ps '+ specname+'.eps')
        # os.system('convert '+specname+'.ps '+ specname+'.pdf')
        # os.system('sips -r 90 {0}.pdf'.fformat(specname))
        f.writelines(r'\includegraphics[angle =0, width = 0.32\textwidth]{%s}'
                    % (path_figure + '{0}/pfold_lc_'.format(mode) + str(IP_name[i]) + '.eps'))
        f.writelines('\n')
        f.writelines(r'\hfill' + '\n')
        f.writelines(r'\includegraphics[angle =0, width = 0.32\textwidth]{%s}'
                    % (path_figure + '{0}/pfold_lc_'.format(mode) + str(IP_name[i]) + '_spin.eps'))
        f.writelines('\n')
        j += 3
        if j % 12 == 0 or i == len(DN_name) - 1:
            f.writelines(r'\end{figure*}' + '\n')
            f.writelines(r'\clearpage' + '\n')

        i += 1

    f.writelines(r'\end{document}' + '\n')
