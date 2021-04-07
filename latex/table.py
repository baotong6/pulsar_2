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
path_table = '/Users/baotong/Desktop/period_Tuc/'
path_out='/Users/baotong/Desktop/aas/pCV_GC/'
type=['47Tuc','terzan5','M28','omg_cen']
for k in range(len(type)):
    result=pd.read_excel(path_table + 'result_0.5_8_all.xlsx', type[k])
    ID=result['seq']
    for i in range(len(ID)):
        print(result['L'][i])
        with open(path_out+'Table_src.tex','a+') as f:
            #f.writelines(r'2 & 267.76657 &	-29.57529 & 3820.83 & 0.99222 & 902 &116.4 & 1.94 &-  &- & IP?' + '\n')
            f.writelines(r'{0} & {1:.5f} & {2:.5f} & {3:.2f} & {4} & {5:.2f} & {6:.2f} &  ${7:.2f}^(+{8:.2f})_({9:.2f})$ &-' .format(
                str(type[k])+' '+str(int(i+1)),result['RA'][i],result['DEC'][i],
            result['P_out'][i],result['counts'][i],result['counts_B'][i],result['VI'][i],
            result['L'][i]/1e31,result['Lmax'][i]/1e31-result['L'][i]/1e31,result['Lmin'][i]/1e31-result['L'][i]/1e31)+'\n'+'\\\\'+'\n')
            #f.writelines(r'\end{document}' + '\n')
    f.close()