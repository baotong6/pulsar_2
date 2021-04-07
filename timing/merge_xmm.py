import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
path='/Users/baotong/xmm/0109110101/txt'
os.chdir(path)
txt_mos1=np.loadtxt('WR46_mos1.txt')
txt_mos2=np.loadtxt('WR46_mos2.txt')
txt_pn=np.loadtxt('WR46_pn.txt')
txt_all=np.vstack((txt_mos1,txt_mos2,txt_pn))
print(txt_all)
txt_all=txt_all[np.lexsort(txt_all[:,::-1].T)]
print(txt_all)
np.savetxt('WR46_all.txt',txt_all)