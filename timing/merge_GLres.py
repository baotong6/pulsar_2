import numpy as np
import matplotlib.pyplot as plt
from tkinter import _flatten
import sys
import os
import string
path='/Users/baotong/Desktop/period_M31XRB/result_ssr/result_HRC/'

# srcid=np.loadtxt('/Users/baotong/Desktop/period_Tuc/result_bright_src_0.5_8/'+'bright_src.txt')
# srcid=srcid.astype(int)
srcid=np.arange(2,215,1)

res=np.loadtxt(path+'result_20ks_{0}.txt'.format(str(srcid[0])))
# res=np.loadtxt(path+'result_10s_{0}.txt'.format(str(srcid[6])))
for i in range(0,len(srcid)):
    if os.path.exists(path+'result_20ks_{0}.txt'.format(str(srcid[i]))):
        res_single=np.loadtxt(path+'result_20ks_{0}.txt'.format(str(srcid[i])))
        res=np.row_stack((res,res_single))
# res=np.column_stack((srcid,res))
np.savetxt(path+ 'all_result_20ks.txt', res, fmt='%10d %10.2f %10.5f %10.10f %10.5f %10d %10.10f %10.10f %10d')

print(res)
