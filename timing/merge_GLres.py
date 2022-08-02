import numpy as np
import matplotlib.pyplot as plt
from tkinter import _flatten
import sys
import os
import string
path='/Users/baotong/Desktop/period_NGC6304/result_GL/'

# srcid=np.loadtxt('/Users/baotong/Desktop/period_Tuc/result_bright_src_0.5_8/'+'bright_src.txt')
# srcid=srcid.astype(int)
srcid=np.arange(1,183,1)

res=np.loadtxt(path+'result_3k_{0}.txt'.format(str(srcid[0])))
# res=np.loadtxt(path+'result_10s_{0}.txt'.format(str(srcid[6])))
for i in range(1,len(srcid)):
    print(i)
    if os.path.exists(path+'result_3k_{0}.txt'.format(str(srcid[i]))):
        res_single=np.loadtxt(path+'result_3k_{0}.txt'.format(str(srcid[i])))
        res=np.row_stack((res,res_single))
# res=np.column_stack((srcid,res))
np.savetxt(path+ 'all_result_3k.txt', res, fmt='%10d %10.2f %10.5f %10.5f %10.5f %10d %10.5f %10.5f %10d')

print(res)
