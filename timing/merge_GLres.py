import numpy as np
import matplotlib.pyplot as plt
from tkinter import _flatten
import sys
import os
import string
path='/Users/baotong/Desktop/period_Tuc/result_bright_src_0.5_8/'

srcid=np.loadtxt('/Users/baotong/Desktop/period_Tuc/result_bright_src_0.5_8/'+'bright_src.txt')
srcid=srcid.astype(int)
# srcid=np.arange(4,888,1)

#print(srcid)
#res=np.loadtxt(path+'result_m03/'+'result_1h_{0}.txt'.format(str(srcid[0])))
res=np.loadtxt(path+'result_3h_10h_{0}.txt'.format(str(srcid[1])))
for i in range(2,len(srcid)):
    if os.path.exists(path+'result_3h_10h_{0}.txt'.format(str(srcid[i]))):
        res_single=np.loadtxt(path+'result_3h_10h_{0}.txt'.format(str(srcid[i])))
        res=np.row_stack((res,res_single))
        np.savetxt(path+ 'all_result_3h_10h.txt', res, fmt='%10d %10.2f %10.5f %15.5f %15.5f %15d %15.5f %15.5f %15.5f')

print(res)
