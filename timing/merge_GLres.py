import numpy as np
import matplotlib.pyplot as plt
from tkinter import _flatten
import sys
import os
import string
path='/Users/baotong/Desktop/period_Tuc/result_startover/result_all_obs_p90_0.5_8/'

srcid=np.loadtxt('/Users/baotong/Desktop/period_Tuc/result_bright_src_0.5_8/'+'bright_src.txt')
srcid=srcid.astype(int)
# srcid=np.arange(1,197,1)

#print(srcid)
#res=np.loadtxt(path+'result_m03/'+'result_1h_{0}.txt'.format(str(srcid[0])))
res=np.loadtxt(path+'result_10s_{0}.txt'.format(str(srcid[6])))
for i in range(1,len(srcid)):
    if os.path.exists(path+'result_10s_{0}.txt'.format(str(srcid[i]))):
        res_single=np.loadtxt(path+'result_10s_{0}.txt'.format(str(srcid[i])))
        res=np.row_stack((res,res_single))

res=np.column_stack((srcid,res))
np.savetxt(path+ 'all_result_10s.txt', res, fmt='%10d %10d %10.5f %10.5f %15.5f %15.5f %15d %15.5f %15.5f %15.5f')

print(res)
