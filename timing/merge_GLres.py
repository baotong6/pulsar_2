import numpy as np
import matplotlib.pyplot as plt
from tkinter import _flatten
import sys
import os
import string
path='/Users/baotong/eSASS/data/raw_data/47_Tuc/result_GL/result_psf50_0.5_5/'

#srcid=np.loadtxt(path+'txt_all_obs_I/'+'cand_other500.txt')
# candid=np.loadtxt(path+'txt_all_obs_0.5_2/'+'cand_id.txt')
# candid=candid.astype(int)
srcid=np.arange(1,889,1)

#print(srcid)
#res=np.loadtxt(path+'result_m03/'+'result_1h_{0}.txt'.format(str(srcid[0])))
res=np.loadtxt(path+'result_0.1h_1h_{0}.txt'.format(str(srcid[0])))
for i in range(1,len(srcid)):
    if os.path.exists(path+'result_0.1h_1h_{0}.txt'.format(str(srcid[i]))):
        res_single=np.loadtxt(path+'result_0.1h_1h_{0}.txt'.format(str(srcid[i])))
        res=np.row_stack((res,res_single))
        np.savetxt(path+ 'all_result_0.1h_1h.txt', res, fmt='%10d %10.2f %10.5f %15.5f %15.5f %15d %15.5f %15.5f %15.5f')

print(res)
