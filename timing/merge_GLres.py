import numpy as np
import matplotlib.pyplot as plt
from tkinter import _flatten
import sys
import os
import string
path='/Users/baotong/Desktop/CDFS/result_all_obs_0.5_8_ep4/'

#srcid=np.loadtxt(path+'txt_all_obs_I/'+'cand_other500.txt')
# candid=np.loadtxt(path+'txt_all_obs_0.5_2/'+'cand_id.txt')
# candid=candid.astype(int)
srcid=np.linspace(1,1055,1055)
srcid=srcid.astype(int)
srcid=srcid
#print(srcid)
#res=np.loadtxt(path+'result_m03/'+'result_1h_{0}.txt'.format(str(srcid[0])))
res=np.loadtxt(path+'result_1h_3h_{0}.txt'.format(str(srcid[2])))
for i in range(3,len(srcid)):
    if os.path.exists(path+'result_1h_3h_{0}.txt'.format(str(srcid[i]))):
        res_single=np.loadtxt(path+'result_1h_3h_{0}.txt'.format(str(srcid[i])))
        res=np.row_stack((res,res_single))
        np.savetxt(path+ 'all_result_1h_3h.txt', res, fmt='%10d %10.2f %10.5f %10.5f %10.5f %10d %10.5f %10.5f %10.5f')

print(res)
