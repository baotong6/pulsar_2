import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import pandas as pd
import sys
import os
from scipy import interpolate
path='/Users/baotong/Desktop/period_Tuc/'
label_all=['47Tuc','terzan5','M28','omg_cen','NGC6397']
pos_all = [[6.0236250, -72.0812833, 3.17 * 60, 3.17 / 8.8 * 60],  # 47Tuc
           [267.0202083, -24.7790556, 0.72 * 60, 0.72 / 3.4 * 60],  # terzan5
           [276.1363750, -24.8702972, 1.97 * 60, 1.97 / 8.2 * 60],  # M28
           [201.69700, -47.47947, 5 * 60, 5 / 2.1 * 60],  # omega_cen
           [265.17539, -53.67433, 2.9 * 60, 2.9 / 58 * 60],  # NGC 6397
           [287.71713, -59.98455, 1.91*60, 1.91 / 11.24 * 60]]  # NGC 6752
# subplt_label=['151','152','153','154','155']
fig, axes = plt.subplots(1, 5, sharey=True, figsize=(20, 10))
figlabel=[[0,0],[0,1],[0,2],[0,3],[0,4]]
for i in range(0,5):
    label=label_all[i];pos=pos_all[i]
    res=pd.read_excel(path+'result_0.5_8_all.xlsx',label)
    # res=np.array(res)[np.where(res['type']=='CV')[0]]
    # res_longP=res[np.where(res['P_out']>3.2*3600)];res_shortP=res[np.where(res['P_out']<3.2*3600)]
    ra = np.array(res['RA'])
    dec = np.array(res['DEC'])
    seq = np.array(res['seq'])
    period = np.array(res['P_out'])
    L = np.array(res['L'])
    Lmin = np.array(res['Lmin'])
    Lmax = np.array(res['Lmax'])
    type = np.array(res['type'])
    distance = (((ra - pos[0]) * 3600) ** 2 + ((dec - pos[1]) * 3600) ** 2) ** 0.5
    k=0
    while k < len(seq):
        if (type[k] != 'CV')&(type[k] != 'shortCV'):
            # if period[i]<4500:
            ra = np.delete(ra, k)
            dec = np.delete(dec, k)
            seq = np.delete(seq, k)
            period = np.delete(period, k)
            type = np.delete(type, k)
            L = np.delete(L, k);
            Lmin = np.delete(Lmin, k);
            Lmax = np.delete(Lmax, k)
            distance=np.delete(distance,k)
        else:
            k += 1;

    # print(distance)
    distance_longP=distance[np.where(period>3*3600)[0]]
    distance_shortP = distance[np.where(period < 3 * 3600)[0]]
    # ax_temp=axes[figlabel[i][0],figlabel[i][1]]
    ax_temp = axes[i]
    # ax_temp.subplot(subplt_label[i])
    # ax_temp.title(label)
    ax_temp.semilogx()
    # ax_temp.xlabel('r (arcsec)')
    # ax_temp.ylabel('CDF')
    ax_temp.set_title(label)
    ax_temp.set_xlabel('r (arcsec)')
    if (i == 0 ): ax_temp.set_ylabel('CDF')
    ax_temp.plot([pos[-1],pos[-1]],[0,1],'--')
    ax_temp.hist(distance_longP,density=True, bins=20,cumulative='True',histtype='step',color='red')
    ax_temp.hist(distance_shortP, density=True, bins=20,cumulative='True', histtype='step',color='green')
    ax_temp.legend(['core-radius','P>3h','P<3h'])
    if i <4: ax_temp.legend_.remove()
plt.subplots_adjust(wspace=0, hspace=0)
plt.show()