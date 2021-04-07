import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import pandas as pd
import sys
import os
from scipy import interpolate

def getlen(a,b):
    length=((a[0]-b[0])**2+(a[1]-b[1])**2)**0.5
    return length

def compare_counterpart(ra_dec1,ra_dec2,seq1,seq2):
    offset=1
    match=[[] for i in range(len(seq1))];
    seq2_all=[];ra2_all=[];dec2_all=[]
    for i in range(len(seq1)):
        for j in range(len(seq2)):
            dis=getlen([ra_dec2[0][j],ra_dec2[1][j]],[ra_dec1[0][i],ra_dec1[1][i]])*3600
            if dis<offset:
                match[i].append([seq1[i],seq2[j],dis])
                seq2_all.append(seq2[j])
                ra2_all.append(ra_dec2[0][j])
                dec2_all.append(ra_dec2[1][j])
    return [match,[ra2_all,dec2_all,seq2_all]]

def plot_color_color(catalog_name,label):
    path='/Users/baotong/Desktop/period_Tuc/'
    file1=np.loadtxt(path+'UV_cat/'+catalog_name)
    file1=file1[np.where((file1[:,5]>-990)&(file1[:,6]>-990))]

    # res = pd.read_excel(path + 'result_0.5_8_all.xlsx', label)
    res=fits.open(path+'xray_properties-592.fits')[1].data

    # ra_xray=res['RA'];dec_xray=res['DEC'];ID_xray=res['seq']
    ra_UV=file1[:,3];dec_UV=file1[:,4];ID_UV=file1[:,0]
    ra_xray = res['RAdeg'];dec_xray = res['DEdeg'];ID_xray = res['Seq']

    f300=file1[:,5]
    f390=file1[:,6]
    [match,info2_all]=compare_counterpart([ra_xray,dec_xray],[ra_UV,dec_UV],ID_xray,ID_UV)
    # print(info2_all)
    plt.scatter(f300-f390,f390,marker='x',color='green',s=10)
    for i in range(len(info2_all[2])):
        index=np.where(ID_UV==info2_all[2][i])
        plt.scatter(f300[index]-f390[index],f390[index],marker='o',color='red',s=50,linewidths=1,edgecolors='black')
    print(len(info2_all[2]))
    plt.gca().invert_yaxis()
    plt.title(label)
    plt.xlabel('f300-f390')
    plt.ylabel('f390')
    plt.savefig(path+'UV_cat/{0}_HST_UV.eps'.format(label))
    plt.show()


plot_color_color('47Tuc_hst_12950_01_wfc3_uvis_multiwave_daophot_trm.cat','47Tuc')