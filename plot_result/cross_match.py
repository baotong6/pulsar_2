#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import pylab as pl
import string
from datetime import datetime
from scipy.interpolate import lagrange
from scipy import interpolate
from scipy.fftpack import fft,ifft
import scipy.signal as ss
import random
import pandas as pd
import read_csv as data

def getlen(a,b):
    length=((a[0]-b[0])**2+(a[1]-b[1])**2)**0.5
    return length

# for globular cluster and pulsar
def get_input(label):
    mode='fits'
    path = "/Users/baotong/Desktop/period_Tuc/"
    if mode=='table':
        type = ['47Tuc', 'terzan5', 'M28', 'omg_cen']
        ##===============from table=================##
        res = pd.read_excel(path + 'result_0.5_8_all.xlsx', label)
        ra1 = np.array(res['RA'])
        dec1 = np.array(res['DEC'])
        seq1 = np.array(res['seq'])
        period = np.array(res['P_out'])
        ##==========================================##
    elif mode=='fits':
        ##===============from fits=================##
        if label=='47Tuc':
            srclist=fits.open(path + 'xray_properties-592.fits')
            ra1 = srclist[1].data['RAdeg']
            dec1= srclist[1].data['DEdeg']
            seq1 = np.linspace(1, len(ra1), len(ra1)).astype(int)
        if label=='terzan5':
            srclist=fits.open(path + 'terzan5_p50_i5_src_1_2_4_8.fits')
            ra1 = srclist[1].data['RA']
            dec1= srclist[1].data['DEC']
            seq1 = np.linspace(1, len(ra1), len(ra1)).astype(int)
        if label=='M28':
            srclist=fits.open(path + 'M28_p50_i5_src_1_2_4_8.fits')
            ra1 = srclist[1].data['RA']
            dec1= srclist[1].data['DEC']
            seq1 = np.linspace(1, len(ra1), len(ra1)).astype(int)
        ##==========================================##

    ra_dec1=[ra1,dec1]
    textfile=[]
    with open(path + 'pulsar_{0}.txt'.format(label), 'r') as file_to_read:
        while True:
            lines = file_to_read.readline()  # 整行读取数据
            textfile.append(lines)
            if not lines:
                break
                pass
    text = textfile[2:-1]
    text=np.array(text)
    seq2=[];RA2=[];DEC2=[]
    for i in range(len(text)):
        #print(text[i])
        seq,NAME,P0,RA,DEC,PB,MINM,MEDM,ASSOC,useless = [i for i in text[i].split(';')]
        seq2.append(seq);RA2.append(RA);DEC2.append(DEC)
    ra_dec2 = [np.array(RA2).astype(float), np.array(DEC2).astype(float)]
    return [ra_dec1,ra_dec2,seq1,seq2]
#print(get_input('47Tuc'))

def get_NGC6397():
    file1='/Users/baotong/Desktop/period_NGC6397/Kaluzny.txt'
    seq1 = np.loadtxt(file1)[:, 0]
    ra1=np.loadtxt(file1)[:,1]
    dec1=np.loadtxt(file1)[:,2]
    ra_dec1 = [ra1, dec1]

    file2=fits.open('/Users/baotong/Desktop/period_NGC6397/ngc6397_catalog.fits')
    ra2 = file2[1].data['RAdeg']
    dec2 = file2[1].data['DEdeg']
    seq2=file2[1].data['Seq']
    ra_dec2 = [ra2, dec2]
    return [ra_dec1, ra_dec2, seq1, seq2]

def compare_counterpart(ra_dec1,ra_dec2,seq1,seq2):

    offset=1
    ##arcsec

    match=[[] for i in range(len(seq1))]
    for i in range(len(seq1)):
        for j in range(len(seq2)):
            dis=getlen([ra_dec2[0][j],ra_dec2[1][j]],[ra_dec1[0][i],ra_dec1[1][i]])*3600
            if dis<offset:
                match[i].append([seq1[i],seq2[j],dis])
    print(np.sort(match))
    # plt.scatter(x1,y1,color='red')
    # plt.scatter(x2,y2,color='green')
    # plt.show()
    # plt.legend(['X-ray','pulsar'])
    # print(match)
#compare_counterpart(get_NGC6397()[0],get_NGC6397()[1],get_NGC6397()[2], get_NGC6397()[3])
# compare_counterpart(get_input('47Tuc')[0],get_input('47Tuc')[1],get_input('47Tuc')[2],get_input('47Tuc')[3])

def compare_counterpart_long(ra_dec1,ra_dec2,seq1,seq2):
    offset=0.05
    match_ID2=[[] for i in range(len(seq1))]
    for i in range(len(seq1)):
        ra1=ra_dec1[0][i];dec1=ra_dec1[1][i]
        dis=((ra_dec2[0]-ra1)**2+(ra_dec2[1]-dec1)**2)**0.5
        match=seq2[np.where(dis<offset)[0]]
        match_ID2[i]=match
        print(i)

    return match_ID2

if __name__=='__main__':
    path='/Users/baotong/Downloads/'
    file1=np.loadtxt(path+'test0512_02.txt')
    file2=np.loadtxt(path+'test0512_03.txt')
    ra_dec1=[file1[:,0],file1[:,1]]
    ra_dec2 = [file2[:, 0], file2[:, 1]]
    seq1=np.arange(1,len(ra_dec1[0])+1,1)
    seq2=np.arange(1,len(ra_dec2[0])+1,1)
    match=compare_counterpart_long(ra_dec1,ra_dec2,seq1,seq2)
    print(seq2)

