#!/bin/bash
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 14:58:00 2023
@author: baotong
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib import ticker
import os
import pandas as pd
import hawkeye as hawk
import NGC104.plot_pXS as plot_pXS
import GC.localCVs as localCVs
import NGC104.CV_model as CV_model
import GC.SDSSCV as SDSS
import re

pos_all = {'Tuc': [6.0236250, -72.0812833, 3.17 * 60, 3.17 / 8.8 * 60],
           'terzan5': [267.0202083, -24.7790556, 0.72 * 60, 0.72 / 3.4 * 60],
           'M28': [276.1363750, -24.8702972, 1.97 * 60, 1.97 / 8.2 * 60],
           'omg': [201.69700, -47.47947, 5 * 60, 5 / 2.1 * 60],
           'NGC6397': [265.17539, -53.67433, 2.9 * 60, 2.9 / 58 * 60],
           'NGC6752': [287.71713, -59.98455, 1.91 * 60, 1.91 / 11.24 * 60],
           'NGC6266': [255.303333, -30.113722, 0.92 * 60, 0.92 / 4.2 * 60],
           'M30': [325.092167, -23.179861, 1.03 * 60, 1.03 / 17.2 * 60]}

def extract_values(input_string):
    pattern = r'(\w+)=(\d+\.\d+)'
    matches = re.findall(pattern, input_string)
    result = {}
    for match in matches:
        name, value = match
        result[name] = float(value)
    return result

def process_dat_file(file_path,text=None):
    # 使用pandas读取dat文件
    df = pd.read_csv(file_path, delimiter=r'\s+')  # 这里假设dat文件是以制表符分隔的，可以根据实际情况进行修改

    # 在这里编写你需要执行的操作，可以根据需求对DataFrame进行处理
    print(f"Processing file: {file_path}")
    # 查看 DataFrame 的形状（大小）
    Xlum1 = df['Xlum1']
    Xlum2 = df['Xlum2']
    Xlum3 = df['Xlum3']
    Xlum4 = df['Xlum4']
    Porb = df['period_days']
    Mwd = df['star_2_mass']
    Mdot = df ['lg_mtransfer_rate']
    MBAML=df['jdot_ml']
    GRAML=df['jdot_gr']
    # plt.plot(Porb*24,Xlum2,'--',label=text)
    plt.scatter(Porb*24,MBAML,label=text)
    print(np.sort(MBAML))
    return (Porb,Xlum1,Xlum3,Mdot,MBAML,text)
    # plt.plot(Porb*24,Xlum4,'--')

def read_dat_files(folder_path):
    # 获取指定文件夹中所有的dat文件
    dat_files = [file for file in os.listdir(folder_path) if file.endswith('.dat')]
    # 逐个处理每个dat文件

    for dat_file in dat_files:
        namef=extract_values(dat_file)
        print(namef)
        labelname=r'$M_{\rm wd}$'+str(namef['MASS2'])+', Z='+str(namef['Z'])
        if namef['Z']>0.0001:
            file_path = os.path.join(folder_path, dat_file)
            process_dat_file(file_path,text=labelname)

def filtered_line(X, Y,Z):
    # 计算相邻数据点的斜率
    # 找到需要绘制的数据点的索引
    # 构建经过筛选的数据
    idex=np.where((Z>=-1e35))[0]
    print('idex=',len(idex))
    flt_X=np.array(X[idex]);flt_Y=np.array(Y[idex])

    # dx = np.diff(np.log(flt_X))
    # dy = np.diff(np.log(flt_Y))
    #
    # distances = np.sqrt(dx ** 2 + dy ** 2)
    # fltindex=np.where(distances<1)[0]
    # filtered_X=flt_X[fltindex];filtered_Y=flt_Y[fltindex]

    return (flt_X,flt_Y)

if __name__ == '__main__':
    # 指定包含dat文件的文件夹路径
    folder_path = '/Users/baotong/Desktop/period_terzan5/Tracks/'
    # 读取并处理dat文件
    # read_dat_files(folder_path)
    read_dat_files(folder_path)
    # plt.loglog()
    # plt.semilogx()
    # plt.semilogy()
    plt.legend()
    plt.show()
