'''
Author: baotong && baotong@smail.nju.edu.cn
Date: 2023-10-31 20:58:31
LastEditors: baotong && baotong@smail.nju.edu.cn
LastEditTime: 2023-11-06 10:48:36
FilePath: /pulsar/GC/test_prob_map.py
Description: 

Copyright (c) 2023 by baotong, All Rights Reserved. 
'''
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import pandas as pd
from astropy import units as u
from astropy import constants as c
from astropy.coordinates import SkyCoord,match_coordinates_sky
from astropy.coordinates import Angle
from scipy.interpolate import griddata
from matplotlib.patches import Circle
font1 = {'family': 'Normal',
         'weight': 'normal',
         'size': 18, }
def find_nearest_grid_point(x, y, grid_ra, grid_dec):
    x = np.clip(x, np.min(grid_ra), np.max(grid_ra))
    y = np.clip(y, np.min(grid_dec), np.max(grid_dec))
    x_idx = np.argmin(np.abs(grid_ra[0] - x))
    y_idx = np.argmin(np.abs(grid_dec[:,0] - y))

    return x_idx, y_idx

def random_match(x,y,filename,sep,save=0,show=1):
    path='/Users/baotong/Desktop/period_terzan5/CSC/'
    cat=pd.read_csv(path+filename)
    ra=np.array(cat['ra'])
    dec=np.array(cat['dec'])
    ra = ra[np.logical_not(np.isnan(ra))]
    dec = dec[np.logical_not(np.isnan(dec))]
    region_ra_min = np.min(ra)
    region_ra_max = np.max(ra)
    region_dec_min = np.min(dec)
    region_dec_max = np.max(dec)
    # N=int((region_dec_max-region_dec_min)*(region_ra_max-region_ra_min)*3600**2/100**2)
    N=int(600/3)  ## 5 arcsec as grid size
    grid_ra, grid_dec = np.meshgrid(np.linspace(region_ra_min, region_ra_max, N), 
                                    np.linspace(region_dec_min, region_dec_max, N))
    # 创建SkyCoord对象，用于处理坐标
    sources_coords = SkyCoord(ra, dec, unit=(u.degree, u.degree), frame='icrs')
    # 使用np.histogram2d计算点源密度
    density, _, _ = np.histogram2d(ra, dec, bins=(grid_ra[0], grid_dec[:, 0]))
    # 计算每个网格单元格的面积
    ra_bin_width = grid_ra[0, 1] - grid_ra[0, 0]
    dec_bin_width = grid_dec[1, 0] - grid_dec[0, 0]
    area = ra_bin_width * dec_bin_width*3600**2
    # 计算点源密度
    density_map = density / area
    probability_map=1-np.exp(-density_map*(3.14*sep**2))
    # 绘制颜色映射
    # 设置颜色范围
    vmin = 0  # 最小值
    vmax = 10 # 最大值
    # 绘制颜色映射，使用暖色调
    # plt.imshow(density_map, cmap='hot', origin='lower', extent=[region_ra_min, region_ra_max, region_dec_min, region_dec_max], vmin=vmin, vmax=vmax)
    xidx,yidx=find_nearest_grid_point(x, y, grid_ra, grid_dec)
    prob_out=probability_map[xidx,yidx]
    # plt.scatter(x,y,s=100,marker='o')
    # plt.imshow(np.log10(density_map+1e-5), cmap='hot', origin='lower', 
    #            extent=[region_ra_min, region_ra_max, region_dec_min, region_dec_max], vmin=np.log10(vmin+1e-5), 
    #            vmax=np.log10(vmax+1e-5))
    # plt.colorbar(label='Log Point Density')
    plt.imshow(np.log10(probability_map+1e-5), cmap='hot', origin='lower', 
               extent=[region_ra_min, region_ra_max, region_dec_min, region_dec_max], vmin=np.log10(vmin+1e-5), 
               vmax=np.log10(vmax+1e-5))
    plt.colorbar(label='random match possibility')
    plt.xlabel('RA',font1)
    plt.ylabel('Dec',font1)
    plt.tick_params(labelsize=16)
    plt.title(f'M28: Possibility of random match within {sep} arcsecond',font1)
    if save:
        plt.savefig('coord.pdf', bbox_inches='tight', pad_inches=0.05)
    if show:
        plt.show()
    return prob_out


def plot_coord(x,y,save=0,show=0):
    fig, ax = plt.subplots(figsize=(8,8))
    circle1 = Circle((x, y), 2/3600, fill=False, color='blue',linewidth=2, linestyle='--',label=r'2 arcsec')
    circle2 = Circle((x, y), 10/60, fill=False, color='blue',linewidth=2, linestyle='--',label=r'10 arcmin')
    ax.add_patch(circle1)
    ax.add_patch(circle2)
    plt.scatter(pmra,pmdec,s=30,marker='.')
    plt.scatter(x,y,s=100,marker='o')
    if save:
        plt.savefig('coord.pdf', bbox_inches='tight', pad_inches=0.05)
    if show:
        plt.show()
# value = nearest_neighbor_interpolation(grid, x, y)
# print("坐标", (x, y), "的估算位置是:", value)
if __name__=='__main__':
    path='/Users/baotong/Desktop/period_terzan5/CSC/'
    xraysrc=pd.read_csv(path+'match_GC_GAIA_3arcsec_new.csv')
    seq=np.array(xraysrc['seq'])
    GCname=np.array(xraysrc['GCname'])
    ra_x=np.array(xraysrc['ra'])
    dec_x=np.array(xraysrc['dec'])
    ra_GAIA=np.array(xraysrc['ra_GAIA'])
    dec_GAIA=np.array(xraysrc['dec_GAIA'])
    sep_GC_GAIA=np.array(xraysrc['sep_GC_GAIA'])
    prob=[]
    # for i in range(len(ra_x)):
    for i in range(7):
        if GCname[i]=='Tuc' or GCname[i]=='omg':
            prob.append([-1])
            continue
        else:
        # if (not np.isnan(ra_GAIA[i])) and (ra_GAIA[i]>0):
            x=ra_x[i];y=dec_x[i];sep=sep_GC_GAIA[i]
            print(seq[i],GCname[i])
            a=random_match(x,y,filename=f'GAIA_{GCname[i]}.csv',sep=0.5,show=1)
            prob.append([a])
    # np.savetxt(path+'prob_res_0.5.txt',prob,fmt='%.4f')



