'''
Author: baotong && baotong@smail.nju.edu.cn
Date: 2023-10-09 08:55:32
LastEditors: baotong && baotong@smail.nju.edu.cn
LastEditTime: 2023-11-02 09:39:19
FilePath: /pulsar/GC/match_CSC_GAIA.py
Description: 

Copyright (c) 2023 by ${git_name_email}, All Rights Reserved. 
'''
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import pandas as pd
import sys
import os
from astropy import units as u
from astropy.coordinates import SkyCoord,match_coordinates_sky
from astroquery.gaia import Gaia
from astropy.coordinates import Angle
Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"

path='/Users/baotong/Desktop/period_terzan5/CSC/'
filename='CSC2.1p_OIR_SDSSspecmatch.csv'
CSC=pd.read_csv(path+filename)
ra_CSC=CSC['ra'];dec_CSC=CSC['dec']
GAIA21P_source_id=CSC['GAIA21P_source_id']
GAIA21P_ra=CSC['GAIA21P_ra']
GAIA21P_dec=CSC['GAIA21P_dec']
Sep_GAIA21P_CSC21P=CSC['Sep_GAIA21P_CSC21P']
GAIA21P_g=CSC['GAIA21P_g']

result_all = pd.read_excel('/Users/baotong/Desktop/period_terzan5/candidate_allGC.xlsx', 'all')
ra_GC = np.array(result_all['ra'])
dec_GC = np.array(result_all['dec'])
seq_GC = np.array(result_all['seq'])
period = np.array(result_all['period_all'])
type = np.array(result_all['GC'])
dist = np.array(result_all['proj_dist'])
counts = np.array(result_all['counts'])
exptime = np.array(result_all['expT'])
L = np.array(result_all['L'])

c1 = SkyCoord(ra=ra_CSC * u.degree, dec=dec_CSC * u.degree)
c2 = SkyCoord(ra=ra_GC * u.degree, dec=dec_GC * u.degree)
idx, d2d, d3d = c2.match_to_catalog_sky(c1)
d2d=Angle(d2d)
d2d=d2d.arcsec
# print(GAIA21P_source_id[idx])
# print(GAIA21P_ra[idx])
# print(GAIA21P_dec[idx])
# print(Sep_GAIA21P_CSC21P[idx])
# print(GAIA21P_g[idx])

ra_GAIA=np.array(GAIA21P_ra[idx]);dec_GAIA=np.array(GAIA21P_dec[idx])
ID_GAIA=np.array(GAIA21P_source_id[idx])
dis_GC_CSC=d2d
ra_CSC=np.array(ra_CSC[idx]);dec_CSC=np.array(dec_CSC[idx])

header = ['seq', 'GCname', 'period','ra','dec','ra_CSC','dec_CSC','sep_GC_CSC','ra_GAIA','dec_GAIA',
          'sep_GC_GAIA','sep_CSC_GAIA','pmra','pmdec','pmra_GC','pmra_GC']
# 将数据转换为 DataFrame
data_list = []
# for i in [1,2]:
for i in range(len(ra_GC)):
    # if np.isnan(ra_GAIA[i]):
    #     row=[seq_GC[i],type[i],period[i],ra_GC[i],dec_GC[i],ra_CSC[i],dec_CSC[i],dis_GC_CSC[i],
    #          0,0,0,0,0,0,0,0]
    #     data_list.append(row)
    #     continue
    # else:
        coord = SkyCoord(ra=ra_GC[i], dec=dec_GC[i], unit=(u.degree, u.degree), frame='icrs')
        width = u.Quantity(3, u.arcsec)
        height = u.Quantity(3, u.arcsec)
        r = Gaia.query_object_async(coordinate=coord, width=width, height=height)
        if len(r)==0:
            row = [seq_GC[i], type[i], period[i], ra_GC[i], dec_GC[i], ra_CSC[i], dec_CSC[i], dis_GC_CSC[i],
                   0, 0, 0, 0, 0, 0, 0, 0]
            data_list.append(row)
            continue
        print(i,len(r),ra_CSC[i],dec_CSC[i])
        c1 = SkyCoord(ra=r[0]['ra'] * u.degree, dec=r[0]['dec'] * u.degree)
        c2 = SkyCoord(ra=ra_GC[i] * u.degree, dec=dec_GC[i] * u.degree)
        cc = SkyCoord(ra_CSC[i] * u.deg, dec_CSC[i] * u.deg, frame='fk5')
        dist1 = c1.separation(c2);dist1=dist1.arcsec
        dist2 = c1.separation(cc);dist2=dist2.arcsec

        row=[seq_GC[i],type[i],period[i],ra_GC[i],dec_GC[i],ra_CSC[i],dec_CSC[i],dis_GC_CSC[i],
             r[0]['ra'],r[0]['dec'],dist1,dist2,r[0]['pmra'],r[0]['pmdec'],0,0]
        data_list.append(row)
        # r.pprint(max_lines=12, max_width=130)
df = pd.DataFrame(data_list, columns=header)
df.to_csv(path+'match_GC_GAIA_3arcsec_new.csv', index=False)





