#!/bin/bash
# -*- coding: utf-8 -*-
'''
Author: baotong && baotong@smail.nju.edu.cn
Date: 2023-03-02 12:51:08
LastEditors: baotong && baotong@smail.nju.edu.cn
LastEditTime: 2024-03-11 14:46:43
FilePath: /pulsar/astrometry/nba.py
Description: 

Copyright (c) 2023 by baotong, All Rights Reserved. 
'''

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astropy.io import fits
import pandas as pd
import sys
import os
from scipy import interpolate
from astropy.wcs import WCS
from astropy.io import fits
from astropy import units as u
plt.rcParams['font.sans-serif']=['Arial Unicode MS'] #用来正常显示中文标签

font1 = {'family': 'Arial Unicode MS',
         'weight': 'normal',
         'size': 18, }

pointinpaint=['MEM','SAS','LAL','OKC','DEN','NOP','CHA',
      'MIN','HOU','ATL','SAC','CLE','TOR','NYK',
      'JTA','WAS','ORL','IND','CHI','POR','PHI',
      'MIA','PHX','DET','BOS','MIL','BKN','GSW',
      'LAC','DAL']
freethrow=['DET','LAL','DAL','ORL','HOU','NYK','POR','SAC','PHI',
           'NOP','MEM','TOR','JTA','CHA','LAC','IND','MIN','OKC',
           'WAS','CLE','MIA','MIL','CHI','DEN','BOS','ATL','PHX',
           'BKN','SAS','GSW']
drives=['OKC','NYK','IND','JTA','SAS','CHA','HOU','SAC','DET','NOP','ATL','CLE',
        'MIN','ORL','BOS','MEM','TOR','LAL','MIA','DAL','LAC','WAS','MIL','PHI',
        'PHX','BKN','CHI','POR','GSW','DEN']

three2023=['DAL','BOS','SAC','MEM','BKN','GSW','IND','JTA','MIL','ATL','SAS','NYK',
       'WAS','HOU','MIA','TOR','CHI','PHX','POR','LAC','OKC','PHI','CLE','DEN','NOP',
       'CHA','MIN','ORL','LAL','DET']
freethow2023=['PHI','ORL','PHX','MIL','ATL','MIN','LAC','LAL','NOP','GSW','NYK','DAL',
              'OKC','SAC','CHA','TOR','JTA','CLE','BOS','MIA','IND','DET','POR','HOU',
              'MEM','CHI','DEN','WAS','BKN','SAS']
fastbreak2023=['IND','WAS','ATL','BKN','TOR','PHI','LAL','CLE','JTA','OKC','NOP','DAL','DET',
               'NYK','MIN','POR','LAC','CHA','BOS','SAS','ORL','MEM','DEN','CHI','HOU','SAC',
               'PHX','MIL','GSW','MIA']
name2024=np.array(['ATL','BOS','BKN','CHA','CHI','CLE','DAL','DEN','DET','GSW',
          'HOU','IND','LAC','LAL','MEM','MIA','MIL','MIN','NOP','NYK',
          'OKC','ORL','PHI','PHX','POR','SAC','SAS','TOR','UTA','WAS'])
paces2024=[105.8,102.4,101.7,101.7,101.0,101.7,103.6,99.7,104.3,104.3,
           103.2,105.7,101.3,104.9,102.4,100.1,105.2,101.3,101.9,100.1,
           103.8,102.1,102.8,102.5,102.5,104.2,105.0,102.6,104.5,106.3]
defense2024=[1.172,1.084,1.133,1.178,1.120,1.075,1.136,1.113,1.171,1.128,
             1.091,1.156,1.110,1.115,1.104,1.104,1.131,1.056,1.099,1.097,
             1.095,1.088,1.110,1.118,1.136,1.136,1.147,1.146,1.149,1.164]
Fouls2024=[18.8,17.1,19.3,19.1,19.3,18.6,19.0,18.5,21.8,20.5,
          21.4,22.2,19.2,16.5,19.6,18.1,20.0,19.5,19.1,18.2,
          19.5,21.0,20.8,18.5,20.5,20.2,18.1,18.4,19.2,20.2]
FTA2024=[25.2,21.8,21.1,19.0,21.2,20.7,23.7,20.6,22.2,22.1,
         23.5,21.3,23.1,24.9,20.3,22.9,25.7,23.7,24.1,23.1,
         22.1,25.9,26.7,24.8,21.9,21.8,20.2,22.2,23.2,20.9]
OPTFTA2024=[22.8,19.1,22.4,22.7,22.2,22.4,22.6,22.6,26.4,24.1,
            25.9,26.8,22.0,19.3,22.4,20.2,22.9,22.0,21.9,20.4,
            24.0,23.6,24.7,21.9,23.7,23.5,22.3,19.7,21.9,24.3]
OPTFouls2024=[20.6,18.4,18.3,17.5,18.9,19.2,21.5,18.7,18.8,19.3,
              20.1,18.8,19.6,20.4,18.6,19.2,20.1,20.3,19.0,19.8,
              19.4,21.8,20.1,20.8,18.5,18.9,18.6,18.4,20.1,18.3]
Pointsinpaint2024=[51.8,46.3,48.7,49.4,47.2,51.8,46.2,53.0,52.3,47.5,
                   50.8,57.0,50.6,55.6,45.1,46.6,48.8,51.0,51.4,50.5,
                   52.9,52.4,52.9,49.0,45.4,51.7,51.3,53.7,54.5,55.7]
OPT_3PT_percent2024=[33.9,35.9,34.9,34.3,38.1,33.0,33.4,30.9,29.4,33.1,
                     33.0,26.4,33.2,36.1,35.2,36.7,30.7,32.1,35.6,35.1,
                     35.8,31.5,30.8,34.1,29.7,34.5,33.7,33.8,35.8,30.2]
point3rate_2024=np.array([40.3,46.9,41.7,37.9,36.6,42.2,44.5,
                          35.1,35.4,43.2,38.6,38.4,38.4,34.8,
                          43.9,38.7,42.9,37.9,36.9,40.4,38.1,
                          36.3,35.9,37.0,37.3,42.7,39.6,37.0,40.8,38.2])
# Opponent Percent of Points from 3 Pointers
pointinpaint=np.array(pointinpaint)
freethrow=np.array(freethrow)
freethow2023=np.array(freethow2023)
three2023=np.array(three2023)
fastbreak2023=np.array(fastbreak2023)
drives=np.array(drives)
paces2024=np.array(paces2024)
defense2024=np.array(defense2024)
Fouls2024=np.array(Fouls2024)
FTA2024=np.array(FTA2024)

fiveftFGA=np.array([31.2,27.4,29.4,29.4,28.7,30.8,25.8,31.6,31.8,27.8,33.1,
                    33.7,28,32.4,27.9,24.8,28.7,29.5,30.9,28.7,30,32.8,31.9,
                    26.1,31.7,28.4,30.1,31.6,30.9,31.9])
fiveftFGpercent=np.array([60.7,67.8,61.7,61.7,60.1,63.5,67.5,63.8,62.9,
                          66.9,60.6,65.2,64.7,68.1,58.4,61.3,68.3,64,63.7,
                          60.3,65,65.1,61,65.1,57.7,66.3,64.8,65.7,64.6,66.2])
x=np.zeros(30)
y=np.zeros(30)
z=np.zeros(30)
for i in range(30):
    name=name2024[i]
    y[i]=defense2024[i]
    # y[i]=(FTA2024[i])/paces2024[i]*100

    x[i]=(FTA2024[i]-OPTFTA2024[i])/paces2024[i]*100
    # y[i]=OPTFTA2024[i]/paces2024[i]*100
    plt.scatter(x[i],y[i],s=50)
    plt.text(x[i],y[i],name,font1)
plt.xlabel('百回合罚球差值',font1)
plt.ylabel('防守效率',font1)
plt.tick_params(labelsize=18)
plt.show()
