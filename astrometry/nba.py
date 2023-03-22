import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import pandas as pd
import sys
import os
from scipy import interpolate
from astropy.wcs import WCS
from astropy.io import fits
from astropy import units as u

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
pointinpaint=np.array(pointinpaint)
freethrow=np.array(freethrow)
drives=np.array(drives)
x=np.zeros(30)
y=np.zeros(30)
z=np.zeros(30)
for i in range(30):
    name=freethrow[i]
    x[i]=i+1
    # y[i]=np.where(freethrow==name)[0]+1
    z[i]=np.where(drives==name)[0]+1
    # plt.scatter(x[i],y[i],label=name)
    plt.scatter(x[i],z[i],label=name)
    plt.text(x[i],z[i],name)
plt.xlabel('Freethrow attemp')
plt.ylabel('Drives')
plt.show()
