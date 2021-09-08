import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import pandas as pd
import sys
import os

path='/Users/baotong/eSASS/data/raw_data/47_Tuc/'
file='pm00_700011_020_EventList_c001.fits'
evt=fits.open(path+file)[1].data['TIME']
TM_NR=fits.open(path+file)[1].data['TM_NR']
flare_START=[];flare_STOP=[]
for i in range(2,9):
    flare_START.extend(fits.open(path+file)[i].data['START'])
    flare_STOP.extend(fits.open(path+file)[i].data['STOP'])
# print(len(flare_START),len(flare_STOP))

GTI_START_all=[];GTI_STOP_all=[]
for i in range(51,57):
    GTI_START_all.extend(fits.open(path+file)[i].data['START'])
    GTI_STOP_all.extend(fits.open(path+file)[i].data['STOP'])
# print(GTI_START_all-np.sort(GTI_START_all))

for i in (1,2,4,5,6,7):
    GTI_START=fits.open(path+file)['FLAREGTI'+str(i)].data['START']
    GTI_STOP=fits.open(path+file)['FLAREGTI'+str(i)].data['STOP']
    TIME_i=evt[np.where(TM_NR==i)]
    goodtime=0
    for j in range(len(GTI_START)):
        identify=(TIME_i-GTI_START[j])*(GTI_STOP[j]-TIME_i)
        rightevt=len(np.where(identify>0)[0])
        goodtime+=rightevt
    print(len(TIME_i))
    print(goodtime)


# plt.hist(evt,bins=1000)
# plt.plot([GTI_START,GTI_START],[np.zeros(len(GTI_START)),np.zeros(len(GTI_START))+500],'--',color='red')
# plt.plot([GTI_STOP,GTI_STOP],[np.zeros(len(GTI_STOP)),np.zeros(len(GTI_STOP))+500],'--',color='green')
# plt.show()