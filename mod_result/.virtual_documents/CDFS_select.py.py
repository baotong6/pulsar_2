get_ipython().run_line_magic("matplotlib", " widget      ")


import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import pandas as pd
import sys
import os
from scipy import interpolate
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits


srcid=np.linspace(1,1055,1055)
srcid=srcid.astype(int)
srcid=srcid


path='/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8/'


def plot_longT_V(data_file,epoch_file):
    epoch_info = np.loadtxt(epoch_file)
    t_start = epoch_info[:, 0]
    t_end = epoch_info[:, 1]
    obsID = epoch_info[:, 2]
    expT = epoch_info[:, 3]
    time = np.loadtxt(data_file)
    cts=[]
    for i in range(len(obsID)):
        cts.append(len(np.where(time[:,2]==obsID[i])[0]))
    cts=np.array(cts)
    CR=cts/expT
#     print(obsID[np.where(CR>0.0006)])
    if np.min(CR)get_ipython().getoutput("=0:")
        VI=np.max(CR)/np.min(CR)
    else:VI=0
    plt.title(dataname[0:-4]+', VI={0}'.format(VI))
    plt.scatter(t_start,CR,marker='+')
    plt.show()
    


def get_bright_id(cts_limit):
    cts_num=[]
    for i in range(len(srcid)):
        dataname='{0}.txt'.format(srcid[i])
        data_file=path+dataname
        time=np.loadtxt(data_file)[:,0]
        cts_num.append(len(time))
    cts_num=np.array(cts_num)
    return srcid[np.where(cts_num>cts_limit)]


get_bright_id(1000)


bright_id=get_bright_id(1000)
for id in bright_id:
    dataname='{0}.txt'.format(id)
    data_file=path+dataname
    print(data_file)
    time=np.loadtxt(data_file)[:,0]
    plot_longT_V(path + dataname, path + 'epoch_src_' + dataname)
    plt.show()









