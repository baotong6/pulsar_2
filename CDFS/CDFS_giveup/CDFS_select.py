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
from tkinter import _flatten
from astropy.stats import poisson_conf_interval
from pwkit.bblocks import tt_bblock as bb

srcid=np.linspace(1,1055,1055)
srcid=srcid.astype(int)

path='/Users/baotong/Desktop/CDFS/'

obs_ID=[581,1431,441,582,2406,2405,2312,1672,2409,2313,2239,8591,9593,9718,8593,8597,8595,8592,8596,
9575,9578,8594,9596,12043,12123,12044,12128,12045,12129,12135,12046,12047,12137,12138,12055,12213,12048,
12049,12050,12222,12219,12051,12218,12223,12052,12220,12053,12054,12230,12231,12227,12233,12232,12234,16183,
16180,16456,16641,16457,16644,16463,17417,17416,16454,16176,16175,16178,16177,16620,16462,17535,17542,16184,
16182,16181,17546,16186,16187,16188,16450,16190,16189,17556,16179,17573,17633,17634,16453,16451,16461,16191,
16460,16459,17552,16455,16458,17677,18709,18719,16452,18730,16185]

path='/Users/baotong/Desktop/CDFS/'


def get_evt_from_epoch(epochnum, source_id):
    src_evt_all = [];
    bkg_evt_all = [];
    epoch_all = []
    for i in range(len(source_id)):
        epoch = np.loadtxt(path + 'txt_all_obs_0.5_8/epoch_src_{0}.txt'.format(source_id[i]))
        useobs = epoch[np.where((epoch[:, 0] > epochnum[0]) & (epoch[:, 0] < epochnum[1]))]
        #         print(useobs)
        use_obsid = useobs[:, 2].astype(int)
        src_evt = [];
        bkg_evt = [];
        for k in range(len(use_obsid)):
            src_evt_k = np.loadtxt(path + 'txt_obs_each/txt_{0}_0.5_8/{1}.txt'.format(use_obsid[k], source_id[i]))
            # bkg_evt_k = np.loadtxt(path + 'txt_obs_each/txt_{0}_0.5_8/{1}_bkg.txt'.format(use_obsid[k], source_id[i]))
            bkg_evt_k = np.loadtxt(path + 'txt_obs_each/txt_{0}_0.5_8/19_bkg.txt'.format(use_obsid[k], source_id[i]))
            if len(src_evt_k) == 0:
                continue
            elif type(src_evt_k[0]) == type(np.array([1.2])[0]):
                src_evt_k = [src_evt_k]

            if len(bkg_evt) == 0:
                bkg_evt = bkg_evt_k

            if len(bkg_evt_k) == 0:
                bkg_evt = bkg_evt
            elif type(bkg_evt_k[0]) == type(np.array([1.2])[0]):
                bkg_evt_k = [bkg_evt_k]
                bkg_evt = np.row_stack((bkg_evt, bkg_evt_k))
            else:
                bkg_evt = np.row_stack((bkg_evt, bkg_evt_k))
            if len(src_evt) == 0:
                src_evt = src_evt_k
            else:
                src_evt = np.row_stack((src_evt, src_evt_k))
        src_evt_all.append(src_evt);
        bkg_evt_all.append(bkg_evt);
        epoch_all.append(useobs)
    return [src_evt_all, bkg_evt_all, epoch_all]

def get_split_txt():
    epoch1=[55e6,100e6]
    epoch2=[300e6,320e6]
    epoch3=[380e6,400e6]
    epoch4=[500e6,580e6]
    EPALL=np.array([epoch1,epoch2,epoch3,epoch4])

    source_id = srcid
    source_id=['XID19_05']
    for i in range(len(EPALL)):
        EVTALL = get_evt_from_epoch(EPALL[i], source_id)
        for j in range(len(source_id)):
            if len(EVTALL[0][j]) == 0:
                np.savetxt(path + 'txt_all_obs_0.5_8_ep{0}/{1}.txt'.format(i + 1, source_id[j]),
                           EVTALL[0][j])
            else:
                np.savetxt(path + 'txt_all_obs_0.5_8_ep{0}/{1}.txt'.format(i + 1, source_id[j]),
                           EVTALL[0][j], fmt="%.7f  %8.3f  %10d")
            if len(EVTALL[1][j]) == 0:
                np.savetxt(path + 'txt_all_obs_0.5_8_ep{0}/{1}_bkg.txt'.format(i + 1, source_id[j]),
                           EVTALL[1][j])
            else:
                np.savetxt(path + 'txt_all_obs_0.5_8_ep{0}/{1}_bkg.txt'.format(i + 1, source_id[j]),
                           EVTALL[1][j], fmt="%.7f  %8.3f  %10d")
            np.savetxt(path + 'txt_all_obs_0.5_8_ep{0}/epoch_src_{1}.txt'.format(i + 1, source_id[j]),
                       EVTALL[2][j], fmt='%15.2f %15.2f %10d %20.2f')

if __name__=='__main__':
    get_split_txt()