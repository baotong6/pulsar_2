import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import string
from datetime import datetime
from scipy.interpolate import lagrange
from scipy import interpolate
from scipy.optimize import curve_fit
import pandas as pd
import linecache
from astropy.stats import poisson_conf_interval
font1 = {'family': 'Normal',
         'weight': 'normal',
         'size': 18, }
font2 = {'family': 'Normal',
         'weight': 'normal',
         'size': 16, }
star_ID=np.array([991,950,925,910,872,780,518,457,398,319,153,64])
catalog=fits.open('/Users/baotong/Desktop/CDFS/7Ms_catalog.fit')
RA=catalog[1].data['RAJ2000'];DEC=catalog[1].data['DEJ2000']
SEQ=catalog[1].data['Seq']
def get_photometry_obs(path,id):
    srcevt = np.loadtxt(path + '{0}.txt'.format(id))
    bkgevt = np.loadtxt(path + '{0}_bkg.txt'.format(id))
    epoch = np.loadtxt(path + 'epoch_src_{0}.txt'.format(id))
    t_start = epoch[:, 0]
    t_end = epoch[:, 1]
    obsID = epoch[:, 2]
    expT = epoch[:, 3]
    ra = RA[np.where(SEQ == id)][0]
    dec = DEC[np.where(SEQ == id)][0]
    cts = [];
    bkg_cts = []
    for k in range(len(obsID)):
        cts.append(len(np.where(srcevt[:, 2] == obsID[k])[0]))
        bkg_cts.append(len(np.where(bkgevt[:, 2] == obsID[k])[0]))
    cts = np.array(cts);
    bkg_cts = np.array(bkg_cts)
    bkg_cts = bkg_cts / 12.
    net_cts = cts - bkg_cts
    net_cts[np.where(net_cts < 0)] = 0
    net_cts_all = np.sum(net_cts)

    return [net_cts,ra,dec,expT,t_start]

def get_obs_lc():
    for i in range(len(star_ID)):
        [net_cts, ra, dec,expT,t_start]=get_photometry_obs(path,star_ID[i])
        net_cts_all=np.sum(net_cts)

        y2_err=np.array(poisson_conf_interval(net_cts,interval='frequentist-confidence'))
        y2_err[0]=net_cts-y2_err[0]
        y2_err[1]=y2_err[1]-net_cts

        CR = net_cts / expT
        err_CR=y2_err/expT
        VI=np.max(CR)/np.min(CR[np.where(CR>0)])

        # plt.title('Star {0};VI {1}'.format(star_ID[i],VI))
        plt.figure(1,(9,6))
        plt.title('Star {0}; net_counts={1},({2},{3})'.format(star_ID[i],int(net_cts_all),ra,dec))
        # plt.plot(t_start, CR, marker='+')
        # plt.ylim(1e-7,5e-3)
        t_start = t_start / 86400 + 2449352.5 - 2400000.5
        plt.errorbar(t_start, CR*1e4, yerr=err_CR*1e4, fmt='.', capsize=1, elinewidth=1, ecolor='red')
        plt.xlabel('MJD',font2)
        plt.ylabel('Counts rate ($10^{-4}$ counts/s)',font2)
        # plt.semilogy()
        plt.tick_params(labelsize=16)
        plt.savefig('/Users/baotong/Desktop/CDFS/fig_LC/'+str(star_ID[i])+'.pdf',bbox_inches='tight', pad_inches=0.0)
        # plt.show()
        plt.close()

def get_epoch_lc():
    for i in range(len(star_ID)):
        figlabel = [[0, 0], [0, 1], [1, 0], [1, 1]]
        fig, axes = plt.subplots(2, 2, figsize=(9, 6))
        for ep in range(4):
            path = '/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep{0}/'.format(ep + 1)
            [net_cts, ra, dec, expT, t_start] = get_photometry_obs(path, star_ID[i])

            # net_cts_all=np.sum(net_cts)
            y2_err = np.array(poisson_conf_interval(net_cts, interval='frequentist-confidence'))
            y2_err[0] = net_cts - y2_err[0]
            y2_err[1] = y2_err[1] - net_cts

            CR = net_cts / expT
            err_CR = y2_err / expT
            VI = np.max(CR) / np.min(CR[np.where(CR > 0)])

            ax_temp = axes[figlabel[ep][0], figlabel[ep][1]]
            print(figlabel[ep][0], figlabel[ep][1])

            # plt.title('Star {0}; Epoch {1})'.format(star_ID[i], ep+1))
            # plt.plot(t_start, CR, marker='+')
            # plt.ylim(1e-7,5e-3)
            t_start = t_start / 86400 + 2449352.5 - 2400000.5
            ax_temp.errorbar(t_start, CR * 1e4, yerr=err_CR * 1e4, fmt='.', capsize=1, elinewidth=1, ecolor='red')
            ax_temp.text(t_start[0],np.max(CR*1e4),'Epoch {0}'.format(ep+1))
            if (ep == 2 or ep== 3): ax_temp.set_xlabel('MJD', font2)
            if (ep == 0): ax_temp.set_ylabel('Counts rate ($10^{-4}$ counts/s)', font2)

            if (ep==0):ax_temp.set_title('Star {0}'.format(star_ID[i]))
            # plt.semilogy()
            ax_temp.tick_params(labelsize=16)
        plt.savefig('/Users/baotong/Desktop/CDFS/fig_LC/' + str(star_ID[i]) + '_epoch.pdf', bbox_inches='tight',
                        pad_inches=0.0)
        # plt.show()
        plt.close()
def get_EP_lc():
    for i in range(len(star_ID)):
        CR_ALL=[];err_CR_ALL=[];T_mean=[]
        for ep in range(4):
            path = '/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep{0}/'.format(ep + 1)
            [net_cts, ra, dec, expT, t_start] = get_photometry_obs(path, star_ID[i])
            net_cts_all=np.sum(net_cts)
            expT_all=np.sum(expT)

            y2_err = np.array(poisson_conf_interval(net_cts_all, interval='frequentist-confidence'))
            y2_err[0] = net_cts_all - y2_err[0]
            y2_err[1] = y2_err[1] - net_cts_all
            t_start_mean=np.mean(t_start)
            t_start_mean = t_start_mean / 86400 + 2449352.5 - 2400000.5

            CR = net_cts_all / expT_all
            err_CR = y2_err / expT_all
            CR*=1e4;err_CR*=1e4
            CR_ALL.append(CR);err_CR_ALL.append(err_CR);T_mean.append(t_start_mean)
        err_CR_ALL=np.array(err_CR_ALL)
        err_CR_ALL=err_CR_ALL.reshape(2,4)
        plt.figure(1, (9, 6))
        plt.title('Star {0}'.format(star_ID[i]))
        plt.errorbar(T_mean,CR_ALL,yerr=err_CR_ALL, fmt='.', capsize=1, elinewidth=1, ecolor='red')
        plt.xlabel('MJD',font2)
        plt.ylabel('Counts rate ($10^{-4}$ counts/s)',font2)
        # plt.semilogy()
        plt.tick_params(labelsize=16)
        plt.savefig('/Users/baotong/Desktop/CDFS/fig_LC/'+str(star_ID[i])+'_EP.pdf',bbox_inches='tight', pad_inches=0.0)
        # plt.show()
        plt.close()

def plot_src_certain():
    time = [];
    bkg_time = [];
    energy = [];
    bkg_energy = []
    id='872'
    path = '/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8_ep4/'
    srcevt = np.loadtxt(path + '{0}.txt'.format(id))
    bkgevt = np.loadtxt(path + '{0}_bkg.txt'.format(id))
    epoch = np.loadtxt(path + 'epoch_src_{0}.txt'.format(id))

    t_start = epoch[:, 0]
    t_end = epoch[:, 1]
    obsID = epoch[:, 2]
    expT = epoch[:, 3]
    useid=obsID[42:43]
    print(useid)
    for id_k in range(len(useid)):
        time = np.concatenate((time, srcevt[:, 0][np.where(srcevt[:, -1] == useid[id_k])]))
        energy = np.concatenate((energy, srcevt[:, 1][np.where(srcevt[:, -1] == useid[id_k])]))
        bkg_time = np.concatenate((bkg_time, bkgevt[:, 0][np.where(bkgevt[:, -1] == useid[id_k])]))
        bkg_energy = np.concatenate((bkg_energy, bkgevt[:, 1][np.where(bkgevt[:, -1] == useid[id_k])]))
    plt.figure(1,(9,6))
    plt.title('Star 872; obsID:17677;binsize=1000s',font2)
    plt.xlabel('Time(s)',font2)
    plt.ylabel('Counts/bin',font2)
    plt.tick_params(labelsize=16)
    plt.hist(time-time[0],bins=int(expT[42:43]/1000),histtype='step',color='red')
    plt.show()
    # plt.savefig()
if __name__=='__main__':
    path = '/Users/baotong/Desktop/CDFS/txt_all_obs_0.5_8/'
    # get_obs_lc()
    # get_epoch_lc()
    # get_EP_lc()
    plot_src_certain()