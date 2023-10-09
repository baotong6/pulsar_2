import numpy as np
import os
import matplotlib.pyplot as plt
import glob
import subprocess
import csv
import pandas as pd
import pandas as pd
from astropy.stats import poisson_conf_interval
import math
from astropy.io import fits

list3=['name','R','Rl','Ru','HR','HRl','HRu']
list4=[]

BEHR=os.path.realpath("/Users/baotong/BEHR/BEHR")


def input_parsfromtxt_lc(path,srcid,ebands,ebandh):
    srcevt=np.loadtxt(path+f'{srcid}.txt')
    bkgevt=np.loadtxt(path+f'{srcid}_bkg.txt')
    # obsidlist=srcevt[:,3]
    obsidlist=[14197., 14198., 13825., 13826., 13827., 13828., 14195., 14196.]
    if srcevt.ndim== 1:srcevt = np.array([srcevt])
    energy_ranges = [ebands ,ebandh]

    def extract_counts(data,energy_ranges,N):
        time_column = data[:, 0]
        energy_column = data[:, 1]

        # 计算每个N秒间隔内落在两个能量区间内的光子数
        time_intervals = np.arange(np.min(time_column), np.max(time_column) + 1, N)
        photon_counts = np.zeros((len(energy_ranges), len(time_intervals) - 1), dtype=int)

        for i, (min_energy, max_energy) in enumerate(energy_ranges):
            for j in range(len(time_intervals) - 1):
                mask = (energy_column >= min_energy) & (energy_column <= max_energy) & (
                            time_column >= time_intervals[j]) & (time_column < time_intervals[j + 1])
                photon_counts[i, j] = np.sum(mask)
        return (photon_counts,time_intervals)

    all_area = np.loadtxt(path + '125_backscale.txt')
    obsid = all_area[:, 0];
    srcarea = all_area[:, 1];
    bkgarea=all_area[:, 2];
    for i in range(len(obsidlist)):
        obsid=obsidlist[i]
        srcevt_use=srcevt[np.where(srcevt[:,2]==obsid)]
        bkgevt_use=bkgevt[np.where(bkgevt[:,2]==obsid)]
        src,src_time=extract_counts(srcevt_use,energy_ranges,3000)
        bkg,bkg_time=extract_counts(bkgevt_use,energy_ranges,3000)
        # x=(src_time[:-1]+src_time[1:])/2
        # y=(src[0]-src[1])/(src[0]+src[1])
        # xerr=np.zeros(len(x))+src_time[1]-src_time[0]
        #
        # plt.errorbar(x,y,xerr=xerr)
        # plt.show()
        softsrc=src[0];softbkg=bkg[0];hardsrc=src[1];hardbkg=bkg[1]
        if len(softbkg)< len(softsrc):
            softbkg=np.concatenate((softbkg,np.zeros(len(softsrc)-len(softbkg))))
        if len(hardbkg)< len(hardsrc):
            hardbkg=np.concatenate((hardbkg,np.zeros(len(hardsrc)-len(hardbkg))))
        softarea = bkgarea[i] / srcarea[i];
        hardarea = softarea
        for k in range(len(softsrc)):
            p1 = "softsrc=" + str(int(softsrc[k]))
            p2 = "hardsrc=" + str(int(hardsrc[k]))
            p3 = "softbkg=" + str(int(softbkg[k]))
            p4 = "hardbkg=" + str(int(hardbkg[k]))
            p5 = "softarea=" + str(softarea)
            p6 = "hardarea=" + str(hardarea)
            # p9="softeff="+softeff
            # p10="hardeff="+hardeff
            os.chdir(path + 'HR/125_BEHR/')
            filename = f'behr_{int(obsid)}_{k}'
            p7 = "output=" + filename
            p8 = "outputHR=True"
            p11 = "level=68.0"

            subprocess.call([BEHR, p1, p2, p3, p4, p5, p6, p7, p8, p11])
            str_tmp = "{0:10d} {1:10d} {2:10f} {3:10f} {4:10f} {5:10f} {6:15.5f} {7:15.1f}".format(int(filename[5:10]),int(filename[11:]), softsrc[k],
                                                                                           hardsrc[k], softbkg[k], hardbkg[k],
                                                                                           softarea, hardarea)
            if not os.path.exists(filename + '.txt'):
                return None
            else:
                with open(f'all_inputpars_{int(obsid)}.txt', 'a+') as f2:
                    f2.writelines(str_tmp + '\n')
                with open(filename + '.txt', "r+") as f:
                    data = f.readlines()
                    R = data[1].split('\t')[2]
                    R_lo = data[1].split('\t')[4]
                    R_up = data[1].split('\t')[5]
                    HR = data[2].split('\t')[2]
                    HR_lo = data[2].split('\t')[4]
                    HR_up = data[2].split('\t')[5]
                    list4.append([filename, R, R_lo, R_up, HR, HR_lo, HR_up])
                    print(filename + ':' + HR + ' ' + HR_lo + ' ' + HR_up)
                test = pd.DataFrame(columns=list3, data=list4)
                test.to_csv(f'{srcid}_HR_sup_{int(obsid)}.csv', index=False)

        # print(softsrc,softbkg,hardsrc,hardbkg)
    # if len(bkgevt)==0 or bkgevt.ndim==1:
    #     softbkg=0;hardbkg=0
    # elif bkgevt.ndim==1:
    #     bkgevt = np.array([bkgevt])
    # else:
    #     softbkg=len(np.where((bkgevt[:,1]<ebands[1])&(bkgevt[:,1]>ebands[0]))[0])
    #     hardbkg=len(np.where((bkgevt[:,1]<ebandh[1])&(bkgevt[:,1]>ebandh[0]))[0])

def plot_125_HR_lc(obsid):
    file=pd.read_csv(path+'HR/125_BEHR/'+'125_HR_sup_14196.csv')
    pars = np.loadtxt(path+'HR/125_BEHR/'+f'all_inputpars_{obsid}.txt')
    S = pars[:, 2] - pars[:, 4] / pars[:, 6];
    H = pars[:, 3] - pars[:, 5] / pars[:, 7]
    seq = pars[:, 0]
    seq = seq.astype('int')
    S[np.isnan(S)] = 0;
    H[np.isnan(H)] = 0
    name=file['name']
    print(len(name))
    useid=[]
    for i in range(len(name)):
        if int(name[i][5:10])==int(obsid):
            useid.append(i)
    HR = file['HR'][useid];
    HRl = file['HRl'][useid];
    HRu = file['HRu'][useid]
    SH = H + S
    SH_1sigma = poisson_conf_interval(SH, interval='frequentist-confidence').T
    SHl = SH_1sigma[:, 0];
    SHu = SH_1sigma[:, 1]

    time=pars[:,1]
    # plt.scatter(HR,SH,label=gc)
    ##==plot HR-S+H diagram==##\
    print(len(HR),len(SH))
    fig, (ax1, ax2) = plt.subplots(ncols=1, nrows=2, sharex='all', gridspec_kw={'height_ratios': [2, 1]},
                                   figsize=(15, 10))
    ax1.errorbar(x=time*3000+1500, y=HR, yerr=[HR - HRl, HRu - HR], xerr=1500, fmt='.', capsize=3, elinewidth=1,
                 label=obsid,color='blue')
    ax2.errorbar(x=time*3000+1500, y=SH, yerr=[SH - SHl, SHu - SH], xerr=1500, fmt='.', capsize=3, elinewidth=1,
                 label=obsid,color='red')
    plt.show()

def run_one_BEHR(srcid):
    file=np.loadtxt(path+'HR/125_pfold_BEHR/'+'hs.txt')
    hardsrc=file[0];hardbkg=file[1];softsrc=file[2];softbkg=file[3];
    softarea=10.53;hardarea=10.53
    for k in range(len(softsrc)):
        p1 = "softsrc=" + str(int(softsrc[k]))
        p2 = "hardsrc=" + str(int(hardsrc[k]))
        p3 = "softbkg=" + str(int(softbkg[k]))
        p4 = "hardbkg=" + str(int(hardbkg[k]))
        p5 = "softarea=" + str(softarea)
        p6 = "hardarea=" + str(hardarea)
        # p9="softeff="+softeff
        # p10="hardeff="+hardeff
        os.chdir(path + 'HR/125_pfold_BEHR/')
        filename = f'behr_{int(srcid)}_{k}'
        p7 = "output=" + filename
        p8 = "outputHR=True"
        p11 = "level=68.0"
        subprocess.call([BEHR, p1, p2, p3, p4, p5, p6, p7, p8, p11])
        str_tmp = "{0:10d} {1:10d} {2:10f} {3:10f} {4:10f} {5:10f} {6:15.5f} {7:15.1f}".format(int(filename[5:8]),
                                                                                               int(filename[9:]),
                                                                                               softsrc[k],
                                                                                               hardsrc[k], softbkg[k],
                                                                                               hardbkg[k],
                                                                                               softarea, hardarea)
        if not os.path.exists(filename + '.txt'):
            return None
        else:
            with open(f'all_inputpars_{int(srcid)}.txt', 'a+') as f2:
                f2.writelines(str_tmp + '\n')
            with open(filename + '.txt', "r+") as f:
                data = f.readlines()
                R = data[1].split('\t')[2]
                R_lo = data[1].split('\t')[4]
                R_up = data[1].split('\t')[5]
                HR = data[2].split('\t')[2]
                HR_lo = data[2].split('\t')[4]
                HR_up = data[2].split('\t')[5]
                list4.append([filename, R, R_lo, R_up, HR, HR_lo, HR_up])
                print(filename + ':' + HR + ' ' + HR_lo + ' ' + HR_up)
            test = pd.DataFrame(columns=list3, data=list4)
            test.to_csv(f'{srcid}_HR_pfold.csv', index=False)


def plot_125_HR_pfold():
    file=pd.read_csv(path+'HR/125_pfold_BEHR/'+'125_HR_pfold.csv')
    pars = np.loadtxt(path+'HR/125_pfold_BEHR/'+f'all_inputpars_125.txt')
    S = pars[:, 2] - pars[:, 4] / pars[:, 6];
    H = pars[:, 3] - pars[:, 5] / pars[:, 7]
    seq = pars[:, 0]
    seq = seq.astype('int')
    S[np.isnan(S)] = 0;
    H[np.isnan(H)] = 0
    name=file['name']
    print(len(name))
    HR = file['HR']
    HRl = file['HRl']
    HRu = file['HRu']
    SH = H + S
    SH_1sigma = poisson_conf_interval(SH, interval='frequentist-confidence').T
    SHl = SH_1sigma[:, 0];
    SHu = SH_1sigma[:, 1]

    time=pars[:,1]
    # plt.scatter(HR,SH,label=gc)
    ##==plot HR-S+H diagram==##\
    print(len(HR),len(SH))
    fig, (ax1, ax2) = plt.subplots(ncols=1, nrows=2, sharex='all', gridspec_kw={'height_ratios': [1, 1]},
                                   figsize=(15, 10))
    HR=np.roll(HR,-8);SH=np.roll(SH,-8);HRl=np.roll(HRl,-8);HRu=np.roll(HRu,-8);SHl=np.roll(SHl,-8);SHu=np.roll(SHu,-8)

    ax1.errorbar(x=time*1/40.+1/80., y=HR, yerr=[HR - HRl, HRu - HR], xerr=1/80., fmt='.', capsize=3, elinewidth=1,color='blue')
    ax2.errorbar(x=time*1/40.+1/80., y=SH, yerr=[SH - SHl, SHu - SH], xerr=1/80., fmt='.', capsize=3, elinewidth=1,color='red')
    plt.show()
if __name__=='__main__':
    obsidlist = [14197., 14198., 13825., 13826., 13827., 13828., 14195., 14196.]
    path='/Users/baotong/Desktop/period_M31XRB/M31ACIS_txt/txt_all_obs_p90/'
    # a=input_parsfromtxt_lc(path,125,ebands=(500,2000),ebandh=(2000,8000))
    # for id in obsidlist:
    #     plot_125_HR_lc(int(id))
    run_one_BEHR(125)
    plot_125_HR_pfold()
# print(a)