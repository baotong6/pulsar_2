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


path="/Users/baotong/Desktop/period_Tuc/"
type=['47Tuc','terzan5','M28','omg_cen']


def add_pos_for_excel(label):
    if label=='omg_cen':dir=label[0:-4]
    else: dir=label
    res=pd.read_excel(path+'result_0.5_8_all.xlsx',label)
    #srclis=fits.open('/Users/baotong/Desktop/period_{0}/{1}_p50_i5_src_1_2_4_8.fits'.format(dir,label))
    srclis=fits.open('/Users/baotong/Desktop/period_Tuc//xray_properties-592.fits'.format(dir))
    #for 47Tuc
    ra=srclis[1].data['RAdeg']
    dec=srclis[1].data['DEdeg']
    seq=np.linspace(1,len(ra),len(ra))
    seq=seq.astype(np.int)
    result_out=np.column_stack((ra[res['seq']-1],dec[res['seq']-1]))
    print(result_out)
    np.savetxt(path+'radec_{0}.txt'.format(label),result_out,fmt="get_ipython().run_line_magic("10.5f", " %10.5f\")")


add_pos_for_excel(type[0])


def make_region(label):
    res=pd.read_excel(path+'result_0.5_8_all.xlsx',label)
    os.system('rm {0}'.format(path)+'all_pCV_{0}.reg'.format(label))
    with open(path+'all_pCV_{0}.reg'.format(label),'a+') as f2:
        f2.writelines('fk5'+'\n')
        
    ra=res['RA']
    dec=res['DEC']
    seq=res['seq']
    for i in range(len(ra)):
        with open(path+'all_pCV_{0}.reg'.format(label),'a+') as f2:
            reg = 'circle(' + str(ra[i]) + ',' + str(dec[i]) + ',' + str('2"') + ')'+" # color=green width=2 text={"+str(seq[i])+"}"
            f2.writelines(reg+'\n')


def make_CDFS_region():
    path='/Users/baotong/Desktop/CDFS/'
    srclist=fits.open(path+'7Ms_catalog.fit')
    ra=srclist[1].data['RAJ2000']
    dec=srclist[1].data['DEJ2000']
    seq=srclist[1].data['Seq']
    os.system('rm {0}'.format(path+'all_pCV_CDFS.reg'))
    with open(path+'all_pCV_CDFS.reg','a+') as f2:
        f2.writelines('fk5'+'\n')
    for i in range(len(ra)):
        with open(path+'all_pCV_CDFS.reg','a+') as f2:
            reg = 'circle(' + str(ra[i]) + ',' + str(dec[i]) + ',' + str('2"') + ')'+" # color=green width=2 text={"+str(seq[i])+"}"
            f2.writelines(reg+'\n')
        
make_CDFS_region()











make_region(type[2])


def print_MJD():
    path='/Volumes/pulsar/M28/merge_data/spectra/aprates/'
    EPOCH=np.loadtxt(path + 'M28_epoch.txt')
    obs_time=(EPOCH[:, 0] + EPOCH[:, 1]) / 2
    time = obs_time / 86400 + 2449352.5 - 2400000.5
    print(time)


print_MJD()


def plot_profile_distb():
    label_all=['47Tuc','terzan5','M28','omg_cen']
    pos_all=[[6.0236250,-72.0812833,3.17*60],  #47Tuc
             [267.0202083,-24.7790556,43.66],  #terzan5
             [276.1363750,-24.8702972,1.97*60], # M28
             [201.69700,-47.47947 , 5*60]]      #omega_cen
    distance_all=[0]
    period_all=[0]
    for i in range(4):
        label=label_all[i]
        pos=pos_all[i]
        res=pd.read_excel(path+'result_0.5_8_all.xlsx',label)
        ra=np.array(res['RA'])
        dec=np.array(res['DEC'])
        seq=np.array(res['seq'])
        period=np.array(res['P_out'])
        i=0
        while i<len(seq):
            if period[i]<4500:
            #if str(seq[i])[-3:]=='001':
                ra=np.delete(ra,i)
                dec=np.delete(dec, i)
                seq=np.delete(seq, i)
                period=np.delete(period, i)
            else:i+=1
        distance=((ra-pos[0])**2+(dec-pos[1])**2)**0.5*3600/pos[-1]
        distance_all=np.concatenate((distance_all,distance))
        period_all=np.concatenate((period_all,period))
    distance_all=distance_all[1:]
    period_all=period_all[1:]
    # P_gap = [7740.0, 11448.0]
    # plt.plot([0,10],[P_gap[0], P_gap[0]], '-', lw=2., color='grey')
    # plt.plot([0,10],[P_gap[1], P_gap[1]], '-', lw=2., color='grey')
    #print(distance)
    period_plot=[]
    num_density=[]
    num_periodbin=[]
    #plt.semilogy()
    #plt.ylim(1e-4,1e2)
    bins=np.linspace(0,10,9)
    bins=np.linspace(0,12,11)
    for i in range(len(bins)-1):
        period_bin=period_all[np.where((distance_all<bins[i+1])&(distance_all>=bins[i]))]
        num_density.append(len(period_bin)/(np.pi*(bins[i+1]**2-bins[i]**2)))
        #print(len(period_bin),np.pi*(bins[i+1]**2-bins[i]**2))
        period_acm=np.mean(period_bin)
        period_plot.append(period_acm)
        num_periodbin.append(len(period_bin))
    x=bins[0:-1]+(bins[1]-bins[0])/2.
    y=num_density
    num_periodbin=np.array(num_periodbin)
    func = interpolate.interp1d(x, y, kind='cubic')
    x_new = np.linspace(min(x), max(x), 50)
    y_new = func(x_new)
    #plt.plot(x, y, color='red')
    #plt.plot(x_new, y_new,color='green',linestyle='--')
    #plt.scatter(x,period_plot,s=num_periodbin*20,marker='o',edgecolors='r')  # period_mean
    plt.scatter(distance_all,period_all)  # all period
    #plt.hist(period_all/distance_all,bins=50)
    #plt.plot([pos[-1], pos[-1]], [1000,50000], '--')

    plt.show()



plot_profile_distb()


def unit_convert(ra_dec):
    test = SkyCoord(ra_dec, unit = (u.deg,u.deg), frame='icrs', distance = 0.48837*u.AU)
    print(test.to_string('hmsdms'))
    


unit_convert('6.0178 -72.08281')


# help(SkyCoord)


def make_region_GCCR():
    path='/Users/baotong/Desktop/period/'
    filename='/Users/baotong/Desktop/paper/NSC/GCCR_tab.txt'
    info_GCCR=[]
    with open(filename,'r') as file_to_read:
        while True:
            lines = file_to_read.readline() # 整行读取数据
            info_GCCR.append(lines)
            if not lines:
                break
                pass
            
    XSRC=np.array([1 ,4 ,14 ,15 ,31 ,32 ,35 ,36 ,37 ,38 ,41 ,47 ,49 ,50 ,56 ,57 ,
          58 ,59 ,61 ,64 ,67 ,70 ,72 ,73 ,74 ,77 ,82 ,85 ,87 ,89 ,92 ,95,
          96 ,97 ,98 ,99 ,100 ,101 ,103 ,106 ,107 ,110])
    
    info_GCCR=info_GCCR[0:-1]  ##去掉末尾空行
    label=[];ra=[];dec=[];
    for i in range(len(info_GCCR)):  
        label_i,ra_i,dec_i=[str(i) for i in info_GCCR[i][0:-1].split(';')]   ##去掉末尾换行符
        print(dec_i)
        label.append(label_i);ra.append(ra_i);dec.append(dec_i)
    with open(path+'all_GCCR.reg','a+') as f2:
        f2.writelines('fk5'+'\n')
    for i in range(len(ra)):
        with open(path+'all_GCCR.reg','a+') as f2:
            if (i+1) in XSRC:
                reg = 'circle(' + str(ra[i]) + ',' + str(dec[i]) + ',' + str('3"') + ')'+" # color=green width=2 text={"+str(label[i][4:])+"}"
                f2.writelines(reg+'\n')
            else:
                reg = 'circle(' + str(ra[i]) + ',' + str(dec[i]) + ',' + str('3"') + ')'+" # color=green width=1 text={"+str(label[i][-4:])+"}"
                f2.writelines(reg+'\n')
    


# make_region_GCCR()



