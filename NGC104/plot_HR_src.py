#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.stats import poisson_conf_interval
from astropy.coordinates import SkyCoord
from astropy import units as u
import hawkeye as hawk
import rocket

ra_center = [6.0223292];
dec_center = [-72.0814444]
inter_radius = 3.17*60
path_in='/Users/baotong/Desktop/period_Tuc/'

def cal_HR_err(cts_s,cts_h):
    ## x and y must be integers
    u=(cts_h-cts_s)/(cts_s+cts_h)
    x=cts_s-cts_h;y=cts_s+cts_h
    if u<0:u=-u
    cts_s_1sigma = poisson_conf_interval(cts_s, interval='frequentist-confidence').T
    cts_h_1sigma = poisson_conf_interval(cts_h, interval='frequentist-confidence').T
    cts_s_low=cts_s_1sigma[0];cts_s_high=cts_s_1sigma[1];cts_h_low=cts_h_1sigma[0];cts_h_high=cts_h_1sigma[1]
    cts_s_err = min(cts_s_high-cts_s,cts_s-cts_s_low)
    cts_h_err = min(cts_h_high - cts_h, cts_h - cts_h_low)
    if cts_s==0: cts_s_err=0
    if cts_h==0: cts_h_err=0

    x_err=cts_s_err+cts_h_err
    y_err=x_err
    u_err=np.sqrt((x_err**2*y**2+y_err**2*x**2)/y**4)
    if np.isnan(u_err):u_err=0
    return (u,u_err)

def input_srcinfo():
    cat = fits.open(path_in + 'Cheng2019.fit')
    srcID_list=cat[1].data['Seq']
    ra = cat[1].data['RAJ2000']
    dec = cat[1].data['DEJ2000']
    L_05_8=cat[1].data['L0_5-8']
    L_05_2=cat[1].data['L0_5-2']
    L_2_8=cat[1].data['L2-8']
    netcts_t=cat[1].data['NetCts-t']
    netcts_s=cat[1].data['NetCts-s']
    netcts_h=cat[1].data['NetCts-h']
    info=np.column_stack((srcID_list,ra,dec,L_05_8,L_05_2,L_2_8,netcts_t,netcts_s,netcts_h))
    return info

def del_twoarray(arr1,arr2):
    for i in range(len(arr2)):
        arr1=np.delete(arr1,np.where(arr1==arr2[i]))
    return arr1

info=input_srcinfo()

inter_srcID = rocket.select_src_bypos(info[:,0], info[:,1], info[:,2], ra_c=ra_center, dec_c=dec_center,
                                      inter_radius=inter_radius)

inter_srcID=inter_srcID.astype('int')
netcts_t=info[:,6]
netcts_s=info[:,7]
netcts_h=info[:,8]
L_05_8=info[:,3]
inter_srcID=inter_srcID-1

def match_id(id_old_in):
    match_result=np.loadtxt('/Users/baotong/Desktop/period_Tuc/match_old_newfits.txt')
    id_new=match_result[:,0];id_old=match_result[:,1]
    list=[]
    for i in range(len(id_old_in)):
        list.append(int(id_new[np.where(id_old==id_old_in[i])][0]))
    return np.array(list)

def input_ID():
    inter_srcID_psrc=np.array([245,261,294,182])
    CV_srcID=np.array([206,321,453,402,462,364,304,350,345,283,258])
    # CV_srcID=np.array([343,321,291,279,245,230,232,418,409,274,384,294])
    AB_srcID=np.array([446,395,381,371,358,347,314,259,246,223,176,464,434,429,394,325,302,305,451,417,273,474,507,392,477,397])
    inter_srcID_psrc=match_id(inter_srcID_psrc)-1
    CV_srcID=match_id(CV_srcID)-1
    AB_srcID=match_id(AB_srcID)-1
    return (inter_srcID_psrc,CV_srcID,AB_srcID)

(inter_srcID_psrc,CV_srcID,AB_srcID)=input_ID()
inter_srcID_psrc=np.intersect1d(inter_srcID,inter_srcID_psrc)
CV_srcID=np.intersect1d(inter_srcID,CV_srcID)
AB_srcID=np.intersect1d(inter_srcID,AB_srcID)
HR=(netcts_h-netcts_s)/(netcts_h+netcts_s)
HR_err=[]
for i in range(len(HR)):
    HR_err.append(cal_HR_err(netcts_s[i],netcts_h[i])[1])
HR_err=np.array(HR_err)

path_out = '/Users/baotong/Desktop/aas/pXS_Tuc/figure/'
def plot_threefig(save=0,show=1):
    plt.figure(1,(9,6))
    plt.semilogy()
    plt.tick_params(labelsize=16)
    plt.ylabel(r'$\rm L_{x(0.5-8 keV)} (erg~s^{-1})$ ',hawk.font1)
    plt.xlabel('HR=(H-S)/(H+S)',hawk.font1)
    plt.scatter(HR[inter_srcID],L_05_8[inter_srcID],marker='.',color='grey',s=50)
    plt.scatter(HR[inter_srcID_psrc],L_05_8[inter_srcID_psrc],marker='*',color='green',s=100)
    plt.scatter(HR[CV_srcID],L_05_8[CV_srcID],marker='v',color='red',s=100)
    plt.scatter(HR[AB_srcID],L_05_8[AB_srcID],marker='o',edgecolors='purple',color='white',s=100)
    plt.legend(['All','LMXB','CV','AB'])
    plt.scatter(HR[197],L_05_8[197],marker='o',color='purple',s=100)
    plt.scatter(HR[282],L_05_8[282],marker='o',color='purple',s=100)
    plt.xlim(-1.6,1.0)
    plt.plot([-0.69,-0.69],[2e29,3e33],'--')
    for i in range(len(inter_srcID_psrc)):
        ##这里psrc只是LMXB
        plt.text(HR[inter_srcID_psrc[i]],L_05_8[inter_srcID_psrc[i]]*1.2,inter_srcID_psrc[i]+1,fontsize=12)
    for i in range(len(CV_srcID)):
        plt.text(HR[CV_srcID[i]],L_05_8[CV_srcID[i]]*1.2,CV_srcID[i]+1,fontsize=12)
    plt.text(HR[197], L_05_8[197] * 1.2, 197+ 1, fontsize=12)
    plt.text(HR[282], L_05_8[282] * 1.2, 282+ 1, fontsize=12)
    filt1 = np.intersect1d(inter_srcID,np.where(L_05_8<1e30))
    filt2 = np.intersect1d(inter_srcID, np.where((L_05_8> 1e30)&(L_05_8 < 1e31)))
    filt3 = np.intersect1d(inter_srcID, np.where((L_05_8> 1e31)&(L_05_8 < 1e32)))
    filt4 = np.intersect1d(inter_srcID, np.where((L_05_8> 1e32)&(L_05_8 < 1e33)))
    x11=HR[filt1];x11_err=HR_err[filt1]
    x12=HR[filt2];x12_err=HR_err[filt2]
    x13=HR[filt3];x13_err=HR_err[filt3]
    x14=HR[filt4];x14_err=HR_err[filt4]
    plt.errorbar(x=0.5,xerr=np.mean(x11_err),y=5.5e29,yerr=4.e29,fmt='.', capsize=5, elinewidth=1.5, ecolor='red',linewidth=1.0)
    plt.errorbar(x=0.5,xerr=np.mean(x12_err),y=5.5e30,yerr=4.e30,fmt='.', capsize=5, elinewidth=1.5, ecolor='green',linewidth=1.0)
    plt.errorbar(x=0.5,xerr=np.mean(x13_err),y=5.5e31,yerr=4.e31,fmt='.', capsize=5, elinewidth=1.5, ecolor='blue',linewidth=1.0)
    plt.errorbar(x=0.5,xerr=np.mean(x14_err),y=5.5e32,yerr=4.e32,fmt='.', capsize=5, elinewidth=1.5, ecolor='orange',linewidth=1.0)
    # for i in range(len(AB_srcID)):
    #     plt.text(HR[AB_srcID[i]], L_05_8[AB_srcID[i]]*1.2, AB_srcID[i] + 1,fontsize=12)
    if save:
        plt.savefig(path_out + '47Tuc_HR.eps', bbox_inches='tight', pad_inches=0.05)
    if show:
        plt.show()

    fig=plt.figure(2,figsize=(9,6))
    ax1=fig.add_subplot(211)
    bins=np.logspace(np.log10(1e29), np.log10(5e32), 19)
    # bins=np.delete(bins,-9)
    ## 19 or 30
    ax1.set_xscale('log')
    ax1.set_xlim(1e29,5e32)
    inter_srcID_filter=inter_srcID
    inter_srcID_filter=del_twoarray(inter_srcID,inter_srcID_psrc)
    # inter_srcID_filter = del_twoarray(inter_srcID_filter, AB_srcID)
    CV_candidate_index=np.intersect1d(inter_srcID_filter,np.where(HR>-0.69))

    const_frac=len(CV_candidate_index)/(len(inter_srcID_filter))

    allnum=ax1.hist(L_05_8[inter_srcID_filter],bins=bins,histtype='step',linewidth=2,color='grey')
    cvnum=ax1.hist(L_05_8[CV_candidate_index],bins=bins,histtype='step',linewidth=2,facecolor='r',
             hatch='/', edgecolor='k',fill=True)
    ax1.legend(['All','HR>-0.69'])
    ax1.plot([4e30,4e30],[0,80],'--',color='r')
    # plt.xlabel('Lx_0.5-8')
    ax1.set_ylabel('Number of sources',hawk.font1)
    ax1.tick_params(labelsize=16)

    ax2 = fig.add_subplot(212)
    ax2.tick_params(labelsize=16)
    ax2.set_xscale('log')
    ax2.set_xlabel(r'$\rm L_{x(0.5-8 keV)} (erg~s^{-1})$ ',hawk.font1)
    ax2.set_ylabel('CV Fraction',hawk.font1)
    ax2.plot(bins[:-1],np.zeros(len(bins[:-1]))+const_frac,'-.',color='c')
    ax2.plot([4e30,4e30],[0,1],'--',color='r')
    numofcv=cvnum[0];numofall=allnum[0]

    # erry=np.sqrt(numofcv)/(numofall+1e-2)
    ax2.step(bins,np.concatenate(([0],numofcv/(numofall+1e-2))),color='k')
    # ax2.errorbar(bins, np.concatenate(([0],numofcv/(numofall+1e-2))), yerr=np.concatenate(([0],erry)), fmt='.', capsize=1, elinewidth=1, ecolor='red')
    ax2.set_xlim(1e29,5e32)
    ax2.tick_params(labelsize=16)
    ax2.get_shared_x_axes().join(ax1, ax2)
    ra=info[:,1];dec=info[:,2]
    c1 = SkyCoord(ra * u.deg, dec * u.deg, frame='fk5')
    c2 = SkyCoord(ra_center * u.deg, dec_center * u.deg, frame='fk5')
    dist = c1.separation(c2)
    dist = dist.arcsec
    if save:
        plt.savefig(path_out + 'CV_2pop.eps', bbox_inches='tight', pad_inches=0.05)
    if show:
        plt.show()

    bright_CV_index=np.intersect1d(np.where(HR>-0.69),np.where(L_05_8>3.6e30))
    faint_CV_index=np.intersect1d(np.where(HR>-0.69),np.where(L_05_8<3.6e30))
    # bright_CV_index=np.where(L_05_8>1e31)
    # faint_CV_index=np.where(L_05_8<1e31)
    bright_CV_index=np.intersect1d(bright_CV_index,inter_srcID)
    faint_CV_index=np.intersect1d(faint_CV_index,inter_srcID)
    print('bright CV:',len(bright_CV_index))
    # print(netcts_t[251])
    print(len(np.where(netcts_t[faint_CV_index]>50)[0]))
    print('faint CV:',len(faint_CV_index))
    rh=3.17*60
    plt.figure(3,(9,6))
    plt.hist(np.log10(dist[bright_CV_index]/rh),bins=30,histtype='step',linewidth=2,cumulative=1,density=1,color='r')
    plt.hist(np.log10(dist[faint_CV_index]/rh),bins=30,histtype='step',linewidth=2,cumulative=1,density=1,color='blue')
    plt.legend(['bright CV','faint CV'],loc='upper left')
    plt.plot([np.log10(21.61/rh),np.log10(21.61/rh)],[0,0.9],'--',color='orange')
    plt.text(np.log10(21.61/rh)-0.18,0.92,'Core radius',fontdict=hawk.font2)
    plt.xlabel(r'$\rm log_{10}(R/R_h)$',hawk.font1)
    plt.ylabel('CDF',hawk.font1)
    plt.tick_params(labelsize=16)
    if save:
        plt.savefig(path_out + 'CDF_2pop.eps', bbox_inches='tight', pad_inches=0.05)
    if show:
        plt.show()

if __name__=="__main__":
    plot_threefig(save=1,show=1)