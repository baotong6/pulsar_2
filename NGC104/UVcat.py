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
from astropy.coordinates import SkyCoord,match_coordinates_sky
path='/Users/baotong/Desktop/period_Tuc/'
font1 = {'family': 'Normal',
         'weight': 'normal',
         'size': 18, }
def var_name(var,all_var=locals()):
    return [var_name for var_name in all_var if all_var[var_name] is var][0]

def trans_coord():
    ##useless
    path='/Users/baotong/Desktop/period_Tuc/ngc0104/'
    filename='hlsp_hugs_hst_wfc3-uvis_ngc0104_f336w_v1_stack-0225s.fits'
    img_file=fits.open(path+filename)
    img_data=img_file[0].data
    w = WCS(path + filename)
    # src_x, src_y = w.all_world2pix(ra, dec, 1)
    ra,dec=w.all_pix2world(5000,5000,1)

def match_two_cat(ra_xray,dec_xray,ra_HST,dec_HST,rad=0.5):
    c_chandra= SkyCoord(ra=ra_xray*u.degree, dec=dec_xray*u.degree)
    c_HST=SkyCoord(ra=ra_HST*u.degree, dec=dec_HST*u.degree)

    idx, d2d, d3d=match_coordinates_sky(matchcoord=c_HST, catalogcoord=c_chandra, nthneighbor=1)
    max_sep = rad * u.arcsec
    sep_constraint = d2d < max_sep
    c_matches = c_HST[sep_constraint]
    catalog_matches = c_chandra[idx[sep_constraint]]
    d2d_matches=d2d[sep_constraint].to(u.arcsec)

    return (c_matches,catalog_matches,d2d_matches,sep_constraint)

def read_uv_cat(filename):
    cat_dat=np.loadtxt(filename)
    ra=cat_dat[:,33];dec=cat_dat[:,34]
    F275W=cat_dat[:,2];F336W=cat_dat[:,8]
    F435W=cat_dat[:,14];F606W=cat_dat[:,20]
    F814W=cat_dat[:,26]
    return (ra,dec,F275W,F336W,F435W,F606W,F814W)

def read_chandra_cat(filename):
    cat_file=fits.open(filename)
    cat_data=cat_file[1].data
    ra=cat_data['RAdeg']
    dec=cat_data['DEdeg']
    pos_err=cat_data['PosErr']

    return (ra,dec,pos_err)

def read_erosita_cat(filename):
    cat=pd.read_excel(filename,header=0)
    ra_hms=cat['RA']
    dec_hms=cat['DEC']
    ra=[];dec=[]
    for i in range(len(dec_hms)):
        skycoord=ra_hms[i]+dec_hms[i]
        c = SkyCoord(skycoord, unit=(u.hourangle, u.deg))
        ra.append(c.ra.value)
        dec.append(c.dec.value)

    return (ra,dec)
def plot_three_cat():
    (ra1,dec1,pos_err)=read_chandra_cat(path+'xray_properties-592.fits')
    (ra2,dec2)=read_erosita_cat(path+'erosita_cat_coord.xlsx')
    (ra3,dec3,F275W,F336W,F435W,F606W,F814W)=read_uv_cat(path+'ngc0104/ngc104_meth1.txt')
    plt.xlabel('RA (deg)')
    plt.ylabel('DEC (deg)')

    plt.scatter(ra1, dec1, marker='o',color='red')
    plt.scatter(ra2, dec2, marker='.', color='green')
    plt.scatter(ra3, dec3, marker='.', color='blue')
    plt.legend(['Chandra','eROSITA','HST_UV'])
    plt.show()
def plot_HR_cat():
    (ra3,dec3,F275W,F336W,F435W,F606W,F814W)=read_uv_cat(path+'ngc0104/ngc104_meth1.txt')
    index=np.where((F606W>-99)&(F814W>-99))
    print(len(ra3),len(index[0]))
    print(np.sort(F606W))
    ra3=ra3[index];dec3=dec3[index];F275W=F275W[index];F336W=F336W[index]
    color=F275W-F336W
    plt.scatter(color,F275W,marker='.', s=0.1,color='green')
    plt.gca().invert_yaxis()
    plt.show()

def plot_close_HR(xrayid,band1,band2,band3,info,save=0,show=1):
    # color=band1-band2
    # plot x=color,y=band3
    # (ra3, dec3, F275W, F336W, F435W, F606W, F814W) = read_uv_cat(path + 'ngc0104/ngc104_meth1.txt')
    (ra1,dec1,pos_err)=read_chandra_cat(path+'xray_properties-592.fits')
    # index=np.where((band1>-99)&(band2>-99)&(band3>-99))
    # ra3=ra3[index];dec3=dec3[index];band1=band1[index];band2=band2[index];band3=band3[index]
    color=info[band1]-info[band2]
    ra3=info['ra3'];dec3=info['dec3']
    id3=np.arange(1,len(ra3)+1,1)
    (c_matches, catalog_matches, d2d_matches,sep_constraint)=match_two_cat([ra1[xrayid-1]],[dec1[xrayid-1]],ra3,dec3,rad=pos_err[xrayid-1]+0.1)
    print(pos_err[xrayid-1])
    print(c_matches,catalog_matches,d2d_matches,id3[sep_constraint])
    sgbindex=np.where((color[sep_constraint]>1.0)&(info[band3][sep_constraint]<16.5))[0]
    sgbid3=id3[sep_constraint][sgbindex]
    close_dist=d2d_matches[sgbindex]
    plt.figure(1,(6,10))
    plt.title(f'{xrayid},rad={pos_err[xrayid-1]} arcsec')
    plt.xlim(-2,8)
    plt.ylim(10,24)
    plt.plot([1.0,1.0],[10,23],'--')
    plt.plot([-1,8],[16.5,16.5], '--')
    plt.scatter(color, info[band3], marker='.', s=0.1, color='green')
    plt.scatter(color[sep_constraint],info[band3][sep_constraint],marker='*',s=200,color='red')
    plt.gca().invert_yaxis()
    plt.xlabel(f'{band1} - {band2}',font1)
    plt.ylabel(f'{band3}',font1)
    plt.tick_params(labelsize=18)
    if save:
        plt.savefig('/Users/baotong/Desktop/period_Tuc/ngc0104/figure/' + f'match_IR_{xrayid}.pdf',bbox_inches='tight', pad_inches=0.01)
    if show:plt.show()
    else:plt.close()

    return (id3[sep_constraint], d2d_matches)
    # return (sgbid3,close_dist)

def plot_pure_HR(sgbid3,band1,band2,band3,info,save=0,show=1,outname=None):
    color = info[band1] - info[band2]
    ra3 = info['ra3'];
    dec3 = info['dec3']
    id3 = np.arange(1, len(ra3) + 1, 1)
    plt.figure(1,(6,10))
    plt.xlim(-4,5)
    plt.ylim(14,28)
    plt.plot([1.0,1.0],[10,23],'--')
    plt.plot([-1,8],[18,18], '--')
    plt.scatter(color, info[band3], marker='.', s=0.1, color='green')
    plt.scatter(color[sgbid3-1],info[band3][sgbid3-1],marker='*',s=200,color='red')
    plt.gca().invert_yaxis()
    plt.xlabel(f'{band1} - {band2}',font1)
    plt.ylabel(f'{band3}',font1)
    plt.tick_params(labelsize=18)
    if save:
        plt.savefig('/Users/baotong/Desktop/period_Tuc/ngc0104/figure/' + f'match_UV_{outname}.pdf',bbox_inches='tight', pad_inches=0.01)
    if show:plt.show()
    else:plt.close()
def astrometry(xrayidlist):
    ra_xray=[];dec_xray=[];
    ra_HST=[];dec_HST=[]
    for xrayid in xrayidlist:
        (ra3,dec3,F275W,F336W,F435W,F606W,F814W)=read_uv_cat(path+'ngc0104/ngc104_meth1.txt')
        (ra1,dec1,pos_err)=read_chandra_cat(path+'xray_properties-592.fits')

        index=np.where((F275W>-99)&(F435W>-99))
        ra3=ra3[index];dec3=dec3[index];F275W=F275W[index];F435W=F435W[index]
        color=F275W-F435W
        id3=np.arange(1,len(ra3)+1,1)
        (c_matches, catalog_matches, d2d_matches,sep_constraint)=match_two_cat([ra1[xrayid-1]],[dec1[xrayid-1]],ra3,dec3,rad=0.1)
        ra_HST.append(c_matches[0].ra.value)
        dec_HST.append(c_matches[0].dec.value)
        ra_xray.append(catalog_matches[0].ra.value)
        dec_xray.append(catalog_matches[0].dec.value)
        print(c_matches,catalog_matches,d2d_matches,id3[sep_constraint])
    ra_xray=np.array(ra_xray);ra_HST=np.array(ra_HST);dec_xray=np.array(dec_xray);dec_HST=np.array(dec_HST)
    d_RA=ra_HST-ra_xray;d_dec=dec_HST-dec_xray
    return (ra_xray,ra_HST,dec_xray,dec_HST,d_RA,d_dec)
if __name__=='__main__':
    (ra3, dec3, F275W, F336W, F435W, F606W, F814W) = read_uv_cat(path + 'ngc0104/ngc104_meth1.txt')
    info={'ra3':ra3,'dec3':dec3,'F275W':F275W,'F336W':F336W,'F435W':F435W,'F606W':F606W,'F814W':F814W}
    xrayid_list=[453,206,402,462,364,304,350,283,321,345,258]
    for xrayid in xrayid_list:
        (sgbid3,close_dist)=plot_close_HR(xrayid=xrayid,band1='F435W',band2='F814W',band3='F814W',info=info,save=1,show=0)
        print(sgbid3,close_dist)
        plot_pure_HR(sgbid3=sgbid3, band1='F275W',band2='F336W',band3='F275W',info=info,save=1,show=0,outname=xrayid)
    # (ra_xray,ra_HST,dec_xray,dec_HST,d_RA,d_dec)=astrometry(xrayidlist=[453,206,402,304,321,258])
    # print(d_RA*3600,d_dec*3600)