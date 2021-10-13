import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
import pandas as pd
import sys
import os
from tkinter import _flatten
from scipy import interpolate
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.stats import poisson_conf_interval
import vg
import funcs_timing as funcs

ra_center=6.022318
dec_center=-72.081443

def angle(v1, v2):
  dx1 = v1[2] - v1[0]
  dy1 = v1[3] - v1[1]
  dx2 = v2[2] - v2[0]
  dy2 = v2[3] - v2[1]
  angle1 = math.atan2(dy1, dx1)
  angle1 = int(angle1 * 180/math.pi)
  # print(angle1)
  angle2 = math.atan2(dy2, dx2)
  angle2 = int(angle2 * 180/math.pi)
  # print(angle2)
  if angle1*angle2 >= 0:
    included_angle = abs(angle1-angle2)
  else:
    included_angle = abs(angle1) + abs(angle2)
    if included_angle > 180:
      included_angle = 360 - included_angle
  return included_angle

def read_erosita_cat(filename):
    cat=pd.read_excel(filename,header=0)
    srcid=cat['NAME']
    srcid=np.array(srcid)
    ra_hms=cat['RA']
    dec_hms=cat['DEC']
    ra=[];dec=[]
    for i in range(len(dec_hms)):
        skycoord=ra_hms[i]+dec_hms[i]
        c = SkyCoord(skycoord, unit=(u.hourangle, u.deg))
        ra.append(c.ra.value)
        dec.append(c.dec.value)

    return (ra,dec,srcid)

def dist_center_twocat(catname1,catname2):
    (ra_eR,dec_eR,srcIDlist)=read_erosita_cat(catname1)
    cat2=fits.open(catname2)
    ra_chand=cat2[1].data['RAdeg']
    dec_chand=cat2[1].data['DEdeg']
    srcIDlist_chand=np.arange(1,len(ra_chand)+1,1)

    c1=SkyCoord(ra_eR*u.deg,dec_eR*u.deg,frame='fk5')
    c2=SkyCoord(ra_chand*u.deg,dec_chand*u.deg,frame='fk5')
    cc=SkyCoord(ra_center*u.deg,dec_center*u.deg,frame='fk5')
    dist1=c1.separation(cc)
    dist2=c2.separation(cc)
    return (dist1,dist2)

def get_err_num_poisson(num1,num2):
    num1_err=np.array(poisson_conf_interval(num1,interval='frequentist-confidence'))
    num1_err[0]=num1-num1_err[0]
    num1_err[1]=num1_err[1]-num1

    num2_err=np.array(poisson_conf_interval(num2,interval='frequentist-confidence'))
    num2_err[0]=num2-num2_err[0]
    num2_err[1]=num2_err[1]-num2

    return (num1_err,num2_err)

def plot_radial_2cat():
    path='/Users/baotong/Desktop/period_Tuc/'
    catname1=path+'erosita_cat_coord.xlsx'
    catname2=path+'xray_properties-592.fits'
    (dist1,dist2)=dist_center_twocat(catname1,catname2)
    dist1=dist1.arcmin;dist2=dist2.arcmin
    bins1=np.arange(int(np.min(dist1)),int(np.max(dist1))+2,1)
    bins2=np.arange(int(np.min(dist2)),int(np.max(dist2))+2,0.5)

    num1=plt.hist(dist1,bins=bins1,histtype='step')[0]
    num2=plt.hist(dist2,bins=bins2,histtype='step')[0]
    plt.close()

    (num1_err,num2_err)=get_err_num_poisson(num1,num2)

    for i in range(len(num1)):
        num1[i]/=np.pi*(bins1[i+1]**2-bins1[i]**2)
        num1_err[0][i] /= np.pi * (bins1[i + 1] ** 2 - bins1[i] ** 2)
        num1_err[1][i] /= np.pi * (bins1[i + 1] ** 2 - bins1[i] ** 2)

    for i in range(len(num2)):
        num2[i]/=np.pi*(bins2[i+1]**2-bins2[i]**2)
        num2_err[0][i] /= np.pi * (bins2[i + 1] ** 2 - bins2[i] ** 2)
        num2_err[1][i] /= np.pi * (bins2[i + 1] ** 2 - bins2[i] ** 2)
    plt.errorbar(x=bins1[:-1]+0.5,y=num1,yerr=num1_err,xerr=np.zeros(len(num1))+0.5,marker='.',linestyle='')
    plt.errorbar(x=bins2[:-1]+0.25,y=num2,yerr=num2_err,xerr=np.zeros(len(num2))+0.25,marker='.',linestyle='')
    plt.xlabel('R (arcmin)',funcs.font2)
    plt.ylabel('Number of source per arcmin^2',funcs.font2)
    plt.tick_params(labelsize=16)
    plt.semilogy()
    plt.show()

    return None

def plot_azimuth_SMC():
    path='/Users/baotong/Desktop/period_Tuc/'
    catname1=path+'erosita_cat_coord.xlsx'
    catname2=path+'xray_properties-592.fits'
    (ra_eR,dec_eR,srcIDlist)=read_erosita_cat(catname1)
    cat2=fits.open(catname2)
    ra_chand=cat2[1].data['RAdeg']
    dec_chand=cat2[1].data['DEdeg']
    srcIDlist_chand=np.arange(1,len(ra_chand)+1,1)
    (dist1,dist2)=dist_center_twocat(catname1,catname2)
    dist1=dist1.arcmin;dist2=dist2.arcmin
    # c3=SkyCoord(ra=13.1583*u.degree,dec=-72.8003*u.degree,distance=62.44*u.kpc)  ##smc coord
    # c2=SkyCoord(ra=ra_center*u.degree,dec=dec_center*u.degree,distance=4.0*u.kpc)
    # c1=SkyCoord(ra=ra_eR*u.degree,dec=dec_eR*u.degree,distance=4.0*u.kpc)
    # c0=SkyCoord(ra=ra_chand*u.degree,dec=dec_chand*u.degree,distance=4.0*u.kpc)
    # dist_chand_smc = c0.separation(c3)
    # dist_eR_smc=c1.separation(c3)
    # dist_ngc104_smc=c2.separation(c3)
    c_center = SkyCoord(ra=ra_center * u.degree, dec=dec_center * u.degree)
    c_smc=SkyCoord(ra=13.1583*u.degree,dec=-72.8003*u.degree)
    c_eR=SkyCoord(ra=ra_eR*u.degree,dec=dec_eR*u.degree)
    c_chand=SkyCoord(ra=ra_chand*u.degree,dec=dec_chand*u.degree)

    v_center=np.array([c_center.cartesian.x.value,c_center.cartesian.y.value,c_center.cartesian.z.value])
    v_smc = np.array([c_smc.cartesian.x.value, c_smc.cartesian.y.value, c_smc.cartesian.z.value])
    v_eR = np.array([c_eR.cartesian.x.value, c_eR.cartesian.y.value, c_eR.cartesian.z.value]).T
    v_chand = np.array([c_chand.cartesian.x.value, c_chand.cartesian.y.value, c_chand.cartesian.z.value]).T
    v1=v_eR-v_center
    v2=v_smc-v_center
    v3=v_chand-v_center
    include_ang=vg.angle(v1,v2)
    include_ang_chandra=vg.angle(v3,v2)
    # include_ang-=90

    index_eR_1 = np.where((dist1 < 20) & (dist1 > 4))[0]
    index_eR_2 = np.where((dist1 < 40) & (dist1 > 20))[0]
    width=10
    bins_fi=np.arange(0,180+width,width)

    num_in=plt.hist(include_ang[index_eR_1],bins=bins_fi,histtype='step')[0]
    num_out=plt.hist(include_ang[index_eR_2], bins=bins_fi, histtype='step')[0]
    num_chandra=plt.hist(include_ang_chandra, bins=bins_fi, histtype='step')[0]

    plt.close()
    (num_in_err, num_out_err) = get_err_num_poisson(num_in, num_out)
    num_chandra_err=np.array(poisson_conf_interval(num_chandra,interval='frequentist-confidence'))
    num_chandra_err[0]=num_chandra-num_chandra_err[0]
    num_chandra_err[1]=num_chandra_err[1]-num_chandra

    for i in range(len(num_in)):
        num_in[i]/=np.pi*(20**2-4**2)*width/180
        num_in_err[0][i] /= np.pi*(20**2-4**2)*width/180
        num_in_err[1][i] /= np.pi*(20**2-4**2)*width/180

    for i in range(len(num_out)):
        num_out[i]/=np.pi*(40**2-20**2)*width/180
        num_out_err[0][i] /= np.pi*(40**2-20**2)*width/180
        num_out_err[1][i] /= np.pi*(40**2-20**2)*width/180



    plt.errorbar(x=bins_fi[:-1]+width/2,y=num_in,yerr=num_in_err,xerr=np.zeros(len(num_in))+width/2,marker='.',linestyle='')
    # plt.errorbar(x=bins_fi[:-1]+width/2,y=num_out,yerr=num_out_err,xerr=np.zeros(len(num_out))+width/2,marker='.',linestyle='')
    # plt.errorbar(x=bins_fi[:-1] + width / 2, y=num_chandra, yerr=num_chandra_err, xerr=np.zeros(len(num_chandra)) + width / 2,
    #              marker='.', linestyle='')

    plt.xlabel('fi (degree)',funcs.font2)
    plt.ylabel('Number of source per arcmin^2',funcs.font2)
    plt.tick_params(labelsize=16)
    plt.semilogy()
    # plt.legend(['in(0-20 arcmin)','out(20-40 arcmin)'])
    plt.show()


    # plt.plot([dist_ngc104_smc.arcmin,dist_ngc104_smc.arcmin],[0,30],'--',color='red')
    # plt.show()


if __name__=='__main__':
    # plot_radial_2cat()
    path='/Users/baotong/Desktop/period_Tuc/'
    plot_azimuth_SMC()