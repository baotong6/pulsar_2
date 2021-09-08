import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
import pandas as pd
import sys
import os
from tkinter import _flatten
import funcs_fits2txt as funcs
from scipy import interpolate
from astropy.coordinates import SkyCoord
from astropy import units as u

def read_erosita_cat(filename):
    cat=pd.read_excel(filename,header=0)
    srcid=cat['NAME']
    ra_hms=cat['RA']
    dec_hms=cat['DEC']
    ra=[];dec=[]
    for i in range(len(dec_hms)):
        skycoord=ra_hms[i]+dec_hms[i]
        c = SkyCoord(skycoord, unit=(u.hourangle, u.deg))
        ra.append(c.ra.value)
        dec.append(c.dec.value)

    return (ra,dec,srcid)

def make_epochfile(path,ID,outpath):
    TSTART=[];TSTOP=[]
    for obsid in ID:
        evtfile='pm00_{0}_020_EventList_c001_bary.fits'.format(obsid)
        TSTART.append(fits.open(path+evtfile)[1].header['TSTART'])
        TSTOP.append(fits.open(path+evtfile)[1].header['TSTOP'])
    exptime=np.array(TSTOP)-np.array(TSTART)
    res=np.column_stack((TSTART,TSTOP,exptime,ID))
    np.savetxt(outpath+'epoch_47Tuc.txt',res,fmt="%20.2f  %20.2f  %20.2f %10d")

def get_singlesrc_evt(evtfile,obsid,ra,dec,radius,outpath,outname):
    evtall=fits.open(evtfile)[1]
    time_all=evtall.data['TIME']
    energy_all=evtall.data['PI']
    RA_all=evtall.data['RA']
    DEC_all=evtall.data['DEC']
    reg=[ra,dec,radius]
    index_out=funcs.where_region(RA_all,DEC_all,reg)
    time=time_all[index_out];energy=energy_all[index_out]
    obsID=np.array([obsid for i in range(len(time))])
    [src_t, src_E, src_ID] = funcs.delete_photon_ID(time, energy, obsID,emin=500,emax=8000)

    src_t = src_t.astype('float')
    src_E = src_E.astype('float')
    src_ID = src_ID.astype('int')
    src_txt = np.column_stack((src_t, src_E, src_ID))
    src_txt = src_txt[src_txt[:, 0].argsort()]
    np.savetxt(outpath + outname +'_'+str(obsid)+ '.txt', src_txt, fmt="%.7f  %10.3f  %10d")

def merge_txt(ID,outpath,outname):
    res_t=[];res_E=[];res_ID=[]
    for obsid in ID:
        res_temp=np.loadtxt(outpath+outname+'_'+str(obsid)+'.txt')
        res_t.append(list(res_temp[:, 0]))
        res_E.append(list(res_temp[:, 1]))
        res_ID.append(list((res_temp[:, 2])))

    res_t = list(_flatten(res_t))
    res_E = list(_flatten(res_E))
    res_ID = list(_flatten(res_ID))
    result = np.column_stack((res_t, res_E, res_ID))
    result = result[result[:, 0].argsort()]

    np.savetxt(outpath + outname+ '_merge.txt', result, fmt="%.7f  %10.3f  %10d")

def get_psfradius(ra,dec,obsid,inpath,outpath,ecf=0.5,if_SN_radius=False,evtname=None):
    ## ra and dec should be array
    ##注意这里的数字位置都是从0开始计数##
    ##erosita现在的psfmap大小都是21x21，注意手动换算##
    ##(ra,dec) to (x,y)##
    image_file='{0}_05_5_img.fits'.format(obsid)
    image=fits.open(inpath+image_file)[0].data
    (xbin,ybin)=image.shape
    w = WCS(path + image_file)
    src_x, src_y = w.all_world2pix(ra, dec, 1)
    img_x=src_x/1024*21;img_y=src_y/1024*21
    psf_x=img_x.astype('int');psf_y=img_y.astype('int')

    psf_filename = 'psfmap_{0}.fits'.format(obsid)
    psffile = fits.open(inpath + psf_filename)
    psfmap_all = psffile[0].data
    ecflist=np.array([0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95])
    def output_psf_intervalue(ecf,psf_x,psf_y,img_x,img_y):
        indexpsf=np.where(ecflist==ecf)[0][0]
        psfmap_ecf=psfmap_all[indexpsf]
        psfmap_ecf=psfmap_ecf.T
        psfmap_ecf[np.where(psfmap_ecf==0.)]+=np.max(psfmap_ecf)
        psf_value=psfmap_ecf[psf_x,psf_y]  ##不太准，所以不用

        psf_bettervalue=[]
        f = interpolate.interp2d(np.arange(0.5,21.5,1),np.arange(0.5,21.5,1),psfmap_ecf,kind='cubic')
        for i in range(len(img_x)):
            psf_bettervalue.append(f(img_x[i],img_y[i])[0])
        psf_bettervalue=np.array(psf_bettervalue)
        return psf_bettervalue

    if if_SN_radius:
        evtname = 'pm00_{0}_020_EventList_c001_bary.fits'.format(obsid)
        for i in range(len(ra)):
            src_reg=[ra,dec,psf_value]

    return (psf_value,psf_bettervalue)


def make_region(ra,dec,radius,obsIDlist,outpath):
    radius=np.array(radius)*4/3600
    for obsid in obsIDlist:
        os.chdir(outpath)
        if os.path.exists('reg_{0}'.format(obsid)):
            os.chdir('reg_{0}/'.format(obsid))
        else:
            os.mkdir('reg_{0}/'.format(obsid))
            os.chdir('reg_{0}/'.format(obsid))
        with open('./all.reg', 'a+') as f2:
            f2.writelines('fk5' + '\n')
        for i in range(len(ra)):
            with open('./{0}.reg'.format(i + 1), 'w+') as f1:
                f1.writelines('fk5'+'\n')
                reg = 'circle(' + str(ra[i]) + ',' + str(dec[i]) + ',' + str(radius[i]) + ')'
                f1.writelines(reg)
            with open('./all.reg', 'a+') as f2:
                f2.writelines(reg + '\n')


def get_all_src_txt(catalog_file,obsIDlist,inpath,outpath):
    ###------read catalog file, need modification---##
    (ra_list,dec_list,id_list)=read_erosita_cat('/Users/baotong/Desktop/period_Tuc/erosita_cat_coord.xlsx')
    ###-------------------------------------------##
    if len(ra_list)!=len(id_list) or len(dec_list)!=len(id_list):
        print('ERROR!')
        return None
    for i in range(len(ra_list)):
        for j in range(obsid):
            evtfile=inpath+'pm00_{0}_020_EventList_c001_bary.fits'.format(obsIDlist[i])
            get_singlesrc_evt(evtfile=evtfile,obsid=obsIDlist[i],ra=ra_list[i],dec=dec_list[i],
                              radius=radius_list[i],outpath=outpath,outname=id_list[i])
        merge_txt(ID=obsIDlist,outpath=outpath,outname=id_list[i])

if __name__=='__main__':
    # path = '/Users/baotong/eSASS/data/47_Tuc/'
    path='/Users/baotong/eSASS/data/raw_data/47_Tuc/'
    ID=[700011,700013,700014,700163,700173,700174,700175]
    # make_epochfile(path=path,ID=ID,outpath=path+'txt/')

    # for obsid in ID:
    #     get_singlesrc_evt(path+'pm00_{0}_020_EventList_c001_bary.fits'.format(obsid),obsid=obsid,
    #                   ra=6.04485,dec=-72.07383,radius=20./3600,outpath='/Users/baotong/eSASS/data/47_Tuc/txt/',
    #                   outname='402')
    #
    # merge_txt(ID,outpath='/Users/baotong/eSASS/data/47_Tuc/txt/',outname='402')

    # get_psfradius(ra=np.array([6.0178,6.1078]),dec=np.array([-72.08281,-72.18281]),obsid=700011,inpath=path,outpath=path)
    (ra,dec,srcID)=read_erosita_cat('/Users/baotong/Desktop/period_Tuc/erosita_cat_coord.xlsx')
    (psf_value,psf_bettervalue)=get_psfradius(ra=ra, dec=dec, obsid=700011, inpath=path,
                  outpath=path,ecf=0.9)
    make_region(ra,dec,radius=psf_bettervalue,obsIDlist=ID,outpath=path)