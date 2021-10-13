import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
import pandas as pd
import sys
import os
from tkinter import _flatten
import funcs_fits2txt as funcs
from funcs_fits2txt import Circle
from scipy import interpolate
from astropy.coordinates import SkyCoord
from astropy import units as u
from bisect import bisect_left,bisect_right

def make_epochfile(path,ID,outpath):
    TSTART=[];TSTOP=[]
    for obsid in ID:
        evtfile='pm00_{0}_020_EventList_c001_bary.fits'.format(obsid)
        TSTART.append(fits.open(path+evtfile)[1].header['TSTART'])
        TSTOP.append(fits.open(path+evtfile)[1].header['TSTOP'])
    exptime=np.array(TSTOP)-np.array(TSTART)
    res=np.column_stack((TSTART,TSTOP,ID,exptime))
    np.savetxt(outpath+'epoch_47Tuc.txt',res,fmt="%20.2f  %20.2f  %20.2f %10d")

def get_singlesrc_evt(evtfile,imgname,obsid,ra,dec,emin,emax,radius,outpath,outname):
    (phy_x,phy_y)=funcs.trans_radec2xy(imgname,ra,dec)
    reg=[phy_x,phy_y,radius]
    [src_t, src_E, src_ID]=funcs.get_evt_srcreg(evtname=evtfile,obsid=obsid,src_reg=reg,emin=emin,emax=emax)
    obsID=np.array([obsid for i in range(len(src_t))])

    src_t = src_t.astype('float')
    src_E = src_E.astype('float')
    src_ID = src_ID.astype('int')
    src_txt = np.column_stack((src_t, src_E, src_ID))
    src_txt = src_txt[src_txt[:, 0].argsort()]
    np.savetxt(outpath + outname +'_'+str(obsid)+ '.txt', src_txt, fmt="%.7f  %10.3f  %10d")

def get_singlebkg_evt(evtfile,imgname,obsid,ra,dec,emin,emax,radius_src,radius_list,ra_list,dec_list,outpath,outname):
    ## 这个radius和radius_list都应该用psf90
    (phy_x,phy_y)=funcs.trans_radec2xy(imgname,ra,dec)
    bkg_reg=[phy_x,phy_y,radius_src,2*radius_src]

    [bkg_t, bkg_E, bkg_ID,bkg_area]=funcs.get_evt_bkgreg(evtname=evtfile, imgname=imgname, obsid=obsid,
                   bkg_reg=bkg_reg, ra_list=ra_list, dec_list=dec_list, radius_list=radius_list,
                   emin=emin, emax=emax)
    bkg_t = bkg_t.astype('float')
    bkg_E = bkg_E.astype('float')
    bkg_ID = bkg_ID.astype('int')
    bkg_txt = np.column_stack((bkg_t, bkg_E, bkg_ID))
    bkg_txt = bkg_txt[bkg_txt[:, 0].argsort()]
    np.savetxt(outpath + outname +'_bkg_'+str(obsid)+ '.txt', bkg_txt, fmt="%.7f  %10.3f  %10d")
    return bkg_area

def merge_txt(obsIDlist,inpath,outpath,outname,epoch_file):
    res_t=[];res_E=[];res_ID=[]
    epoch_ID=[];epoch_start=[];epoch_stop=[];epoch_expt=[]

    epoch_all = np.loadtxt(inpath + epoch_file)
    if epoch_all.ndim == 1:
        epoch_all = np.array([epoch_all])
    obs_tstart=epoch_all[:,0]
    obs_tstop = epoch_all[:,1]
    obs_ID_all=epoch_all[:,2]
    obs_expt=epoch_all[:,3]
    obs_ID_all = obs_ID_all.astype(int)

    for i in range(len(obs_ID_all)):
        obsid=obs_ID_all[i]
        # txt_filename=inpath+'txt_psf50_{0}/'.format(obsid)+outname+'_'+str(obsid)+'.txt'
        txt_filename=inpath+'txt_psf50_{0}/'.format(obsid)+outname+'_bkg_'+str(obsid)+'.txt'
        if os.path.exists(txt_filename):
            res_temp=np.loadtxt(txt_filename)
            if res_temp.size==0:
                print('Empty')
                continue
            if res_temp.ndim == 1:
                res_temp=np.array([res_temp])
            if res_temp.ndim == 2:
                res_t.append(list(res_temp[:, 0]))
                res_E.append(list(res_temp[:, 1]))
                res_ID.append(list((res_temp[:, 2])))

                epoch_ID.append(obs_ID_all[i])
                epoch_start.append(obs_tstart[i])
                epoch_stop.append(obs_tstop[i])
                epoch_expt.append(obs_expt[i])

    res_t = list(_flatten(res_t))
    res_E = list(_flatten(res_E))
    res_ID = list(_flatten(res_ID))
    result = np.column_stack((res_t, res_E, res_ID))
    result = result[result[:, 0].argsort()]

    epoch_info = np.column_stack((epoch_start, epoch_stop, epoch_ID, epoch_expt))
    # np.savetxt(outpath + 'epoch_src_' + str(outname) + '.txt', epoch_info,
    #            fmt='%15.2f %15.2f %10d %20.2f')
    np.savetxt(outpath + 'epoch_bkg_' + str(outname) + '.txt', epoch_info,
               fmt='%15.2f %15.2f %10d %20.2f')
    # np.savetxt(outpath + outname+ '.txt', result, fmt="%.7f  %10.3f  %10d")
    np.savetxt(outpath + outname+ '_bkg.txt', result, fmt="%.7f  %10.3f  %10d")

def get_psfradius(srcID,ra,dec,obsid,inpath,outpath,ecf=0.5,if_SN_radius=False):
    ## ra and dec should be array
    ##注意这里的数字位置都是从0开始计数##
    ##erosita现在的psfmap大小都是21x21，注意手动换算##
    ##返回值均为physical坐标##
    ##(ra,dec) to (x,y)##
    image_file='{0}_05_5_img.fits'.format(obsid)
    image=fits.open(inpath+image_file)[0].data
    (xbin,ybin)=image.shape
    w = WCS(path + image_file)
    src_x, src_y = w.all_world2pix(ra, dec, 1)
    # print(src_x,src_y)
    img_x=src_x/xbin*21;img_y=src_y/ybin*21

    ##先筛掉不在图像中的
    expmap_file='expmap_05_5_{0}.fits'.format(obsid)
    expmap=fits.open(inpath+expmap_file)[0].data
    not_include_index=np.zeros(len(srcID))+1
    for i in range(len(ra)):
        if src_x[i]<=0 or src_x[i]>=xbin-1 or src_y[i]<=0 or src_y[i]>=xbin-1:
            not_include_index[i]=0
            src_x[i]=0;img_x[i]=0;src_y[i]=0;img_y[i]=0
        elif expmap[int(src_x[i]),int(src_y[i])]==0:
            not_include_index[i]=0

    psf_x=img_x.astype('int');psf_y=img_y.astype('int')

    psf_filename = 'psfmap_{0}.fits'.format(obsid)
    psffile = fits.open(inpath + psf_filename)
    psfmap_all = psffile[0].data

    ecflist=np.array([0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95])

    def output_psf_intervalue(ecf,psf_x,psf_y,img_x,img_y,if_ecflist=False):
        if if_ecflist==False:
            indexpsf=np.where(ecflist==ecf)[0][0]
            psfmap_ecf = psfmap_all[indexpsf]
            psfmap_ecf = psfmap_ecf.T
            psfmap_ecf[np.where(psfmap_ecf == 0.)] += np.max(psfmap_ecf)
            psf_value = psfmap_ecf[psf_x, psf_y]  ##不太准，所以不用

            psf_bettervalue = []
            f = interpolate.interp2d(np.arange(0.5, 21.5, 1), np.arange(0.5, 21.5, 1), psfmap_ecf, kind='cubic')
            for i in range(len(img_x)):
                psf_bettervalue.append(f(img_x[i], img_y[i])[0])
            psf_bettervalue = np.array(psf_bettervalue)*80

            return psf_bettervalue
        else:
            psf_bettervalue_list=[]
            for i in range(len(ecflist)):
                psfmap_ecf=psfmap_all[i]
                psfmap_ecf=psfmap_ecf.T
                psfmap_ecf[np.where(psfmap_ecf==0.)]+=np.max(psfmap_ecf)
                psf_value=psfmap_ecf[psf_x,psf_y]  ##不太准，所以不用

                psf_bettervalue=[]
                f = interpolate.interp2d(np.arange(0.5,21.5,1),np.arange(0.5,21.5,1),psfmap_ecf,kind='cubic')
                for i in range(len(img_x)):
                    psf_bettervalue.append(f(img_x[i],img_y[i])[0])
                # psf_bettervalue=np.array(psf_bettervalue)
                psf_bettervalue_list.append(psf_bettervalue)
            psf_bettervalue_list=np.array(psf_bettervalue_list)*80

            return psf_bettervalue_list
    if not if_SN_radius:
        psf_bettervalue = output_psf_intervalue(ecf, psf_x, psf_y, img_x, img_y, if_ecflist=False)
        return (psf_bettervalue,not_include_index)

    if if_SN_radius:
        ##此时输入参数ecf是无用的
        evtname = 'pm00_{0}_020_EventList_c001_bary.fits'.format(obsid)
        imgname = '{0}_05_5_img.fits'.format(obsid)
        psf_bv_list=output_psf_intervalue(ecf,psf_x,psf_y,img_x,img_y,if_ecflist=True)
        psf_bestSNR_no_overlap=[];ECF_allsrc=[];SNR_allsrc=[]
        (phy_x,phy_y)=funcs.trans_radec2xy(inpath+imgname,ra,dec)
        for i in range(len(phy_x)):
            dist=np.sqrt((phy_x-phy_x[i])**2+(phy_y-phy_y[i])**2)
            close_dist=np.sort(dist)[1]

            SNR=[]
            bkg_reg = [phy_x[i], phy_y[i], 2 * psf_bv_list[:, i][10], 4 * psf_bv_list[:, i][10]]
            ##默认annulus用2-4倍psf90
            [bkg_t, bkg_E, bkg_ID, bkg_area] = funcs.get_evt_bkgreg(inpath + evtname,inpath+ imgname,obsid, bkg_reg,
                                                                    ra_list=ra, dec_list=dec,
                                                                    radius_list=psf_bv_list[10],
                                                                    emin=500, emax=5000)
            # print(len(bkg_t), bkg_area)
            if bkg_area == 0:
                SNR.append(0)
                continue

            overlap_ecf=[]
            for j in range(len(ecflist)):
                src_reg=[phy_x[i],phy_y[i],psf_bv_list[:,i][j]]
                [src_t, src_E, src_ID] = funcs.get_evt_srcreg(inpath+evtname,obsid,src_reg,emin=500,emax=5000)

                src_cts=len(src_t)
                bkg_cts=len(bkg_t)
                src_area=np.pi*src_reg[2]**2
                SNR.append((src_cts-bkg_cts*src_area/bkg_area)/np.sqrt(src_cts))

            src_bestSNR_psf_radius = psf_bv_list[:, i][np.argmax(SNR)]  ##用不到，先留着
            ##根据overlap的值再缩小半径##
            if close_dist < 2 * psf_bv_list[:, i][10]:
                #按照2倍psf75半径为overlap
                overlap_dist=2 * psf_bv_list[:, i]-close_dist
                overlap_index=funcs.find_last_negative(overlap_dist)

            else:
                overlap_index=np.argmax(SNR)
                if overlap_index==11:overlap_index=0
            src_final_SNR_no_overlap=psf_bv_list[:,i][overlap_index]

            print(src_bestSNR_psf_radius,src_final_SNR_no_overlap,ecflist[overlap_index],SNR[overlap_index])

            psf_bestSNR_no_overlap.append(src_final_SNR_no_overlap)
            ECF_allsrc.append(ecflist[overlap_index])
            SNR_allsrc.append(SNR[overlap_index])
        psf_bestSNR_no_overlap=np.array(psf_bestSNR_no_overlap)
        ECF_allsrc=np.array(ECF_allsrc);SNR_allsrc=np.array(SNR_allsrc)

        return (psf_bestSNR_no_overlap,ECF_allsrc,SNR_allsrc,not_include_index)

def make_region(srcID,ra,dec,radius,ecf,obsIDlist,outpath):
    if ecf=='SNR':
        ecfstr='SNRpsf'
    else:
        ecfstr=str(int(ecf*100))
    radius=np.array(radius)/20/3600
    for obsid in obsIDlist:
        os.chdir(outpath)
        if os.path.exists('reg_{0}/region_{1}'.format(obsid,ecfstr)):
            os.chdir('reg_{0}/region_{1}'.format(obsid,ecfstr))
        else:
            os.mkdir('reg_{0}/region_{1}'.format(obsid,ecfstr))
            os.chdir('reg_{0}/region_{1}'.format(obsid,ecfstr))
        with open('./all.reg', 'a+') as f2:
            f2.writelines('fk5' + '\n')
        for i in range(len(ra)):
            with open('./{0}.reg'.format(srcID[i]), 'w+') as f1:
                f1.writelines('fk5'+'\n')
                reg = 'circle(' + str(ra[i]) + ',' + str(dec[i]) + ',' + str(radius[i]) + ')'
                f1.writelines(reg)
            with open('./all.reg', 'a+') as f2:
                f2.writelines(reg +' # text={'+str(srcID[i])+'}'+ '\n')

def write_cat_psfinfo(srcID,ra,dec,radius,ecf,SNR,not_include_index,obsIDlist,outpath):
    ##输入的radius单位统一为image的pixel，乘4才是角秒##
    os.chdir(outpath)
    for obsid in obsIDlist:
        if os.path.exists('reg_{0}/region_SNRpsf'.format(obsid)):
            print('OK')
        else:
            os.mkdir('reg_{0}/region_SNRpsf'.format(obsid))
        info=np.column_stack((srcID,ra,dec,radius/20.,ecf,SNR,not_include_index))
        np.savetxt('reg_{0}/region_SNRpsf/src_info_{0}.txt'.format(obsid),info,fmt='%8d %10.5f %10.5f %10.5f %10.2f %10.5f %5d')

def get_all_src_txt(catalog_file,obsIDlist,inpath):
    ###------read catalog file, need modification---##
    (ra_list,dec_list,srcIDlist)=funcs.read_erosita_cat('/Users/baotong/Desktop/period_Tuc/erosita_cat_coord.xlsx')
    # ###-------------------------------------------##
    if len(ra_list)!=len(srcIDlist) or len(dec_list)!=len(srcIDlist):
        print('ERROR!')
        return None
    for j in range(len(obsIDlist)):
        srcarea=[]
        (psf_bettervalue, not_include_index) = get_psfradius(srcID=srcIDlist, ra=ra, dec=dec, obsid=obsIDlist[j], inpath=path,
                                                             outpath=path, ecf=0.90, if_SN_radius=False)
        outpath = inpath + 'txt/txt_psf90_{0}/'.format(obsIDlist[j])
        evtfile = inpath + 'pm00_{0}_020_EventList_c001_bary.fits'.format(obsIDlist[j])
        imgfilename =inpath+ '{0}_05_5_img.fits'.format(obsIDlist[j])
        for i in range(len(ra_list)):
            evtfile=inpath+'pm00_{0}_020_EventList_c001_bary.fits'.format(obsIDlist[j])
            get_singlesrc_evt(evtfile=evtfile,imgname=imgfilename,obsid=obsIDlist[j],ra=ra_list[i],dec=dec_list[i],
                              radius=psf_bettervalue[i], emin=200,emax=5000,outpath=outpath,outname=str(srcIDlist[i]))
            srcarea.append(np.pi*psf_bettervalue[i]**2)
        np.savetxt(outpath + 'src_area.txt', np.column_stack((srcIDlist, np.array(srcarea))))
    # for i in range(len(ra_list)):
    #     merge_txt(obsIDlist=obsIDlist,inpath=inpath+'txt/',outpath=inpath+'txt/txt_merge_psf50_0.2_5/',outname=str(srcIDlist[i]),epoch_file='epoch_47Tuc.txt')

def get_all_src_txt_SNRpsf(catalog_file,obsIDlist,inpath):
    (ra_list, dec_list, srcIDlist) =funcs.read_erosita_cat(catalog_file)
    if len(ra_list)!=len(srcIDlist) or len(dec_list)!=len(srcIDlist):
        print('ERROR!')
        return None
    # for j in range(len(obsIDlist)):
    #     print(obsIDlist[j])
    #     evtfile = inpath + 'pm00_{0}_020_EventList_c001_bary.fits'.format(obsIDlist[j])
    #     imgfilename =inpath+ '{0}_05_5_img.fits'.format(obsIDlist[j])
    #     src_info=np.loadtxt(inpath+'reg_{0}/region_SNRpsf/src_info_{0}.txt'.format(obsIDlist[j]))
    #     radius_list=src_info[:,3];not_include_index=src_info[:,6]
    #     radius_list*=20.
    #     outpath=inpath+'txt/txt_SNR_{0}/'.format(obsIDlist[j])
    #     for i in range(len(ra_list)):
    #         if not_include_index[i]==1:
    #             get_singlesrc_evt(evtfile=evtfile,imgname=imgfilename, obsid=obsIDlist[j], ra=ra_list[i], dec=dec_list[i],
    #                               radius=radius_list[i], emin=500,emax=5000,outpath=outpath, outname=str(srcIDlist[i]))

    for i in range(len(ra_list)):
        merge_txt(obsIDlist=obsIDlist[0:4],inpath=inpath+'txt/',outpath=inpath+'txt/txt_merge_psf75_0.5_5/',outname=str(srcIDlist[i]),epoch_file='epoch_47Tuc.txt')

def get_all_bkg_txt(catalog_file,obsIDlist,inpath):
    ###------read catalog file, need modification---##
    (ra_list,dec_list,srcIDlist)=funcs.read_erosita_cat('/Users/baotong/Desktop/period_Tuc/erosita_cat_coord.xlsx')
    ###-------------------------------------------##
    if len(ra_list)!=len(srcIDlist) or len(dec_list)!=len(srcIDlist):
        print('ERROR!')
        return None
    # for j in range(len(obsIDlist)):
    #     bkgarea_all=[]
    #     (psf_bettervalue, not_include_index) = get_psfradius(srcID=srcIDlist, ra=ra, dec=dec, obsid=obsIDlist[j], inpath=path,
    #                                                          outpath=path, ecf=0.9, if_SN_radius=False)
    #     (psf_bettervalue_src, not_include_index_2) = get_psfradius(srcID=srcIDlist, ra=ra, dec=dec, obsid=obsIDlist[j], inpath=path,
    #                                                          outpath=path, ecf=0.5, if_SN_radius=False)
    #     outpath = inpath + 'txt/txt_psf75_{0}/'.format(obsIDlist[j])
    #     evtfile = inpath + 'pm00_{0}_020_EventList_c001_bary.fits'.format(obsIDlist[j])
    #     imgfilename =inpath+ '{0}_05_5_img.fits'.format(obsIDlist[j])
    #     for i in range(len(ra_list)):
    #         evtfile=inpath+'pm00_{0}_020_EventList_c001_bary.fits'.format(obsIDlist[j])
    #         bkgarea_temp=get_singlebkg_evt(evtfile=evtfile,imgname=imgfilename,obsid=obsIDlist[j],ra=ra_list[i],dec=dec_list[i],
    #                           emin=200,emax=5000,radius_src=psf_bettervalue[i],radius_list=psf_bettervalue_src,
    #                           ra_list=ra_list,dec_list=dec_list,outpath=outpath,outname=str(srcIDlist[i]))
    #         bkgarea_all.append(bkgarea_temp)
    #     np.savetxt(outpath+'bkg_area.txt',np.column_stack((srcIDlist,np.array(bkgarea_all))))


    for i in range(len(ra_list)):
        merge_txt(obsIDlist=obsIDlist, inpath=inpath + 'txt/', outpath=inpath + 'txt/txt_merge_psf50_0.2_5/',
                  outname=str(srcIDlist[i]), epoch_file='epoch_47Tuc.txt')

def filter_obs(src_evt,useid):
    src_evt_use = src_evt[np.where(src_evt[:-1] == useid[0])[0]]
    i=1
    while i < len(useid):
        id=useid[i]
        src_evt_use_temp=src_evt[np.where(src_evt[:-1]==id)[0]]
        src_evt_use = np.concatenate((src_evt_use, src_evt_use_temp))
        i+=1
    return src_evt_use

def merge_4reg_txt(srcIDlist):
    path = '/Users/baotong/eSASS/data/raw_data/47_Tuc/txt/'
    obs_reg1 = np.array([700011, 700173])
    obs_reg2 = np.array([700163, 700174])
    obs_reg3 = np.array([700013, 700175])
    obs_reg4 = np.array([700014])
    obs_04= np.array([700011,700163,700013,700014])
    obsIDlist=[obs_reg1,obs_reg2,obs_reg3,obs_reg4]
    # for k in range(len(obsIDlist)):
    #     for i in range(len(srcIDlist)):
    #         merge_txt(obsIDlist[k], inpath=path, outpath=path + 'txt_reg{0}_psf50_0.2_5/'.format(k+1),
    #                   outname=str(srcIDlist[i]), epoch_file='epoch_47Tuc_reg{0}.txt'.format(k+1))
    for i in range(len(srcIDlist)):
        merge_txt(obs_04, inpath=path, outpath=path + 'txt_obs04_psf50_0.2_5/',
                  outname=str(srcIDlist[i]), epoch_file='epoch_47Tuc_obs04.txt')
if __name__=='__main__':
    # path = '/Users/baotong/eSASS/data/47_Tuc/'
    path='/Users/baotong/eSASS/data/raw_data/47_Tuc/'
    ID=[700011,700013,700014,700163,700173,700174,700175]
    catalog_file='/Users/baotong/Desktop/period_Tuc/erosita_cat_coord.xlsx'
    # make_epochfile(path=path,ID=ID,outpath=path+'txt/')

    # for obsid in ID:
    #     get_singlesrc_evt(path+'pm00_{0}_020_EventList_c001_bary.fits'.format(obsid),obsid=obsid,
    #                   ra=6.04485,dec=-72.07383,radius=20./3600,outpath='/Users/baotong/eSASS/data/47_Tuc/txt/',
    #                   outname='402')
    #
    # merge_txt(ID,outpath='/Users/baotong/eSASS/data/47_Tuc/txt/',outname='402')

    # get_psfradius(ra=np.array([6.0178,6.1078]),dec=np.array([-72.08281,-72.18281]),obsid=700011,inpath=path,outpath=path)
    (ra,dec,srcIDlist)=funcs.read_erosita_cat(catalog_file)

    # get_all_src_txt(catalog_file, obsIDlist=ID, inpath=path)
    get_all_bkg_txt(catalog_file, obsIDlist=ID, inpath=path)

    # get_all_src_txt_SNRpsf(catalog_file, obsIDlist=ID, inpath=path)
    # merge_4reg_txt(srcIDlist)

    # for obsID in ID[0:]:
    #     (psf_bestSNR_no_overlap,ECF_allsrc,SNR_allsrc,not_include_index)=get_psfradius(srcID=srcIDlist[0:],ra=ra[0:], dec=dec[0:], obsid=obsID, inpath=path,
    #                   outpath=path, ecf=0.9,if_SN_radius=True)
    #     write_cat_psfinfo(srcID=srcIDlist[0:],ra=ra[0:],dec=dec[0:],radius=psf_bestSNR_no_overlap,
    #                       ecf=ECF_allsrc,SNR=SNR_allsrc,not_include_index=not_include_index,obsIDlist=[obsID],outpath=path)
    #     make_region(srcID=srcIDlist[0:],ra=ra[0:],dec=dec[0:],ecf='SNR',radius=psf_bestSNR_no_overlap,obsIDlist=[obsID],outpath=path)

    # ##make all region## ##
    # for obsID in ID:
    #     (psf_bettervalue,not_include_index)=get_psfradius(srcID=srcIDlist,ra=ra, dec=dec, obsid=obsID, inpath=path,
    #                   outpath=path, ecf=0.90,if_SN_radius=False)
    #
    #     make_region(srcID=srcIDlist,ra=ra,dec=dec,ecf=0.90,radius=psf_bettervalue,obsIDlist=[obsID],outpath=path)