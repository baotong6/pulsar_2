#!/bin/bash
# -*- coding: utf-8 -*-
#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
from astropy.io import fits
from astropy.wcs import WCS

# ------choose obsID-------------
obsList=["0108060401","0108060501","0108060601","0108060701","0108061801","0108061901","0108062101",
         "0108062301","0555780101","0555780201","0555780301","0555780401","0555780501","0555780601",
         "0555780701","0555780801","0555780901","0555781001","0555782301","0604960101","0604960201",
         "0604960301","0604960401","0604961101","0604961201","0604960701","0604960501","0604961301",
         "0604960601","0604960801","0604960901","0604961001","0604961801"]
# obsList=["0604960301"]
# ------detName List-------------
det1 = "mos1"
det2 = "mos2"
det3 = "mos"
# -------------------------------
#
# ------choose det --------------
detList = [det1,det2,det3]

path = '/Users/baotong/Desktop/CDFS/xmm_CDFS/'
os.chdir(path)
filename = 'XMM_CDFS.fit'
catalog = fits.open(filename)
seq = catalog[1].data['recno']
ra = catalog[1].data['RAJ2000']
dec = catalog[1].data['DEJ2000']
counts = catalog[1].data['counts']
source_info = np.column_stack((seq, ra, dec, counts))
source_info = source_info[np.lexsort(-source_info.T)]
# ra = source_info[:, 1];
# dec = source_info[:, 2]
# seq = source_info[:, 0]
# seq = seq.astype('int')

##----single source-------##
ra=[53.143475]
dec=[-27.653576]
seq=['XID643']
def make_txt_from_evt_coord(evtfile,imgfile,epochfile,outpath,seq,ra,dec):
    w = WCS(imgfile)
    img_x, img_y = w.all_world2pix(ra, dec, 1)
    print(img_x)
    phy_x=img_x*87-43
    phy_y=img_y*87-43

    hdul_evt=fits.open(evtfile)
    x=hdul_evt[1].data['X']
    y=hdul_evt[1].data['Y']

    energy=hdul_evt[1].data['PI']
    time=hdul_evt[1].data['TIME']


    def where_region(x,y,reg):
        r=np.array((x-reg[0],y-reg[1]))
        len_r=np.sqrt(r[0]**2+r[1]**2)
        temp=len_r-reg[2]
        return np.where(temp<=0)

    def make_region():
        for i in range(len(phy_x)):
            # with open('./region_90/{0}.reg'.format(i + 1), 'w+') as f1:
            with open('/Volumes/pulsar/xmmCDFS/merge_evt_txt_fromfits/region_10arcsec/{0}.reg'.format(seq[i]), 'w+') as f1:
                reg = 'circle(' + str(phy_x[i]) + ',' + str(phy_y[i]) + ',' + str(200) + ')'
                f1.writelines(reg)
            with open('/Volumes/pulsar/xmmCDFS/merge_evt_txt_fromfits/region_10arcsec/all.reg', 'a+') as f2:
                f2.writelines(reg + '\n')

            with open('/Volumes/pulsar/xmmCDFS/merge_evt_txt_fromfits/region_10arcsec/all_name.reg', 'a+') as f2:
                reg = 'circle(' + str(phy_x[i]) + ',' + str(phy_y[i]) + ',' + str(200) + ')'+'   # text={'+str(seq[i])+'}'
                f2.writelines(reg + '\n')

            with open('/Volumes/pulsar/xmmCDFS/merge_evt_txt_fromfits/region_10arcsec/wcs_all.reg', 'a+') as f3:
                reg = 'circle(' + str(ra[i]) + ',' + str(dec[i]) + ',' + str(10)+'"' + ')'
                f3.writelines(reg + '\n')
        return None
    # make_region()
    epoch = np.loadtxt(epochfile)
    TSTART=epoch[:,0]
    TSTOP=epoch[:,1]
    obs_ID=epoch[:,2]

    for i in range(len(ra)):
        reg = [phy_x[i], phy_y[i], 240]
        src_cts_index = where_region(x, y, reg)
        src_cts_t = time[src_cts_index]
        src_cts_E = energy[src_cts_index]
        src_obsID=np.zeros(len(src_cts_t))
        obsindex_del = []
        for k in range(len(obs_ID)):
            src_obsID[np.where((src_cts_t>TSTART[k])&(src_cts_t<TSTOP[k]))]+=int(obs_ID[k])
            counts = len(np.where((src_cts_t > TSTART[k]) & (src_cts_t < TSTOP[k]))[0])
            if counts == 0: obsindex_del.append(k)
        epoch_out = np.delete(epoch, obsindex_del, axis=0)
        source_info=np.column_stack((src_cts_t,src_cts_E,src_obsID))
        print('run source: {0}'.format(i+1))
        np.savetxt(outpath + 'NID'+str(seq[i]) + '_mos_12sec.txt', source_info, fmt="%20.7f  %20.3f  %20d")

        np.savetxt(outpath + 'epoch_' + 'NID'+str(seq[i]) + '_mos_12sec.txt', epoch_out, fmt="%15.2f  %15.2f %15d %15.3f")

    return None


def filter_evt(detList):
    ##可能哪里写错了，别用！！##
    # w = WCS(path + 'image_ep123456_all_merge_filt_time.fits')
    # src_x, src_y = w.all_world2pix(ra[0], dec[0], 1)
    # img_x = src_x
    # img_y = src_y
    path = "/Volumes/pulsar/xmmCDFS"
    for i in range(102,len(ra)):
        srcName = "NID{0}".format(seq[i])
        srcReg = "circle({0},{1},0.002777777777777778)".format(ra[i],dec[i])
        bkgReg = "circle(52.9715056,-27.8814980,0.002777777777777778)"
        for det in detList:
            datapath = path + "/"+'merge_evt/'
            os.chdir(datapath)
            print(datapath)
            # os.environ['SAS_CCF'] = path + "/" + obsID + "/cal/ccf.cif"
            # os.chdir(datapath)
            # with open("SAS.txt", 'r') as f:
            #     sasname = f.readline()
            # print(sasname)
            # os.environ['SAS_ODF'] = path + "/" + obsID + "/cal/" + sasname[0:-1]

            # print("running obsservation" + " " + obsID)
            print("running detector" + " " + det)

            cmd = "evselect table=ep123456_all_" + det + "_filt_time.fits energycolumn=PI " \
                                                "expression='#XMMEA_EP && (PATTERN<=4) && ((RA,DEC) IN " + srcReg + ")" + " &&(PI in [200:12000])' withfilteredset=yes filteredset=" \
                      + det + '_' + srcName + "_src_filt_time_10sec.fits filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes"
            print(" ")
            print("1 extract src event list")
            print(cmd)
            os.system(cmd)

            # if det == "mos":
            #     cmd = "evselect table=" + det + "_filt_time.fits energycolumn=PI " \
            #                                     "expression='#XMMEA_EP && (PATTERN<=4) && ((RA,DEC) IN " + bkgReg + ")" + " &&(PI in [200:12000])' withfilteredset=yes filteredset=" \
            #           + det + '_' + srcName + "_bkg_filt_time_20sec.fits filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes"
            # else:
            #     cmd = "evselect table=" + det + "_filt_time.fits energycolumn=PI " \
            #                                     "expression='#XMMEA_EM && (PATTERN<=12) && ((RA,DEC) IN " + bkgReg + ")" + " &&(PI in [200:12000])' withfilteredset=yes filteredset=" \
            #           + det + '_' + srcName + "_bkg_filt_time_20sec.fits filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes"
            # print(" ")
            # print("2 extract bkg event list")
            # print(cmd)
            # os.system(cmd)
    return None


def get_txt_mergeevt(evtfile,srcName,outpath,epochfile,band=[200,10000]):
    if os.path.exists(evtfile):
        hdul_evt= fits.open(evtfile)
        epoch=np.loadtxt(epochfile)

        energy=hdul_evt[1].data['PI']
        time=hdul_evt[1].data['TIME']

        src_t=time;src_E=energy
        obsindex_del=[]
        for i in range(len(epoch)):
            TSTART=epoch[i][0]
            TSTOP=epoch[i][1]

            counts=len(np.where((time>TSTART)&(time<TSTOP))[0])
            if counts==0:obsindex_del.append(i)
        epoch_out=np.delete(epoch,obsindex_del,axis=0)



        def delete_photon_ID(time,energy,band):
            i=0
            while i < len(energy):
                if energy[i]>band[1] or energy[i]<band[0]:
                    energy=np.delete(energy,i)
                    time=np.delete(time,i)
                    i=i-1
                i=i+1
            return [time,energy]

        [src_t,src_E]=delete_photon_ID(src_t,src_E,band)

        src_txt = np.column_stack((src_t,src_E))
        src_txt = src_txt[src_txt[:,0].argsort()]

        np.savetxt(outpath  +srcName+'_mos_10sec.txt', src_txt, fmt="%.7f  %.3f ")
        # np.savetxt(outpath + obsID + '/txt/' + srcName + '_' + mode + '_bkg_10sec.txt', bkg_txt, fmt="%.7f  %.3f ")
        np.savetxt(outpath + 'epoch_'+srcName+'_mos_10sec.txt', epoch_out,fmt="%.2f  %.2f %15d %15.3f")
    else:
        print("No such file: {0}".format(evtfile))

    return None

if __name__ == "__main__":

    # epochfile='epoch_mos_xmmobs.txt'
    # for k in range(91,len(ra)):
    #     srcName = "NID{0}".format(seq[k])
    #     evtfile='/Volumes/pulsar/xmmCDFS/merge_evt/'+'mos_NID{0}_src_filt_time_10sec.fits'.format(seq[k])
    #     outpath='/Volumes/pulsar/xmmCDFS/merge_evt_txt/'
    #     get_txt_mergeevt(evtfile,srcName,outpath,outpath+epochfile,[500,8000])
    path='/Volumes/pulsar/xmmCDFS/merge_evt/'
    outpath='/Volumes/pulsar/xmmCDFS/merge_evt_txt_fromfits/epoch5_txt/'
    evtfile=path+'ep5_all_mos_filt_time.fits'
    imgfile=path+'ep123456_all_mos_filt_time.img'
    epochfile=path+'epoch_mos_xmmobs.txt'

    make_txt_from_evt_coord(evtfile,imgfile,epochfile,outpath,seq,ra,dec)

