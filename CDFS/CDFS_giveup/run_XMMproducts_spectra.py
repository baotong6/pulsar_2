#!/usr/bin/env python3
#written by by Tong
# written on 2018-6-13 to run XMM products on SgrA PWN for 5 XMM obsID
# modified on 2019-3-9 to be the pipeline for point source (dummpy users version,not recommended for research)
# modified on 2020-4-7 to add the extraction of light curve
# modified on 2020-4-8 to add the particle background filtering
# point source or extended source: change flag in afrgen: extendedsource=yes/no; currently set to no for point source
import sys
import os
import string
from pathlib import Path
from astropy.io import fits
# ------obsID List---------------
# path="/Users/baotong/xmm/M28_LMXB"
path="/Volumes/pulsar/xmmCDFS"
# obsID1 = "0506440101"
# -------------------------------
# ------choose obsID-------------
obsList=["0108060401","0108060501","0108060601","0108060701","0108061801","0108061901","0108062101",
         "0108062301","0555780101","0555780201","0555780301","0555780401","0555780501","0555780601",
         "0555780701","0555780801","0555780901","0555781001","0555782301","0604960101","0604960201",
         "0604960301","0604960401","0604961101","0604961201","0604960701","0604960501","0604961301",
         "0604960601","0604960801","0604960901","0604961001","0604961801"]
# obsList=["0604960301"]
# 0555780501
# obsList=["0108060601","0108060701","0108061801","0108061901","0108062101",
#          "0108062301","0555780101","0555780201","0555780301","0555780401","0555780501","0555780601",
#          "0555780701","0555780801","0555780901","0555781001","0555782301","0604960101","0604960201",
#          "0604960301","0604960401","0604960501","0604960601","0604960701","0604960801","0604960901",
#          "0604961001","0604961101","0604961201","0604961301","0604961801"]
# obsList=["0555780601","0555780701","0555780801","0555780901","0555781001","0555782301","0604960101","0604960201",
#          "0604960301","0604960401","0604960501","0604960601","0604960701","0604960801","0604960901",
#          "0604961001","0604961101","0604961201","0604961301","0604961801"]
# obsList=['0108062101']
# -------------------------------
# ------detName List-------------
det1 = "mos1"
det2 = "mos2"
det3 = "pn"
# -------------------------------
#
# ------choose det --------------
detList = [det1,det2,det3]
# -------------------------------
process=0
filt_particle_bkg=0
lc=0
filter_evt=0
spectra=0
combine_spec=0
merge_evt=0

ra=53.1457;
dec=-27.8252

for obsID in obsList:
   os.chdir(path+"/"+obsID)
   mypath=Path("./cal")
   mypath=Path("./cal")

   if process:
      if mypath.is_dir():
         print("continued")
      else:
         os.system("mkdir cal")
         os.system("gunzip ODF/*.gz")
         #------------ set environment -----------------
         os.chdir("./cal")
         #cmd = "alias SAS_ODF=" + path +"/"+ obsID + "/ODF"
         os.environ['SAS_ODF'] = path +"/"+ obsID + "/ODF"
         os.system("cifbuild")
         os.environ['SAS_CCF'] = path+"/"+obsID+"/cal/ccf.cif"
         os.system("odfingest")
         os.system("rm SAS.txt")
         os.system("ls *.SAS >SAS.txt")
         with open("SAS.txt",'r') as f:
            sasname=f.readline()
         print(sasname)
         os.environ['SAS_ODF'] = path + "/" + obsID + "/cal/"+sasname[0:-1]
         # # ---------------------------------------------
         # -----------processing for evt2 file----------
         os.system("emchain")
         os.system("epchain")
         # # ---------------------------------------------
         #------------rename file-----------------------
         os.system("rm *FIT.txt")
         v1=os.popen("ls *M1*EVLI*.FIT")
         mos1name=v1.read()[0:-1]
         v2=os.popen("ls *M2*EVLI*.FIT")
         mos2name=v2.read()[0:-1]
         v3=os.popen("ls *PI*EVLI*.FIT")
         pnname=v3.read()[0:-1]

         cmd="cp "+mos1name+" mos1.fits"
         os.system(cmd)
         cmd="cp "+mos2name+" mos2.fits"
         os.system(cmd)
         cmd="cp "+pnname+" pn.fits"
         os.system(cmd)
         # # ---------------------------------------------
         # #-----------barycen----------------------------
         for det in detList:
            cmd="cp "+det+".fits "+det+"_bary.fits"
            print(cmd)
            os.system(cmd)
            cmd="barycen withtable=yes table="+det+"_bary.fits"+": withsrccoordinates=yes srcra="+str(ra)+" srcdec="+str(dec)
            print(cmd)
            os.system(cmd)
         # # ---------------------------------------------
   #---------choose region-------------------
   #You should open pn.fits and choose the region, saved in physical coordinate.
   # # ---------------------------------------------
   #---------reg def------------------------------------
   # ---------------------------------------------
   srcName = "XID19"
   srcReg = "circle(52.9469910,-27.8870110,0.005555555555555556)"
   bkgReg = "circle(52.9715056,-27.8814980,0.005555555555555556)"

   if filt_particle_bkg:
      pn_threshold=0.5;mos_threshold=0.4
      for det in detList:
         print("running obsservation"+" "+obsID)
         print("running detector"+" "+det)

         datapath = path+"/"+obsID+"/cal/"
         print(datapath)
         os.environ['SAS_CCF'] = path + "/" + obsID + "/cal/ccf.cif"
         os.chdir(datapath)

         if det == "pn":
             cmd = "evselect table=pn_bary.fits withrateset=Y rateset=rateEPIC_pn.fits maketimecolumn=Y " \
                   "timebinsize=100 makeratecolumn=Y expression='#XMMEA_EP && (PI>10000&&PI<12000) && (PATTERN==0)'"

         else:
            cmd = "evselect table={0}_bary.fits withrateset=Y rateset=rateEPIC_{0}.fits maketimecolumn=Y " \
                   "timebinsize=100 makeratecolumn=Y expression='#XMMEA_EM && (PI>10000) && (PATTERN==0)'".format(det)
         print(" ")
         print("1 make lc")
         print(cmd)
         os.system(cmd)

         if det == "pn":
             cmd = "tabgtigen table=rateEPIC_pn.fits expression='RATE<={0}' gtiset=gti_pn.fits".format(pn_threshold)

         else:
            cmd = "tabgtigen table=rateEPIC_{0}.fits expression='RATE<={1}' gtiset=gti_{0}.fits".format(det,mos_threshold)
         print(" ")
         print("2 make GTI file")
         print(cmd)
         os.system(cmd)


         if det == "pn":
             cmd = "evselect table=pn_bary.fits withfilteredset=yes expression='#XMMEA_EP && (PATTERN<=4) && gti(gti_pn.fits,TIME) && (PI>200)' " \
                   "filteredset=pn_filt_time.fits filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes"
         else:
            cmd = "evselect table={0}_bary.fits withfilteredset=yes expression=" \
                  "'#XMMEA_EM && (PATTERN<=12) && gti(gti_{0}.fits,TIME) && (PI>200)' filteredset={0}_filt_time.fits filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes".format(det)

         print(" ")
         print("3 filter GTI fits")
         print(cmd)
         os.system(cmd)

   if lc:
      # you should run this step multiple times to determine the best tmin and tmax
      lenbin=100;tmin=0;tmax=0
      ## you can also specify the tmin and tmax (read from lc)
      ## otherwise they would be total duration
      for det in detList:
         print("running obsservation"+" "+obsID)
         print("running detector"+" "+det)

         datapath = path+"/"+obsID+"/cal/"
         print(datapath)
         os.environ['SAS_CCF'] = path + "/" + obsID + "/cal/ccf.cif"
         os.chdir(datapath)
         with open("SAS.txt",'r') as f:
            sasname=f.readline()
         print(sasname)
         os.environ['SAS_ODF'] = path + "/" + obsID + "/cal/"+sasname[0:-1]

         # ##------barycen-----##
         # cmd = "cp " + det + "_filt_time.fits " + det + "_filt_time_bary.fits"
         # print(cmd)
         # os.system(cmd)
         # cmd = "barycen withtable=yes table=" + det + "_filt_time_bary.fits" + ": withsrccoordinates=yes srcra=" + str(
         #    ra) + " srcdec=" + str(dec)
         # print(cmd)
         # os.system(cmd)
         # print("0 barycenter correction for cleaned events")

         ##------------------##
         if tmin==0 and tmax==0:
            event=fits.open(det+"_filt_time.fits")
            tmin=event[1].header['TSTART']
            tmax=event[1].header['TSTOP']

         if det == "pn":
             cmd = "evselect table="+det+"_filt_time.fits energycolumn=PI " \
                                                  "expression='#XMMEA_EP && (PATTERN<=4) && ((RA,DEC) IN "+srcReg+ ")"+" &&(PI in [500:8000])' withrateset=yes rateset="\
                   +det+'_'+srcName+"_src_lc_bin{0}_0.5_8.lc timebinsize={0} maketimecolumn=yes makeratecolumn=yes timemin={1} timemax={2}".format(lenbin,tmin,tmax)
         else:
            cmd = "evselect table=" + det + "_filt_time.fits energycolumn=PI " \
                                                       "expression='#XMMEA_EM && (PATTERN<=12) && ((RA,DEC) IN " + srcReg + ")"+" &&(PI in [500:8000])' withrateset=yes rateset=" \
                  + det +'_'+srcName+ "_src_lc_bin{0}_0.5_8.lc timebinsize={0} maketimecolumn=yes makeratecolumn=yes timemin={1} timemax={2}".format(lenbin, tmin, tmax)
         print(" ")
         print("1 extract src light curve")
         print(cmd)
         os.system(cmd)

         if det == "pn":
             cmd = "evselect table="+det+"_filt_time.fits energycolumn=PI " \
                                                  "expression='#XMMEA_EP && (PATTERN<=4) && ((RA,DEC) IN "+bkgReg+")"+" &&(PI in [500:8000])' withrateset=yes rateset="\
                   +det+'_'+srcName+"_bkg_lc_bin{0}_0.5_8.lc timebinsize={0} maketimecolumn=yes makeratecolumn=yes timemin={1} timemax={2}".format(lenbin,tmin,tmax)
         else:
            cmd = "evselect table=" + det + "_filt_time.fits energycolumn=PI " \
                                                       "expression='#XMMEA_EM && (PATTERN<=12) && ((RA,DEC) IN " + bkgReg + ")"+" &&(PI in [500:8000])' withrateset=yes rateset=" \
                  + det + '_'+srcName+"_bkg_lc_bin{0}_0.5_8.lc timebinsize={0} maketimecolumn=yes makeratecolumn=yes timemin={1} timemax={2}".format(lenbin, tmin, tmax)
         print(" ")
         print("2 extract bkg light curve")
         print(cmd)
         os.system(cmd)


         cmd ="epiclccorr srctslist={0}_{2}_src_lc_bin{1}_0.5_8.lc eventlist={0}_bary.fits outset={0}_{2}_lccorr_bin{1}_0.5_8.lc bkgtslist={0}_{2}_bkg_lc_bin{1}_0.5_8.lc withbkgset=yes applyabsolutecorrections=yes".format(det,lenbin,srcName)

         print(" ")
         print("3 extract corrected light curve")
         print(cmd)
         os.system(cmd)
   if filter_evt:
      for det in detList:

         datapath = path+"/"+obsID+"/cal/"
         print(datapath)
         os.environ['SAS_CCF'] = path + "/" + obsID + "/cal/ccf.cif"
         os.chdir(datapath)
         with open("SAS.txt",'r') as f:
            sasname=f.readline()
         print(sasname)
         os.environ['SAS_ODF'] = path + "/" + obsID + "/cal/"+sasname[0:-1]

         print("running obsservation"+" "+obsID)
         print("running detector"+" "+det)

         if det == "pn":
            cmd = "evselect table=" + det + "_filt_time.fits energycolumn=PI " \
                                            "expression='#XMMEA_EP && (PATTERN<=4) && ((RA,DEC) IN " + srcReg + ")" + " &&(PI in [200:12000])' withfilteredset=yes filteredset=" \
                  + det + '_' + srcName + "_src_filt_time_20sec.fits filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes"
         else:
            cmd = "evselect table=" + det + "_filt_time.fits energycolumn=PI " \
                                            "expression='#XMMEA_EM && (PATTERN<=12) && ((RA,DEC) IN " + srcReg + ")" + " &&(PI in [200:12000])' withfilteredset=yes filteredset=" \
                  + det + '_' + srcName + "_src_filt_time_20sec.fits filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes"
         print(" ")
         print("1 extract src event list")
         print(cmd)
         os.system(cmd)

         if det == "pn":
            cmd = "evselect table=" + det + "_filt_time.fits energycolumn=PI " \
                                            "expression='#XMMEA_EP && (PATTERN<=4) && ((RA,DEC) IN " + bkgReg + ")" + " &&(PI in [200:12000])' withfilteredset=yes filteredset=" \
                  + det + '_' + srcName + "_bkg_filt_time_20sec.fits filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes"
         else:
            cmd = "evselect table=" + det + "_filt_time.fits energycolumn=PI " \
                                            "expression='#XMMEA_EM && (PATTERN<=12) && ((RA,DEC) IN " + bkgReg + ")" + " &&(PI in [200:12000])' withfilteredset=yes filteredset=" \
                  + det + '_' + srcName + "_bkg_filt_time_20sec.fits filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes"
         print(" ")
         print("2 extract bkg event list")
         print(cmd)
         os.system(cmd)




   if spectra:
      for det in detList:
         print("running obsservation"+" "+obsID)
         print("running detector"+" "+det)
         filt_label='_filt_time'

         datapath = path+"/"+obsID+"/cal/"
         print(datapath)
         os.environ['SAS_CCF'] = path + "/" + obsID + "/cal/ccf.cif"

         if det == "pn":
             cmd = "evselect table="+datapath+det+filt_label+".fits withfilteredset=yes expression='(PATTERN <= 4)&&(PI in [200:15000])&&#XMMEA_EP' filteredset="+datapath+det+filt_label+"_filt.fits filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes"
         else:
             cmd = "evselect table="+datapath+det+filt_label+".fits withfilteredset=yes expression='(PATTERN <= 12)&&(PI in [200:12000])&&#XMMEA_EM' filteredset="+datapath+det+filt_label+"_filt.fits filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes"
         print(" ")
         print("1 filter by energy")
         print(cmd)
         os.system(cmd)
         cmd = "evselect table="+datapath+det+filt_label+"_filt.fits withimageset=yes imageset="+datapath+det+filt_label+"_filt.img xcolumn=X ycolumn=Y imagebinning=imageSize ximagesize=600 yimagesize=600"
         print(" ")
         print("2 create det image")
         print(cmd)
         os.system(cmd)
         cmd = "evselect table='"+datapath+det+filt_label+"_filt.fits:EVENTS' destruct=false withfilteredset=true withimageset=true imageset="+datapath+det+filt_label+"_detmap.ds xcolumn=DETX ycolumn=DETY #withxranges=true ximagemin=-1500 ximagemax=1500 withyranges=true #yimagemin=-1500 yimagemax=1500 imagebinning='imageSize' #ximagesize=20 yimagesize=20 #writedss=true updateexposure=true"
         print(" ")
         print("3 create det map")
         print(cmd)
         os.system(cmd)

         if det == "pn":
             cmd = "evselect table='"+datapath+det+filt_label+"_filt.fits:EVENTS' withspectrumset=yes " \
                                                              "spectrumset="+datapath+det+"_"+srcName+".pha energycolumn=PI spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=20479 expression='#XMMEA_EP && (PATTERN<=4) && ((X,Y) IN "+srcReg+")'"
         else:
             cmd = "evselect table='"+datapath+det+filt_label+"_filt.fits:EVENTS' withspectrumset=yes " \
                                                              "spectrumset="+datapath+det+"_"+srcName+".pha energycolumn=PI spectralbinsize=15 withspecranges=yes specchannelmin=0 specchannelmax=11999 expression='#XMMEA_EM && (PATTERN<=12) && ((X,Y) IN "+srcReg+")'"
         print(" ")
         print("4 extract source spectrum")
         print(cmd)
         os.system(cmd)
         if det == "pn":
             cmd = "evselect table='"+datapath+det+filt_label+"_filt.fits:EVENTS' withspectrumset=yes " \
                                                              "spectrumset="+datapath+det+"_BKG_for"+srcName+".pha energycolumn=PI spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=20479 expression='#XMMEA_EP && (PATTERN<=4) && ((X,Y) IN "+bkgReg+")'"
         else:
             cmd = "evselect table='"+datapath+det+filt_label+"_filt.fits:EVENTS' withspectrumset=yes " \
                                                              "spectrumset="+datapath+det+"_BKG_for"+srcName+".pha energycolumn=PI spectralbinsize=15 withspecranges=yes specchannelmin=0 specchannelmax=11999 expression='#XMMEA_EM && (PATTERN<=12) && ((X,Y) IN "+bkgReg+")'"
         print(" ")
         print("5 extract background spectrum")
         print(cmd)
         os.system(cmd)

         cmd = "backscale spectrumset="+datapath+det+"_"+srcName+".pha badpixlocation="+datapath+det+filt_label+"_filt.fits"
         print(" ")
         print("6 create source backscale keyword")
         print(cmd)
         os.system(cmd)
         cmd = "backscale spectrumset="+datapath+det+"_BKG_for"+srcName+".pha badpixlocation="+datapath+det+filt_label+"_filt.fits"
         print(" ")
         print("7 create background backscale keyword")
         print(cmd)
         os.system(cmd)

         cmd = "rmfgen spectrumset="+datapath+det+"_"+srcName+".pha rmfset="+datapath+det+"_"+srcName+".rmf"
         print(" ")
         print("8 create rmf")
         print(cmd)
         os.system(cmd)
         cmd = "arfgen spectrumset="+datapath+det+"_"+srcName+".pha arfset="+datapath+det+"_"+srcName+".arf withrmfset=yes " \
                                                                                                      "rmfset="+datapath+det+"_"+srcName+".rmf badpixlocation="+datapath+det+filt_label+"_filt.fits extendedsource=yes detmaptype=dataset detmaparray="+datapath+det+filt_label+"_detmap.ds"
         print(" ")
         print("9 create arf")
         print(cmd)
         os.system(cmd)
         cmd="fparkey " +det+"_BKG_for"+srcName+".pha " +datapath+det+"_"+srcName+".pha " +"BACKFILE add=yes"
         print("10 add key")
         print(" ")
         print(cmd)
         os.system(cmd)
         cmd="fparkey "+det+"_"+srcName+".rmf " +datapath+det+"_"+srcName+".pha " +"RESPFILE add=yes"
         print(" ")
         print(cmd)
         os.system(cmd)
         cmd="fparkey " +det+"_"+srcName+".arf " + datapath + det + "_" + srcName + ".pha " + "ANCRFILE add=yes"
         print(" ")
         print(cmd)
         os.system(cmd)
         print(" ")
         print(" ")



   if combine_spec:
      datapath = path + "/" + obsID + "/cal/"
      os.chdir(datapath+'add_2obs_spec/')
      srcName='vcc1154_polygon'
      # cmd="epicspeccombine pha="+'"mos1_{0}.pha mos2_{0}.pha mos1_{0}_02.pha mos2_{0}_02.pha"'.format(srcName) \
      #         +" bkg="+'"mos2_BKG_for{0}.pha mos2_BKG_for{0}.pha mos1_BKG_for{0}_02.pha mos2_BKG_for{0}_02.pha"'.format(srcName) \
      #         +" rmf="+'"mos1_{0}.rmf mos2_{0}.rmf mos1_{0}_02.rmf mos2_{0}_02.rmf"'.format(srcName) \
      #         +" arf="+'"mos1_{0}.arf mos2_{0}.arf mos1_{0}_02.arf mos2_{0}_02.arf"'.format(srcName) \
      #         +" filepha="+'"src_spec_grp_mos_2obs.pha"'\
      #         +" filebkg="+'"bkg_spec_grp_mos_2obs.pha"'\
      #         +" filersp="+'"response_grp_mos_2obs.rmf"'

      cmd="epicspeccombine pha="+'"mos1_{0}.pha mos2_{0}.pha mos1_{0}_obs2.pha mos2_{0}_obs2.pha "'.format(srcName) \
              +" bkg="+'"mos1_BKG_for{0}.pha mos2_BKG_for{0}.pha mos2_BKG_for{0}_obs2.pha mos2_BKG_for{0}_obs2.pha "'.format(srcName) \
              +" rmf="+'"mos1_{0}.rmf mos2_{0}.rmf mos1_{0}_obs2.rmf mos2_{0}_obs2.rmf"'.format(srcName) \
              +" arf="+'"mos1_{0}.arf mos2_{0}.arf mos1_{0}_obs2.arf mos2_{0}_obs2.arf"'.format(srcName) \
              +" filepha="+'"merge_2obs_mos_src.pha"'\
              +" filebkg="+'"merge_2obs_mos_bkg.pha"'\
              +" filersp="+'"response_merge_2obs_mos.rmf"'

      # cmd="epicspeccombine pha="+'"pn_{0}.pha pn_{0}_obs2.pha "'.format(srcName) \
      #         +" bkg="+'"pn_BKG_for{0}.pha pn_BKG_for{0}_obs2.pha "'.format(srcName) \
      #         +" rmf="+'"pn_{0}.rmf pn_{0}_obs2.rmf"'.format(srcName) \
      #         +" arf="+'"pn_{0}.arf pn_{0}_obs2.arf"'.format(srcName) \
      #         +" filepha="+'"merge_2obs_pn_src.pha"'\
      #         +" filebkg="+'"merge_2obs_pn_bkg.pha"'\
      #         +" filersp="+'"response_merge_2obs_pn.rmf"'

      print(cmd)
      os.system(cmd)
      cmd = "fparkey " + "response_merge_2obs_mos.rmf " + "merge_2obs_mos_src.pha" + " RESPFILE add=yes"
      print(cmd)
      os.system(cmd)
      cmd = "fparkey " + "merge_2obs_mos_bkg.pha " + "merge_2obs_mos_src.pha" + " BACKFILE add=yes"
      print(cmd)
      os.system(cmd)
      cmd = "fparkey " + "None " + "merge_2obs_mos_src.pha" + " ANCRFILE add=yes"
      print(cmd)
      os.system(cmd)

   if merge_evt:
      datapath = path + "/" + obsID + "/cal/"
      os.chdir(datapath)

      # cmd='merge set1=mos1_filt_time.fits set2=mos2_filt_time.fits outset=mos_all_filt_time.fits withradec=Y ra=53.1476 dec=-27.7368'
      # print(cmd)
      # os.system(cmd)
      #
      # cmd='merge set1=mos_all_filt_time.fits set2=pn_filt_time.fits outset={0}_merge_filt_time.fits withradec=Y ra=53.1476 dec=-27.7368'.format(obsID)
      # print(cmd)
      # os.system(cmd)

      # cmd='cp {0}_merge_filt_time.fits /Volumes/pulsar/xmmCDFS/merge_evt/{0}_merge_filt_time.fits'.format(obsID)
      # cmd = 'cp pn_filt_time.fits /Volumes/pulsar/xmmCDFS/merge_evt/{0}_pn_filt_time.fits'.format(obsID)
      # print(cmd)
      # os.system(cmd)

      cmd = 'cp mos_all_filt_time.fits /Volumes/pulsar/xmmCDFS/merge_evt/{0}_mos_filt_time.fits'.format(obsID)
      print(cmd)
      os.system(cmd)
# ##====merge different obs=======##

epoch=[0,2,8,12,19,23,33]
i=0
while i in range(len(epoch)):
   EP=i+1
   #epoch 标号
   j=epoch[i]
   #第i个观测号
   while j <epoch[i+1]-1:
      obsID1=obsList[j]
      obsID2=obsList[j+1]

      os.chdir(path+'/merge_evt/')
      if j==epoch[i]:
         # cmd='cp {0}_merge_filt_time.fits ep{1}_temp{2}_merge_filt_time.fits'.format(obsID1,EP,j)
         cmd = 'cp {0}_mos_filt_time.fits ep{1}_temp{2}_mos_filt_time.fits'.format(obsID1, EP, j)

         print(cmd)
         os.system(cmd)
      ##先把第一个文件cp为temp_merge_filt_time.fits

      cmd='merge set1=ep{0}_temp{1}_mos_filt_time.fits set2={2}_mos_filt_time.fits outset=ep{0}_temp{3}_mos_filt_time.fits withradec=Y mergedifferentobs=yes ra=53.1476 dec=-27.7368'.format(EP,j,obsID2,j+1)
      print(cmd)
      os.system(cmd)
      j+=1
   i+=1
