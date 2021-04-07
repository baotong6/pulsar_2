import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import string

def get_and_plot():
    srcID='202'
    obsID=['5934','6362','6365','9500','9501','9502',
           '9503','9504','9505','9854','9855','9892','9893']
    for id in obsID:
        path='/Volumes/pulsar/LimWin_damage/'
        os.chdir(path)
        # regfile='region_{0}/region_90/{1}.reg'.format(id,srcID)
        # with open(path+'merge_data/timing/'+regfile,'r') as f:
        #     reg=f.readline()
        reg='circle(17:51:38.3376,-29:42:40.032,7.63)'
        os.chdir(path + id + '/cal')

        asolname = 'aspect_1_bcc.fits'
        cmd='dither_region infile={0} region="{1}" maskfile=maskfile.fits ' \
            'wcsfile=evt2file_bcc.fits outfile=fracarea_vs_time_{2}_src{3}.fits clobber=yes verbose=5'.format(asolname,reg,id,srcID)
        print(cmd)
        #os.system(cmd)
        filename = 'fracarea_vs_time_{0}_src{1}.fits'.format(id,srcID)
        hdul = fits.open(path+id+'/cal/' + filename)
        time = hdul[1].data.field(0)
        time = time - time[0]
        frac = hdul[1].data.field(4)
        plt.figure()
        plt.plot(time, frac)
        plt.savefig('frac_{0}_src{1}.eps'.format(id,srcID))
        plt.show()
#get_and_plot()
def compute_index():
    path = '/Volumes/pulsar/LimWin_damage/'
    srcID=['324','118','16','206','40','20','88'
        ,'141','521','229','42','513'
        ,'371','196','152','194','99','191'
        ,'500','153','114','244','202']
    #srcID=['324']
    obsID=['5934','6362','6365','9500','9501','9502',
           '9503','9504','9505','9854','9855','9892','9893']
    dither_info_time=[]
    dither_info_frac=[]
    index_total = []
    for i in range(len(srcID)):
        index_single_obs=[]
        obstime=[]
        for id in obsID:
            os.chdir(path + id + '/cal')
            filename = 'fracarea_vs_time_{0}_src{1}.fits'.format(id, srcID[i])
            hdul = fits.open(path + id + '/cal/' + filename)
            time = hdul[1].data.field(0)
            frac = hdul[1].data.field(4)

            time=np.array(time)
            frac=np.array(frac)
            time.astype('float')
            frac.astype('float')

            dither_info_time=np.concatenate((dither_info_time,time))
            dither_info_frac=np.concatenate((dither_info_frac,frac))

            time_2=np.concatenate((time[1:],[0]))
            exptime=time_2[0:-1]-time[0:-1]

            frac_2=np.concatenate((frac[1:],[0]))
            exp_frac=frac_2[0:-1]/2+frac[0:-1]/2
            index_single_obs.append(np.sum(exp_frac*exptime)/np.sum(exptime))
            obstime.append(np.sum(exptime))

        temp_index=np.sum(np.array(obstime)*np.array(index_single_obs))/np.sum(np.array(obstime))
        index_total.append(temp_index)

        dither_info=np.column_stack((dither_info_time,dither_info_frac))
        #print(dither_info)
        np.savetxt(path+'merge_data/dither/'+'{0}_info.txt'.format(srcID[i]),dither_info,fmt='%10.5f  %10.5f')

    print(index_total)

compute_index()


# path='/Volumes/pulsar/SgrA/3393/cal/'
# os.chdir(path)
# filename='fracarea_vs_time_3393.fits'
# hdul=fits.open(path+filename)
# time=hdul[1].data.field(0)
# time=time-time[0]
# frac=hdul[1].data.field(4)
# plt.plot(time,frac)
# plt.show()
