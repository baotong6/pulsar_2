import numpy as np
from astropy.wcs import WCS
import os
from astropy.io import fits

def read_region(regname):
    reg_file=[]
    with open(regname, 'r') as file_to_read:
        while True:
            lines = file_to_read.readline() # 整行读取数据
            reg_file.append(lines)
            if not lines:
                break
                pass
    region=reg_file[-2][7:-2]
    reg_x,reg_y,reg_r=[float(i) for i in region.split(',')]
    return [reg_x,reg_y,reg_r]

def make_fk5_reg(srcid,ra,dec,psfradii,outpath):
    ## psfradii must be arcsec
    for i in range(len(srcid)):
        with open(outpath+'{0}.reg'.format(srcid[i]), 'w+') as f1:
            f1.writelines('fk5' + '\n')
            reg = 'circle(' + str(ra) + ',' + str(dec) + ',' + str(psfradii) + '"' + ')'
            f1.writelines(reg)
        with open(outpath+'all.reg', 'a+') as f2:
            reg = 'circle(' + str(ra) + ',' + str(dec) + ',' + str(psfradii) + '"' + ')' + '# text=' + '{' + str(
                srcid[i]) + '}'
            f2.writelines(reg + '\n')

def make_phy_reg(srcid,x,y,psfradii,outpath,inpath=None,imgfile=None,coordtype='physical',singlename=None):
    ## coordtype=['physical','WCS']
    ## only for chandra now
    if coordtype=='WCS':
        if not os.path.exists(inpath + imgfile):
            print('ERROR: Missing or Could not find inpath+imgfile in WCS coord mode')
            return None

        w = WCS(inpath + imgfile)
        src_x, src_y = w.all_world2pix(x, y, 1)
        binx, biny = np.shape(fits.open(inpath + imgfile)[0].data)
        phy_x = src_x + 4096 - binx / 2
        phy_y = src_y + 4096 - biny / 2
        src_radius = psfradii*2.03252

        reg_x=phy_x;reg_y=phy_y;reg_r=src_radius

    elif coordtype=='physical':
        reg_x = x;reg_y = y;reg_r = psfradii
        # print(reg_r,reg_x,reg_y)
        if singlename and len(srcid)==1:
            for i in range(len(x)):
                with open(outpath+'{0}.reg'.format(singlename), 'w+') as f1:
                    reg = 'circle(' + str(reg_x[i]) + ',' + str(reg_y[i]) + ',' + str(reg_r[i]) + ')'
                    f1.writelines(reg)
            return None
        os.system(f'rm {outpath}all.reg')
        os.system(f'rm {outpath}all_bkg.reg')
        for i in range(len(srcid)):
            with open(outpath + '{0}.reg'.format(i + 1), 'w+') as f1:
                reg = 'circle(' + str(reg_x[i]) + ',' + str(reg_y[i]) + ',' + str(reg_r[i]) + ')'
                f1.writelines(reg)
            with open(outpath+'{0}_bkg.reg'.format(i + 1), 'w+') as f2:
                bkgreg='annulus(' + str(reg_x[i]) + ',' + str(reg_y[i]) + ',' + str(2*reg_r[i]) +',' + str(4*reg_r[i])+ ')'
                f2.writelines(bkgreg)
            with open(outpath + 'all.reg', 'a+') as f3:
                f3.writelines(reg + '# text=' + '{' + str(srcid[i]) + '}' '\n')
            with open(outpath + 'all_bkg.reg', 'a+') as f3:
                f3.writelines(bkgreg + '# text=' + '{' + str(srcid[i]) + '}' '\n')

    else:
        print('EROOR: Unkown/Ungiven coordtype')
        return None