import numpy as np
from astropy.wcs import WCS

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

def make_phy_reg(srcid,ra,dec,psfradii,imgfile,inpath,outpath):
    ## only for chandra now
    w = WCS(inpath + imgfile)
    src_x, src_y = w.all_world2pix(ra, dec, 1)
    phy_x = src_x + 2896
    phy_y = src_y + 2896
    src_radius = psfradii*2.03252
    for i in range(len(srcid)):
        with open(outpath + '{0}.reg'.format(i + 1), 'w+') as f1:
            reg = 'circle(' + str(phy_x[i]) + ',' + str(phy_y[i]) + ',' + str(src_radius[i]) + ')'
            f1.writelines(reg)
        with open(outpath + 'all.reg', 'a+') as f2:
            f2.writelines(reg + '# text=' + '{' + str(srcid[i]) + '}' '\n')
