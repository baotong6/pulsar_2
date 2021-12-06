#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u

def select_src_bypos(srcID_list,ra,dec,ra_c,dec_c,inter_radius,outpath,outname):
    path = '/Volumes/pulsar/47Tuc/merge_data/timing/txt_startover/txt_all_obs_p50/'
    c1 = SkyCoord(ra * u.deg, dec * u.deg, frame='fk5')
    c2 = SkyCoord(ra_c * u.deg, dec_c * u.deg, frame='fk5')
    dist = c1.separation(c2)
    dist = dist.arcsec
    inter_srcID = srcID_list[np.where(dist < inter_radius)[0]]
    np.savetxt(outpath + '{0}'.format(outname), inter_srcID, fmt='%10d')
