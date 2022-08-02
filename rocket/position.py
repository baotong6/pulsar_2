#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
import math
from dataclasses import dataclass
from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import SkyCoord

@dataclass
class Circle:
    x: float
    y: float
    r: float

    @property
    def coord(self):
        return self.x, self.y

def select_src_bypos(srcID_list,ra,dec,ra_c,dec_c,inter_radius,outpath=None,outname=None,save=0):
    # path = '/Volumes/pulsar/47Tuc/merge_data/timing/txt_startover/txt_all_obs_p50/'
    c1 = SkyCoord(ra * u.deg, dec * u.deg, frame='fk5')
    c2 = SkyCoord(ra_c * u.deg, dec_c * u.deg, frame='fk5')
    dist = c1.separation(c2)
    dist = dist.arcsec
    inter_srcID = srcID_list[np.where(dist < inter_radius)[0]]
    counts=[]
    for i in range(len(inter_srcID)):
        evt=np.loadtxt(outpath+f'{inter_srcID[i]}.txt')
        counts.append(len(evt))
    print(counts)
    if save:
        if not(outpath) or not(outname):
            print('please provide outpath and outname if save')
            return inter_srcID
        outevt=np.column_stack((inter_srcID,counts,dist[np.where(dist < inter_radius)[0]]))
        np.savetxt(outpath + '{0}'.format(outname), outevt, fmt='%10d  %10d %10.5f')
    return inter_srcID

def where_region_annulus(x, y, reg):
    an=reg
    dist = np.sqrt((x-reg[0])**2+(y-reg[1])**2)
    temp_in= dist - an[2]
    temp_out=dist-an[3]
    return np.where((temp_in >= 0)&(temp_out <= 0))


def find_intersection(c1: Circle, c2: Circle)-> float:
    """Finds intersection area of two circles.

    Returns intersection area of two circles otherwise 0
    """
    d = math.dist(c1.coord, c2.coord)
    rad1sqr = c1.r ** 2
    rad2sqr = c2.r ** 2

    if d == 0:
        # the circle centers are the same
        return math.pi * min(c1.r, c2.r) ** 2

    angle1 = (rad1sqr + d ** 2 - rad2sqr) / (2 * c1.r * d)
    angle2 = (rad2sqr + d ** 2 - rad1sqr) / (2 * c2.r * d)

    # check if the circles are overlapping
    if (-1 <= angle1 < 1) or (-1 <= angle2 < 1):
        theta1 = math.acos(angle1) * 2
        theta2 = math.acos(angle2) * 2

        area1 = (0.5 * theta2 * rad2sqr) - (0.5 * rad2sqr * math.sin(theta2))
        area2 = (0.5 * theta1 * rad1sqr) - (0.5 * rad1sqr * math.sin(theta1))

        return area1 + area2
    elif angle1 < -1 or angle2 < -1:
        # Smaller circle is completely inside the largest circle.
        # Intersection area will be area of smaller circle
        # return area(c1_r), area(c2_r)
        return math.pi * min(c1.r, c2.r) ** 2
    return 0