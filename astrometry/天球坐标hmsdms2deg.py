"""

Author: HU Zhensong
该文档用来把string形式的hms dms转成deg

"""
import re
def HMS2deg(hra):
    rera = re.findall('[-+]?\d*\.\d+|[-+]?\d+', hra)
    if (float(rera[0])<0):
        if (len(rera)==3):
            ra = float(rera[0])*15 -  float(rera[1])*15/60 -  float(rera[2])*15/3600
        if (len(rera)==2):
            ra = float(rera[0])*15 -  float(rera[1])*15/60
    else:
        if (len(rera)==3):
            ra = float(rera[0])*15 +  float(rera[1])*15/60 + float(rera[2])*15/3600
        if (len(rera)==2):
            ra = float(rera[0])*15 +  float(rera[1])*15/60
    return ra
def DMS2deg(ddec):
    redec = re.findall('[-+]?\d*\.\d+|[-+]?\d+', ddec)  
    if (float(redec[0]) < 0):
        if (len(redec) == 3):
            dec = float(redec[0]) -  float(redec[1])/60 -  float(redec[2])/3600
        if (len(redec) == 2):
            dec = float(redec[0]) -  float(redec[1])/60
    else:
        if (len(redec) == 3):
            dec = float(redec[0]) +  float(redec[1])/60 +  float(redec[2])/3600
        if (len(redec) == 2):
            dec = float(redec[0]) +  float(redec[1])/60

    return dec