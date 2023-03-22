#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 19:33:00 2023

@author: helin
"""

import os 
import glob
import subprocess
import csv
import pandas as pd
import math
BEHR=os.path.realpath("/Users/helin/Desktop/BEHR/BEHR")

behr_path = '/Users/helin/Documents/test/chang/HR/'
os.chdir(behr_path)

file = open('/Users/helin/Documents/test/chang/93_BEHR_par.csv')
fileReader = csv.reader(file)
filedata = list(fileReader)

list3=['name','HR','HRl','HRu']
list4=[]
         
for i in filedata[1:]:
    if float(i[filedata[0].index('softeff')]) == 0.0:
    
        softsrc=str(round(float(i[filedata[0].index('softsrc')])))
        hardsrc=str(round(float(i[filedata[0].index('hardsrc')])))
        softbkg=str(round(float(i[filedata[0].index('softbkg')])))
        hardbkg=str(round(float(i[filedata[0].index('hardbkg')])))
        softarea=str(round(float(i[filedata[0].index('softarea')])))
        hardarea=str(round(float(i[filedata[0].index('hardarea')])))
        
        ### if effective area correction is needed! ###
        ###=========================================###
        # softeff=str(round(float(i[filedata[0].index('softeff')])))
        # hardeff=str(round(float(i[filedata[0].index('hardeff')])))
           
        p1="softsrc="+softsrc
        p2="hardsrc="+hardsrc
        p3="softbkg="+softbkg
        p4="hardbkg="+hardbkg
        p5="softarea="+softarea
        p6="hardarea="+hardarea
        # p9="softeff="+softeff
        # p10="hardeff="+hardeff
        filename=i[filedata[0].index('name')]+'_93'
        p7="output="+filename
        p8="outputHR=True"
        p11="level=68.0"
        subprocess.call([BEHR,p1,p2,p3,p4,p5,p6,p7,p8,p11])
        # subprocess.call([BEHR,p1,p2,p3,p4,p5,p6,p9,p10,p7,p8,p11])
        
        with open(behr_path+filename+'.txt',"r") as f:
            data=f.readlines()
            HR=data[2].split('\t')[2]
            HR_lo=data[2].split('\t')[4]
            HR_up=data[2].split('\t')[5]
            list4.append([filename,HR,HR_lo,HR_up])
            print(filename+':'+HR+' '+HR_lo+' '+HR_up)
                   
test = pd.DataFrame(columns=list3,data=list4)
test.to_csv('/Users/helin/Documents/test/chang/93_HR_sup.csv',index=False)
