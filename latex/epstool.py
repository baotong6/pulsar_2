#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import sys,os
import read_csv as data
path='/Users/baotong/Desktop/aas/V63/figure/LW'
ID=data.ID_LW
for i in range(len(ID)):
    os.chdir(path)
    file=str(ID[i])+'_lc.eps'
    file_out=str(ID[i])+'_lc_cut.eps'
    command='epstool --copy --bbox '+file +' '+file_out
    print(command)
    os.system(command)

