#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import string

path='/Users/baotong/ciao-4.10/BEHR/idl/'
input_1=np.loadtxt(path+'HR_1671_info.txt')
input_2=np.loadtxt(path+'HR_1671_info_com.txt')
net_soft_1=input_1[:,0]-input_1[:,1]/input_1[:,4]
net_hard_1=input_1[:,2]-input_1[:,3]/input_1[:,5]
net_soft_2=input_2[:,0]-input_2[:,1]/input_2[:,4]
net_hard_2=input_2[:,2]-input_2[:,3]/input_2[:,5]

plt.scatter(net_soft_1,net_hard_1,color='red')
plt.scatter(net_soft_2,net_hard_2,color='green')
plt.legend(['before','norm'])
plt.xlabel('soft')
plt.ylabel('hard')
plt.show()