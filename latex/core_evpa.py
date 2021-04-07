#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 14:05:06 2020

@author: stary
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def RM(x,a,b):
    return a+b*x*x

# core sample
#in deg.    
pa1=np.array([12.42,6.46,0.06,0.86]) 
pa2=np.array([9.20,3.70,1.04,1.18])
pa3=np.array([11.83,5.81,1.57,2.05])
pa4=np.array([14.56,7.92,1.00,1.59])
pa5=np.array([13.39,7.38,1.43,0.25])
pa=np.array([pa1,pa2,pa3,pa4,pa5])
#pcov=np.zeros(5,5)
popt=np.zeros((5,2))

X=np.array([8.6268e10,8.8268e10,9.8328e10,1.00268e11])  #freq.
x=3e8/X #wavelength=lambda
ydata=np.zeros((5,50))
xdata=np.linspace(x[0],x[3], 50)

for i in range(5):
    #print(i)
    Y=pa[i]
    ycore=Y*np.pi/180 #in rad
    popt[i], pcov = curve_fit(RM, x,ycore)#,sigma=yerr*np.pi/180, absolute_sigma=True) #absolute true 绝对sigma值;Default method is ‘lm’ for unconstrained 
    ydata[i]=np.array([RM(j, popt[i][0],popt[i][1]) for j in xdata])
    print("PA0 & RM:", popt[i])
    #paramater err
    perr = np.sqrt(np.diag(pcov))
    print("perr:",perr)

plt.figure()
flag=['darkorange','g','black','royalblue','red']
id=[r'1: $6.0 \times 10^{4}$',r'2: $3.7 \times 10^{4}$',r'3: $4.8 \times 10^{4}$',r'4: $6.7 \times 10^{4}$',r'5: $6.4 \times 10^{4}$']
for i in range(5):
    plt.plot(1e6*xdata**2,ydata[i]*180/np.pi,'--',color=flag[i]) #x xais in mm
    plt.plot(1e6*x**2,pa[i],'o', color=flag[i], label=id[i],markersize=5)
    plt.xlabel(r' $\rm Wavelength^{2}$ ($\rm mm^{2}$)')
    plt.ylabel(r'EVPA (deg)')
    plt.title('RM')
    #plt.savefig('core-sample.eps'）
plt.legend(loc='upper left')
plt.savefig('core-sample.eps')


plt.figure()
plt.subplot(321)
plt.plot(1e6*xdata**2,ydata[0]*180/np.pi,'--') #x xais in mm
plt.plot(1e6*x**2,pa[0],'o', markersize=5)
#plt.tick_params(axis='x', labelsize=)
#plt.xlabel(r' $\rm Wavelength^{2}$ ($\rm mm^{2}$)')
plt.ylabel(r'EVPA (deg)')
plt.text(9,10,r'1')
#plt.savefig('core-sample-1.eps')

plt.subplot(322)
plt.plot(1e6*xdata**2,ydata[1]*180/np.pi,'--') #x xais in mm
plt.plot(1e6*x**2,pa[1],'o', markersize=5)
plt.xlabel(r' $\rm Wavelength^{2}$ ($\rm mm^{2}$)')
#plt.ylabel(r'EVPA (deg)')
plt.text(9,7.5,r'2')
#plt.savefig('core-sample-2.eps')

plt.subplot(323)
plt.plot(1e6*xdata**2,ydata[3]*180/np.pi,'--') #x xais in mm
plt.plot(1e6*x**2,pa[3],'o', markersize=5)
plt.xlabel(r' $\rm Wavelength^{2}$ ($\rm mm^{2}$)')
plt.ylabel(r'EVPA (deg)')
plt.text(9,12,r'4')
#plt.savefig('core-sample-4.eps')

plt.subplot(324)
plt.plot(1e6*xdata**2,ydata[4]*180/np.pi,'--') #x xais in mm
plt.plot(1e6*x**2,pa[4],'o', markersize=5)
plt.xlabel(r' $\rm Wavelength^{2}$ ($\rm mm^{2}$)')
#plt.ylabel(r'EVPA (deg)')
plt.text(9,11,r'5')
#plt.savefig('core-sample-4.eps')

plt.subplot(313)
plt.plot(1e6*xdata**2,ydata[2]*180/np.pi,'--') #x xais in mm
plt.plot(1e6*x**2,pa[2],'o', markersize=5)
plt.xlabel(r' $\rm Wavelength^{2}$ ($\rm mm^{2}$)')
plt.ylabel(r'EVPA (deg)')
plt.text(9,10,r'3')
plt.savefig('core-sample-3.eps')