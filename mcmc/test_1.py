import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
from scipy.optimize import curve_fit
def judge(xy,ab,ij):
    l1=np.sqrt((xy[0]-ij[0])**2+(xy[1]-ij[1])**2)
    l2= np.sqrt((ab[0] - ij[0]) ** 2 + (ab[1] - ij[1]) ** 2)
    if l1>l2:
        return [0,1]
    elif l1==l2:
        return [0,0]
    elif l1<l2:
        return [1,0]

def get_all(xy,ab):
    r1=r2=0
    for i in range(3):
        for j in range(3):
            ij=[i,j]
            r1+=judge(xy,ab,ij)[0]
            r2+=judge(xy,ab,ij)[1]

    return (r1,r2)




x=np.linspace(0,2,1000)
y=np.linspace(0,2,1000)
P=np.zeros(1000)
#P=np.zeros([100,100])
for j in range(1000):
    # for k in range(100):
        out_P=0
        xy=[x[j],y[j]]
        for i in range(100):
            ab=[np.random.random(1)[0]*2,np.random.random(1)[0]*2]
            (r1,r2)=get_all(xy,ab)
            #print(r1,r2)
            if r1>r2:
                out_P+=1
        #P[j][k]=out_P/100.
        P[j]=out_P/100.
x=x-1
def f(x):
    if np.abs(x)<=1:
        return 1-np.pi/4*x**2
    else:
        return 1-np.pi/4*x**2+x**2*np.arccos(1/np.abs(x))-np.sqrt(x**2-1)

x1=np.linspace(-1.414,1.414,1000)
y1=[f(x1[i]) for i in range(len(x1))]
# n = len(x)                          #the number of data
# mean = sum(x*P)/n                   #note this correction
# sigma = sum(P*(x-mean)**2)/n        #note this correction
#
# def gaus(x,a,x0,sigma):
#        return a*np.exp(-(x-x0)**2/(2*sigma**2))
#
# popt,pcov = curve_fit(gaus,x,P,p0=[1,mean,sigma])
# print(popt,pcov)
plt.plot(x*1.414,P,'b+:',label='data')
plt.plot(x1,y1,'ro:',label='fit')
#plt.legend()
plt.show()
# plt.contourf(x,y,P)
# plt.colorbar()