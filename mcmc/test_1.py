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

def baoyan(N=48,m=12,n=4,x=20):
    ## N为总人数；m为宿舍数；n为每个宿舍的人数；x为保研总人数
    NUM=10000
    count=0
    for k in range(NUM):
        a=np.arange(0,48,1)
        b=np.random.choice(a,20,replace=False,p=None)
        c=np.sort(b)
        i=0
        while i <= len(c)-4:
            if c[i]%4==0 and c[i+1]==c[i]+1 and c[i+2]==c[i]+2 and c[i+3]==c[i]+3:
                count+=1
                break
            else:i += 1
    print(count/NUM)
# baoyan()

def pingpang(a=0.9,n=11,N=5000):
    ## a为胜率，n为要求连胜回合数
    clist=[]
    for i in range(N):
        count=1
        c=np.zeros(11)
        while True:
            c=np.concatenate((c[1:],np.random.choice([0,1],1,p=[0.1,0.9])))
            if np.sum(c)==11:
                break
            else:
                count+=1
        clist.append(count)
    d=np.average(np.array(clist))
    plt.hist(clist,bins=np.linspace(11,100,90),histtype='step',density=True)
    plt.show()
    print(d)
    return d
pingpang()
# x=np.linspace(0,2,1000)
# y=np.linspace(0,2,1000)
# P=np.zeros(1000)
# #P=np.zeros([100,100])
# for j in range(1000):
#     # for k in range(100):
#         out_P=0
#         xy=[x[j],y[j]]
#         for i in range(100):
#             ab=[np.random.random(1)[0]*2,np.random.random(1)[0]*2]
#             (r1,r2)=get_all(xy,ab)
#             #print(r1,r2)
#             if r1>r2:
#                 out_P+=1
#         #P[j][k]=out_P/100.
#         P[j]=out_P/100.
# x=x-1
# def f(x):
#     if np.abs(x)<=1:
#         return 1-np.pi/4*x**2
#     else:
#         return 1-np.pi/4*x**2+x**2*np.arccos(1/np.abs(x))-np.sqrt(x**2-1)
#
# x1=np.linspace(-1.414,1.414,1000)
# y1=[f(x1[i]) for i in range(len(x1))]
# # n = len(x)                          #the number of data
# # mean = sum(x*P)/n                   #note this correction
# # sigma = sum(P*(x-mean)**2)/n        #note this correction
# #
# # def gaus(x,a,x0,sigma):
# #        return a*np.exp(-(x-x0)**2/(2*sigma**2))
# #
# # popt,pcov = curve_fit(gaus,x,P,p0=[1,mean,sigma])
# # print(popt,pcov)
# plt.plot(x*1.414,P,'b+:',label='data')
# plt.plot(x1,y1,'ro:',label='fit')
# #plt.legend()
# plt.show()
# # plt.contourf(x,y,P)
# # plt.colorbar()