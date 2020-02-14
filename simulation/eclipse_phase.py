import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import string
import pylab as pylab
import datetime
from scipy.interpolate import lagrange
from scipy import interpolate
from scipy.fftpack import fft,ifft
import scipy.signal as ss
import scipy.stats as stats
import sympy
from scipy.optimize import root
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import numpy as np
# from sympy.abc import x,y
#note 1为伴星，2为白矮星
# P=13581.
# r1=0.368*6.5e5*1e3
# a=1.241*6.5e5*1e3
# M1=0.287
# r2=7.8*1e6*((1.44/M1)**0.6667-(M1/1.44)**0.6667)**0.5


# r1=0.368*6.5e5*1e3
# M1=np.linspace(0.01,1,1000)
# M=M1+M2
#a=(P/31536000)**(2./3.)*(1.5e8*1e3)/(M**0.33)
#r1=M1**0.7*6.5e5*1e3
#a=(P/31536000)**(2./3.)*(1.5e8*1e3)/((M2+(r1/(6.5e5*1e3))**(10/7))**0.33)

plt.semilogy()
# i=np.linspace(74,90,1000)/180*np.pi
# i=75/180.*np.pi
x,y=sympy.symbols('x y')
M2=0.75
P=44365.57232
r2=7.8*1e6*((1.44/M2)**0.6667-(M2/1.44)**0.6667)**0.5/(6.5e5*1e3)
# dt12_model=0.04*P
# dt23_model=0.04*P
r1=np.linspace(0.1,1,500)
i=np.linspace(80,90,500)/180*np.pi
dt12=P/(2*np.pi)*((((r1+r2)/(1.31163/ (((r1**(10/7.))+M2) ** 0.33)))**2-np.cos(i)**2)-(((r1-r2)/(1.31163/ (((r1**(10/7.))+M2) ** 0.33)))**2-np.cos(i)**2))
dt23=2*P/(2*np.pi)*(((r1-r2)/(1.31163/ (((r1**(10/7.))+M2) ** 0.33)))**2-np.cos(i)**2)**0.5
where_are_nan = np.isnan(dt12)
dt12[where_are_nan] = 0.
where_are_nan = np.isnan(dt23)
dt23[where_are_nan] = 0.

dt12_model=np.zeros(len(dt12))+0.03*P
dt23_model=np.zeros(len(dt12))+0.4*P
r1, i = np.meshgrid(r1, i)
ax = plt.axes(projection='3d')
ax.plot3D(r1, i, np.log10(np.abs(dt12-dt12_model)+1.), 'red')
ax.plot3D(r1, i, np.log10(np.abs(dt23-dt23_model)+1.), 'green')
ax.plot3D(r1, i, np.log10(dt12_model), 'yellow')
ax.plot3D(r1, i, np.log10(dt23_model), 'yellow')
print(np.log10(np.abs(dt23-dt23_model)+1.))
plt.show()





#
# print(M2,P,r2,dt12_model,dt23_model)
# def my_func(X):
#     M1,r1,a,i=X
#
#     f=[2*P/(2*np.pi)*(((r1-r2)/a)**2-np.cos(i)**2)**0.5-dt23_model,
#        P/(2*np.pi)*((((r1+r2)/a)**2-np.cos(i)**2)**0.5-(((r1-r2)/a)**2-np.cos(i)**2)**0.5)-dt12_model,
#        a-1.31163/ ((M1+M2) ** 0.33),
#     r1-M1**0.7]
#
#     return f
# M1=np.linspace(0.01,1,1000)
#
#
# a=(P / 31536000) ** (2. / 3.) * (230.) / ((M1+M2) ** 0.33)
# r1=M1**0.7
# #print(a)
# print(r1)
#sol2 = root(my_func, [0.5,0.2,20.,1.14])
#print(sol2)

# aa=sympy.solve([2*P/(2*np.pi)*(((x-r2)/((P/31536000)**(2./3.)*(1.5e8*1e3)/((M2+(x/(6.5e5*1e3))**(10/7))**0.33)))**2-sympy.cos(y)**2)**0.5-dt23_model],
#          [P/(2*np.pi)*((((x+r2)/((P/31536000)**(2./3.)*(1.5e8*1e3)/((M2+(x/(6.5e5*1e3))**(10/7))**0.33)))**2-sympy.cos(y)**2)**0.5-
#                        (((x-r2)/((P/31536000)**(2./3.)*(1.5e8*1e3)/((M2+(x/(6.5e5*1e3))**(10/7))**0.33)))**2-sympy.cos(y)**2)**0.5)-dt12_model]
#          ,[x,y])
# for k in aa:
#     print(k)
# i_all=np.linspace(70,90,100)/180*np.pi
# err=[]
# for i in i_all:
#     j=0
#     r1_use=r1
#     a_use=a
#     M_use=M
#     M1_use=M1
#     while j <len(M1_use):
#         temp = (r1[j] - r2) / a[j] ** 2 - np.cos(i) ** 2
#         if temp<0:
#             r1_use=np.delete(r1_use,j)
#             a_use=np.delete(a_use, j)
#             M_use=np.delete(M_use, j)
#             M1_use=np.delete(M1_use, j)
#         else:
#             j+=1
#
#     dt23 = 2 * P / (2 * np.pi) * (((r1_use - r2) / a_use) ** 2 - np.cos(i) ** 2) ** 0.5
#     dt12 = P / (2 * np.pi) * ((((r1_use + r2) / a_use) ** 2 - np.cos(i) ** 2) ** 0.5 - (((r1_use - r2) / a_use) ** 2 - np.cos(i) ** 2) ** 0.5)
#
#     dt23_er=np.abs(dt23-dt23_model)
#     dt12_er = np.abs(dt12 - dt12_model)
#     if len(dt23_er)==0:
#         continue
#     else:
#         err.append(M1_use[np.where(dt23_er==np.min(dt23_er))]-M1_use[np.where(dt12_er==np.min(dt12_er))])
# print(err)
#
# plt.title('i={0}'.format(i*180/np.pi))
# dt14=2*P/(2*np.pi)*(((r1+r2)/a)**2-np.cos(i)**2)**0.5
# plt.plot(M1,np.zeros(len(M1))+dt12_model,'--')
# plt.plot(M1,np.zeros(len(M1))+dt23_model,'--')
# dt23=2*P/(2*np.pi)*(((r1-r2)/a)**2-np.cos(i)**2)**0.5
# dt12=P/(2*np.pi)*((((r1+r2)/a)**2-np.cos(i)**2)**0.5-(((r1-r2)/a)**2-np.cos(i)**2)**0.5)
# # # dt23=dt12
# plt.plot(M1,dt12)
# plt.plot(M1,dt23)
# plt.xlabel('M1')
# plt.ylabel('time')
# plt.legend(['dt12_obs','dt23_obs','dt12','dt23'])
# plt.show()


# print(dt14)
# print(dt23)
# print(dt12)