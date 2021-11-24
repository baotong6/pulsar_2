import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from astropy.io import fits
import pandas as pd
import sys
import os
from scipy import interpolate
from scipy.interpolate import interp1d,griddata
from scipy.optimize import curve_fit
import astropy.units as u
import astropy.constants as c
from scipy import interpolate
import scipy.stats
import matplotlib.font_manager as font_manager
from scipy import integrate
path='/Users/baotong/Downloads/xyz/'
os.chdir(path)
def make_3d():
    x_2d=np.array(pd.read_csv('2D/x_2d.csv'))
    y_2d=np.array(pd.read_csv('2D/y_2d.csv'))
    z_2d=np.array(pd.read_csv('2D/z_2d.csv'))
    xyz=np.zeros((50,50,3),dtype='float')
    # print(xyz)
    for i in range(len(xyz)):
        for j in range(len(xyz[0])):
            xyz[i,j]=np.array([x_2d[i,j],y_2d[i,j],z_2d[i,j]])
    np.save(path+'xyz_3d.ndy',xyz)

x_1d=pd.read_csv('1D/x_1d.csv').iloc[:,0]
y_1d=pd.read_csv('1D/y_1d.csv').iloc[:,0]
z_1d=pd.read_csv('1D/z_1d.csv').iloc[:,0]
zlog_1d=np.log10(-z_1d+1)
print(z_1d.max(),z_1d.min())

xnew = np.linspace(x_1d.min(),x_1d.max(),100)
ynew = np.linspace(y_1d.min(),y_1d.max(),100)
X,Y=np.meshgrid(xnew,ynew)
points_grid=np.zeros((100,100,2),dtype='float')
points=[]
for i in range(len(x_1d)):
    points.append([x_1d[i],y_1d[i]])
for i in range(len(xnew)):
    for j in range(len(ynew)):
        points_grid[i,j]=[xnew[i],ynew[j]]
print(points_grid)
points=np.array(points)
points_grid=np.array(points_grid)
inter_logz=griddata(points, zlog_1d, points_grid, method='nearest')
print(inter_logz)


# newfunc = interpolate.interp2d(x_1d, y_1d, zlog_1d, kind='quintic')
fig1=plt.figure(1)
ax1=Axes3D(fig1)
ax1.scatter3D(x_1d,y_1d,zlog_1d,cmap=plt.get_cmap('plasma'))

x_2d,y_2d=np.meshgrid(xnew,ynew)


fig2=plt.figure(2)
ax2=Axes3D(fig2)
ax2.scatter3D(y_2d,x_2d,inter_logz,cmap=plt.get_cmap('plasma'))
ax2.plot_surface(y_2d,x_2d,inter_logz,cmap='rainbow')
ax2.contourf3D(y_2d,x_2d,inter_logz,cmap='rainbow')

plt.show()

